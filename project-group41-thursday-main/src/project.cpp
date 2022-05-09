#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "RF_front_end.h"
#include <thread>
#include <iostream>
#include "mono_stereo.h"
#include <mutex>
#include <queue>
#include <condition_variable>
int generated=0;

#define pi 3.14159265358979323846
void demodulator(std::vector<float> &fm_demod, std::vector<float> &i_data, std::vector<float> &q_data, float &last_i, float &last_q){
	fm_demod.clear();
	fm_demod.resize(i_data.size(),0.0);

	float der_i;
	float der_q;

	for(unsigned int k = 0; k < i_data.size(); k++){
		der_q = q_data[k] - last_q;
		der_i = i_data[k] - last_i;

		last_q = q_data[k];
		last_i = i_data[k];

		fm_demod[k] = (i_data[k] * der_q - q_data[k] * der_i)/(pow(i_data[k],2) + pow(q_data[k],2));
	}
}

int modeSelection(int argc, char *argv[]){
	int mode = 0;
	if (argc < 2){
		std::cerr << "no mode selected, will run in mode 0 by default..." << std::ends;
	}
	else if (argc == 2){
		mode = atoi(argv[1]);
		if (mode > 3){
			std::cerr << "invalid mode number detected, terminating..." << std::ends;
			exit(1);
		}
	}
	else{
		std::cerr << "please read document to properly select mode" << std::ends;
	}
	return mode;
}



/*void my_compute(std::vector<float>& elem){
  float hash=0.0;
  for (int k=0;k<(int)elem.size();k++){
  hash=elem[k];
  //std::cout<<"Received "<<hash<<std::endl;
  }
  //std::cerr<<"Received"<<hash<<std::endl;


}
*/

void ms_thread(int &mode, std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar){
  int rf_decim;
	int audio_decim;
	float i_Fs;
	int audio_Fc = 16e3;
	int mono_us;
	int counter=0;
	int num_taps = 101;

	std::vector<float> audio_state;
	std::vector<float> state_mono_allpass;
	std::vector<float> stereo_carrier_state;
	std::vector<float> extraction_state;
	std::vector<float> stereo_state;

	std::vector<float> audio_block_data_allpass;
	std::vector<float> audio_block_data;
	std::vector<short int> audio_data;
	std::vector<float> carrier_recovery_data;
	std::vector<float> stereo_mix;
	std::vector<float> stereo_filt;

	std::vector<float> my_coeff;
	std::vector<float> mono_coeff;
	std::vector<float> stereo_carrier_coeff;
	std::vector<float> extraction_coff;
	std::vector<float> mono_us_coeff;

	std::vector<float> stereo_recovery;
	std::vector<float> extraction_data,data_comb;

	if(mode == 0){
		i_Fs = 2.4e5;
		rf_decim = 10;
		audio_decim = 5;
		mono_us = 1;
		impulseResponseLPF(i_Fs, audio_Fc, num_taps, mono_coeff);
	} else if(mode == 1){
		i_Fs = 2.4e5;
		rf_decim = 8;
		audio_decim = 5;
		mono_us = 1;
		impulseResponseLPF(i_Fs, audio_Fc, num_taps, mono_coeff);
	}else if(mode == 2){
		i_Fs = 2.4e5;
		rf_decim = 10;
		audio_decim = 800;
		mono_us = 147;
		impulseResponseLPF(i_Fs*mono_us, audio_Fc, num_taps*mono_us, mono_coeff);
	}else{
		i_Fs = 2.88e5;
		rf_decim = 4;
		audio_decim = 320;
		mono_us = 49;
		impulseResponseLPF(i_Fs*mono_us, audio_Fc, num_taps*mono_us, mono_coeff);
	}
  unsigned int block_size = 1024 * rf_decim  * 2;

	unsigned int block_count = 1;




	audio_state.resize(num_taps-1,0.0);
	stereo_carrier_state.resize(num_taps-1,0.0);
	extraction_state.resize(num_taps-1,0.0);
	stereo_state.resize(num_taps-1,0.0);
	state_mono_allpass.resize(block_size/2,0.0);

	float state_integrator = 0.0,state_phaseEst = 0.0,state_feedbackI = 1.0,state_feedbackQ = 0.0,state_ncoLast = 1.0;
	int state_trigOffset = 0;

  impulseResponseBPF(i_Fs, 18500, 19500, num_taps, stereo_carrier_coeff);
	impulseResponseBPF(i_Fs, 22000, 54000, num_taps, extraction_coff);

  int k=0;
  while(1){

    //////////////////////////////////////////////
    std::unique_lock<std::mutex> my_lock(my_mutex);
    while(my_queue.empty()){
      my_cvar.wait(my_lock);
    }
    std::vector<float> fm_demod=my_queue.front();
    my_queue.pop();
    my_cvar.notify_one();
    my_lock.unlock();
    k++;
    //std::cout<<"extracted block  "<<k<<std::endl;
    //my_compute(elem);
    //////////////////////////////////////////////

		if(std::isnan(fm_demod[0])) {
			fm_demod[0] = 0;
		}



    if(mode == 0 || mode == 1){
			allpass(fm_demod ,state_mono_allpass, audio_block_data_allpass);
			filter_ds(mono_coeff, audio_block_data_allpass, audio_state, audio_block_data, audio_decim);

			//stereo carrier recovery

			//block_processing_with_state(stereo_carrier_coeff, fm_demod, stereo_carrier_state, carrier_recovery_data);
			filter_ds(stereo_carrier_coeff, fm_demod, stereo_carrier_state, carrier_recovery_data,1);
			pll(carrier_recovery_data, 19e3, i_Fs, 2.0, 0.0, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoLast, state_trigOffset, 0.01, stereo_recovery);


			//extraction
			//block_processing_with_state(extraction_coff, fm_demod, extraction_state, extraction_data);
			filter_ds(extraction_coff, fm_demod, extraction_state, extraction_data,1);

			//mixer
			stereo_mix.resize(stereo_recovery.size());
			for(unsigned int i = 0; i < stereo_recovery.size(); i++){
		    stereo_mix[i] = stereo_recovery[i] * extraction_data[i] * 2;
		  }



		  filter_ds(mono_coeff, stereo_mix, stereo_state, stereo_filt, audio_decim);

			//combiner
			data_comb.clear();
			data_comb.resize(audio_block_data.size()*2,0.0);

			for(unsigned int i = 0; i < audio_block_data.size(); i++){
				data_comb[2*i] = (audio_block_data[i] + stereo_filt[i]);
				data_comb[2*i + 1] = (audio_block_data[i] - stereo_filt[i]);
			}


		} else{//mode 2 & 3
			//mono
			allpass(fm_demod ,state_mono_allpass, audio_block_data_allpass);
			filter_ds_with_us(mono_coeff, audio_block_data_allpass, audio_state, audio_block_data, audio_decim, mono_us);

			//stereo carrier recovery

			//block_processing_with_state(stereo_carrier_coeff, fm_demod, stereo_carrier_state, carrier_recovery_data);
			filter_ds_with_us(stereo_carrier_coeff, fm_demod, stereo_carrier_state, carrier_recovery_data,1,1);
			pll(carrier_recovery_data, 19e3, i_Fs, 2.0, 0.0, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoLast, state_trigOffset, 0.01, stereo_recovery);


			//extraction
			//impulseResponseBPF(rf_Fs, 22000, 54000, num_taps, extraction_coff);
			//block_processing_with_state(extraction_coff, fm_demod, extraction_state, extraction_data);
			filter_ds_with_us(extraction_coff, fm_demod, extraction_state, extraction_data,1,1);

			//mixer
			stereo_mix.resize(stereo_recovery.size());
			for(unsigned int i = 0; i < stereo_recovery.size(); i++){
		    stereo_mix[i] = stereo_recovery[i] * extraction_data[i] * 2;
		  }



		  filter_ds_with_us(mono_coeff, stereo_mix, stereo_state, stereo_filt, audio_decim, mono_us);
			/*
			std::cerr << "test" << '\n';
			std::cerr << fm_demod.size() << '\n';
			std::cerr << audio_block_data.size() << '\n';
			std::cerr << stereo_filt.size() << '\n';
			*/
			//combiner
			data_comb.clear();
			data_comb.resize(audio_block_data.size()*2,0.0);

			for(unsigned int i = 0; i < audio_block_data.size(); i++){
				data_comb[2*i] = (audio_block_data[i] + stereo_filt[i]);
				data_comb[2*i + 1] = (audio_block_data[i] - stereo_filt[i]);
			}

		}
		counter++;
		/*
		std::cerr << "test" << '\n';
		std::cerr << carrier_recovery_data[2] << '\n';
		std::cerr << carrier_recovery_data[3] << '\n';
		std::cerr << audio_block_data_allpass[0] << '\n';
		std::cerr << audio_block_data_allpass[1] << '\n';
		*/
		/*std::cerr << "abc" << '\n';
		std::cerr << mono_coeff[0] << '\n';
		std::cerr << audio_block_data[1] << '\n';
		std::cerr << stereo_filt[2] << '\n';
		std::cerr << stereo_filt[3] << '\n';
		*/


		audio_data.resize(data_comb.size());
		 //std::cerr << "/* bbb */" << '\n';
		for(unsigned int i=0; i<data_comb.size(); i++){
			//If block is null set as zero
			if(std::isnan(data_comb[i])) {
				audio_data[i] =0;
			} else { //Cast the value to short int
				audio_data[i] = static_cast<short int>(data_comb[i] *16384);
			}

		}




		//Write to standard out

		fwrite(&audio_data[0],sizeof(short int),audio_data.size(), stdout);

		block_count++;
    if(my_queue.empty() && (std::cin.rdstate()) != 0){
			//std::cerr << "ccc" << '\n';
      break;
    }




  }

}


void rf_thread(int &mode,std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar){
  int rf_decim;
	float rf_Fc = 1e5;
	float rf_Fs;
	int num_taps = 501;
	//int audio_Fs;
	//int audio_Fc = 16e3;
	//int mono_us;

	if(mode == 0){
		rf_Fs = 2.4e6;
		rf_decim = 10;
		//audio_Fs = 48e3;
		//mono_us = 1;
	} else if(mode == 1){
		rf_Fs = 1.92e6;
		rf_decim = 8;
		//audio_Fs = 48e3;
	//mono_us = 1;
	}else if(mode == 2){
		rf_Fs = 2.4e6;
		rf_decim = 10;
		//audio_Fs = 441e2;
		//mono_us = 147;
	}else{
		rf_Fs = 1.152e6;
		rf_decim = 4;
		//audio_Fs = 441e2;
		//mono_us = 49;
	}




	unsigned int block_size = 1024 * rf_decim  * 2;

	std::vector<float> iq_data,i_data,q_data,i_filter,q_filter;
	std::vector<float> i_state;
	std::vector<float> q_state;

	std::vector<float> my_coeff;

	std::vector<float> phase_IQ;


	i_state.resize(num_taps-1,0.0);
	q_state.resize(num_taps-1,0.0);
	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);
	i_filter.clear();
	i_filter.resize(block_size,0.0);
	q_filter.clear();
	q_filter.resize(block_size,0.0);

	float last_i = 0.0;
	float last_q = 0.0;

	unsigned int block_count = 1;

	impulseResponseLPF(rf_Fs, rf_Fc, num_taps, my_coeff);

  while(1){

    //generated++;

    readStdInBlock(block_size, block_count, iq_data);

		for(unsigned int i = 0; i < (block_size)/2; i++){
			i_data[i] = iq_data[i*2];
			q_data[i] = iq_data[1+i*2];
		}


		//unsigned int queue_entry = block_count % QUEUE_BLOCKS;
		//Demoadulate data

		std::vector<float> queue_block;

		//filter_ds_rf(my_coeff, i_data, i_state, i_filter, q_data, q_state, q_filter,rf_decim);
		RF_front_end(my_coeff, i_data, q_data, i_state, q_state, block_size, rf_Fc, rf_Fs, rf_decim, num_taps, i_filter, q_filter);

		//demodulator(queue_block, i_filter, q_filter, last_i, last_q);
		fmDemodArctan(queue_block, i_filter, q_filter, phase_IQ);

		/*std::cerr << queue_block[0] << '\n';
		std::cerr << queue_block[1] << '\n';
		std::cerr << queue_block[2] << '\n';
		std::cerr << queue_block[3] << '\n';

		std::cerr << queue_block.size() << '\n';
		std::cerr << "fm_demod " << '\n';
		*/
    ///////////////////////////////////////////////////////////
    std::unique_lock<std::mutex> my_lock(my_mutex);
    while(my_queue.size()>=4){

      my_cvar.wait(my_lock);
    }
		//std::cerr << "queue_block[0]" << '\n';
		//std::cerr << queue_block[0] << '\n';

    my_queue.push(queue_block);
    my_cvar.notify_one();
    my_lock.unlock();
		//std::cerr << queue_block[0] << '\n';

  /////////////////////////
  block_count++;
  if( (std::cin.rdstate()) != 0)
			{
				//std::cerr << "bbb" << '\n';
				break;
			}
    }
		//std::cerr << "test" << '\n';

}

void rrc(float Fs, float num_taps, std::vector<float> &impulseResponseRRC){
  float T_symbol = 1/2375.0;
  float beta = 0.90;

  impulseResponseRRC.clear();
  impulseResponseRRC.resize(num_taps);

  for(unsigned int k = 0; k < num_taps; k++){
    float t = (float)(k-num_taps/2)/Fs;
    if (t == 0.0){
      impulseResponseRRC[k] = 1.0 + beta*(4/pi-1);
    } else if (t == -T_symbol/(4*beta) || t == T_symbol/(4*beta)){
      impulseResponseRRC[k] = (beta/sqrt(2))*(((1+2/pi)* \
					(sin(pi/(4*beta)))) + ((1-2/pi)*(cos(pi/(4*beta)))));
    } else {
      impulseResponseRRC[k] = (sin(pi*t*(1-beta)/T_symbol) +  \
					4*beta*(t/T_symbol)*cos(pi*t*(1+beta)/T_symbol))/ \
					(pi*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);
    }
  }
}



void rds_thread(){

	int rf_decim;
	int audio_decim;
	float rf_Fc = 1e5;
	float rf_Fs;
	int num_taps = 501;
	int audio_Fs;
	int audio_Fc = 16e3;
	int mono_us;
	int mode = 0;

	int rational_us = 209;
	int rational_ds = 960;

	if(mode == 0){
		rf_Fs = 2.4e6;
		rf_decim = 10;
		audio_decim = 5;
		audio_Fs = 48e3;
		mono_us = 1;
	} else if(mode == 1){
		rf_Fs = 1.92e6;
		rf_decim = 8;
		audio_decim = 5;
		audio_Fs = 48e3;
		mono_us = 1;
	}else if(mode == 2){
		rf_Fs = 2.4e6;
		rf_decim = 10;
		audio_decim = 800;
		audio_Fs = 441e2;
		mono_us = 147;
	}else if(mode == 3){
		rf_Fs = 1.152e3;
		rf_decim = 4;
		audio_decim = 320;
		audio_Fs = 441e2;
		mono_us = 49;
	}

	std::vector<float> rds_extraction_coeff;
	std::vector<float> rds_carrier_coeff;
	std::vector<float> rds_demod_lpf_coeff;
	std::vector<float> rrc_coeff;

	std::vector<float> rds_extraction_state;
	std::vector<float> rds_carrier_state;
	std::vector<float> rds_extraction_data_allpass_state;
	std::vector<float> rds_demod_lpf_state;
	std::vector<float> rrc_state;

	std::vector<float> rds_extraction_data;
	std::vector<float> rds_extraction_data_sn;
	std::vector<float> rds_carrier_data;
	std::vector<float> rds_carrier_data_pll;
	std::vector<float> rds_carrier_data_pll_Q;
	std::vector<float> rds_extraction_data_allpass;
	std::vector<float> rds_data_mix;
	std::vector<float> rds_demod_lpf_data;
	std::vector<float> rds_rational_resample_data;
	std::vector<float> rrc_data;

  std::vector<bool> info_decode;
  std::vector<bool> diff_decode;


	std::vector<float> fm_demod;

	std::vector<float> symbol;
	unsigned int block_count = 0, temp_high, offset;
	bool prev_bit;

	impulseResponseBPF(rf_Fs, 54.0e3, 60.0e3, num_taps,rds_extraction_coeff);
	impulseResponseBPF(rf_Fs, 113.5e3, 114.5e3, num_taps,rds_carrier_coeff);
	impulseResponseLPF(rf_Fs, 3e3, num_taps, rds_demod_lpf_coeff);

	float state_integrator = 0.0,state_phaseEst = 0.0,state_feedbackI = 1.0,state_feedbackQ = 0.0,state_ncoLast = 1.0,state_ncoQLast = 1.0;
	int state_trigOffset = 0;
	unsigned int count_01 = 0, count_12 = 0, start;



	while(true){
		//RDS channel extraction
		block_processing_with_state(rds_extraction_coeff, fm_demod, rds_extraction_state, rds_extraction_data);

		//RDS carrier recovery
		rds_extraction_data_sn.resize(rds_extraction_data.size(),0.0);
		//squaring nonlinearity
		for(unsigned int i = 0;i<rds_extraction_data.size();i++){
			rds_extraction_data_sn[i] = rds_extraction_data[i] * rds_extraction_data[i];
		}
		//BPF
		block_processing_with_state(rds_carrier_coeff, rds_extraction_data_sn, rds_carrier_state, rds_carrier_data);

		//pll
		rds_carrier_data_pll.resize(rds_carrier_data.size(),0.0);
		pllIQ(rds_carrier_data, 114e3, audio_Fs, 0.5, 0.0, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoLast, state_ncoQLast, state_trigOffset, 0.01, rds_carrier_data_pll, rds_carrier_data_pll_Q);

		//APF
		allpass(rds_extraction_data, rds_extraction_data_allpass_state, rds_extraction_data_allpass);

		//miixng
		rds_data_mix.resize(rds_carrier_data_pll.size(),0.0);
		for (unsigned int i = 0; i < rds_carrier_data_pll.size();i++){
			rds_data_mix[i] = rds_carrier_data_pll_Q[i] * rds_carrier_data_pll[i] * rds_extraction_data_allpass[i] * 2;
		}

		//LPF
		block_processing_with_state(rds_demod_lpf_coeff, rds_data_mix, rds_demod_lpf_state, rds_demod_lpf_data);

		//rational resampler
		rds_rational_resample_data.clear();
		rds_rational_resample_data.resize(rds_demod_lpf_data.size() * rational_us,0.0);
		for(unsigned int i=0; i < rds_demod_lpf_data.size(); i++){
			rds_rational_resample_data[i*rational_us] = rds_demod_lpf_data[i];
		}
		for(unsigned int i=0; i < rds_rational_resample_data.size()/rational_ds; i++){
			rds_rational_resample_data[i] = rds_rational_resample_data[i * rational_ds];
		}
		rds_rational_resample_data.resize(rds_rational_resample_data.size()/rational_ds);

		//rrc
		rrc(22*2375, num_taps, rrc_coeff);
		rrc_data.resize(rds_rational_resample_data.size(),0.0);
		block_processing_with_state(rrc_coeff, rds_rational_resample_data, rrc_state, rrc_data);

		symbol.clear(); symbol.resize(rrc_data.size()/22, 0.0);
		//clock and data recovery
		if (block_count == 0){
			temp_high = rrc_data[0];
			offset = 0;
			for (unsigned int i = 0; i < 22; i++){ // 22, constraint for multiples
				if(rrc_data[i] > temp_high){
					temp_high = rrc_data[i];
					offset = i;
				}
			}
    }
		for(unsigned int k = 0; k < symbol.size(); k++){
			symbol[k] = rrc_data[22*k+offset];
		}

		for(unsigned int m; m < symbol.size()/4; m++){
			if ((symbol[2*m] > 0 && symbol[2*m+1] > 0 )|| (symbol[2*m] < 0 && symbol[2*m+1] < 0)){ // compare if both of the signs are High or Low
				count_01 += 1;
			} else if ((symbol[2*m+1] > 0 && symbol[2*m+2] > 0) || (symbol[2*m+1] < 0 && symbol[2*m+2]<0)) {
				count_12 += 1;
			}
		}
		if (count_01 > count_12){
			start = 1;
		} else {
			start = 0;
		}

    info_decode.resize((int)(symbol.size()/2)-start,0);

    // Decode the symbol
    for(unsigned int i = 0; i < info_decode.size();i++){
      // break when reaching the array boundary
      if (2*i + 1 + start > symbol.size()-1){
        break;
      }
      if(symbol[2*i+1] > symbol[2*i+1+start]){
        info_decode[i] = true; // HL
      } else {
        info_decode[i] = false; // LH
      }
    }

    diff_decode.resize(info_decode.size(),0.0);

    // Differential Decoding
    // XOR for differential Decoding
    if(block_count == 0){
      prev_bit = info_decode[0];
      diff_decode[0] = info_decode[0];
      for(unsigned int i = 1; i < info_decode.size(); i++){
        diff_decode[i] = prev_bit ^ info_decode[i];
        prev_bit = info_decode[i];
      }
    } else{
      for(unsigned int i = 0; i < info_decode.size(); i++){
        diff_decode[i] = prev_bit ^ info_decode[i];
        prev_bit = info_decode[i];
      }
    }

	}
}

int main(int argc, char* argv[]){
  int mode = modeSelection(argc, argv);
  std::queue<std::vector<float>> my_queue;
  std::mutex my_mutex;
  std::condition_variable my_cvar;

  std::thread rf=std::thread(ms_thread,std::ref(mode),std::ref(my_queue),std::ref(my_mutex),std::ref(my_cvar));
  std::thread ms=std::thread(rf_thread,std::ref(mode),std::ref(my_queue),std::ref(my_mutex),std::ref(my_cvar));

  rf.join();
  ms.join();


  return 0;

}
