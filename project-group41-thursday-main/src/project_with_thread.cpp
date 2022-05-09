









#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "RF_front_end.cpp"
#include "Mono.cpp"
#include "mono_stereo.cpp"
#include <thread>
#include <iostream>
#include <mutex>
#include <queue>
#include <condition_variable>


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


#define QUEUE_BLOCKS 5
void ms_thread(int &mode, std::queue<void *> &my_queue, std::mutex &mutex_lock, std::condition_variable& my_cvar){
  //int mode=0;
	int rf_decim;
	int audio_decim;
	float rf_Fc = 1e5;
	float rf_Fs;
	int num_taps = 501;
	int audio_Fs;
	int audio_Fc = 16e3;
	int mono_us;

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
  unsigned int block_size = 1024 * rf_decim * audio_decim * 2;

	unsigned int block_count = 1;


	std::vector<float> audio_state,state_mono_allpass,stereo_carrier_state,extraction_state,stereo_state;

	std::vector<float> audio_block_data_beforeallpass;
	std::vector<float> audio_block_data;
	std::vector<short int> audio_data;
	std::vector<float> carrier_recovery_data,stereo_mix,stereo_filt;

	std::vector<float> my_coeff,mono_coeff,stereo_carrier_coeff,extraction_coff,mixer_coeff;

	std::vector<float> stereo_recovery,extraction_data,data_comb;

	audio_state.resize(num_taps-1,0.0);
	stereo_carrier_state.resize(num_taps-1,0.0);
	extraction_state.resize(num_taps-1,0.0);
	stereo_state.resize(num_taps-1,0.0);

	float state_integrator = 0.0,state_phaseEst = 0.0,state_feedbackI = 1.0,state_feedbackQ = 0.0,state_ncoLast = 1.0;
	int state_trigOffset = 0;

  impulseResponseBPF(audio_Fs, 18500, 19500, num_taps, stereo_carrier_coeff);
	impulseResponseBPF(audio_Fs, 22000, 54000, num_taps, extraction_coff);
	impulseResponseLPF(audio_Fs, 16e3, num_taps, mixer_coeff);
	impulseResponseLPF(audio_Fs, audio_Fc, num_taps, mono_coeff);
  while(true){

    ///////////////////
    std::unique_lock<std::mutex> my_lock(mutex_lock);
    while (my_queue.empty()){
      my_cvar.wait(my_lock);
    }

    int *ptr=(int *)my_queue.front();
    my_queue.pop();
    my_cvar.notify_one();
    my_lock.unlock();


		if(mode == 0 || mode == 1){
			//mono
			filter_ds_ptr(mono_coeff, ptr, audio_state, audio_block_data_beforeallpass, audio_decim,block_size);

			allpass(audio_block_data_beforeallpass ,state_mono_allpass, audio_block_data);

			//stereo carrier recovery

			//block_processing_with_state(stereo_carrier_coeff, fm_demod, stereo_carrier_state, carrier_recovery_data);
			filter_ds_ptr(stereo_carrier_coeff, ptr, stereo_carrier_state, carrier_recovery_data,1,block_size);
			pll(carrier_recovery_data, 19e3, audio_Fs, 2.0, 0.0, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoLast, state_trigOffset, 0.01, stereo_recovery);

			//extraction

			//block_processing_with_state(extraction_coff, fm_demod, extraction_state, extraction_data);
			filter_ds_ptr(extraction_coff, ptr, extraction_state, extraction_data,1,block_size);

			//mixer
			stereo_mix.resize(stereo_recovery.size());
			for(unsigned int i = 0; i < stereo_recovery.size(); i++){
		    stereo_mix[i] = stereo_recovery[i] * extraction_data[i] * 2;
		  }

		  filter_ds(mixer_coeff, stereo_mix, stereo_state, stereo_filt, audio_decim);

			//combiner
			data_comb.clear();
			data_comb.resize(audio_block_data.size()*2,0.0);

			for(unsigned int i = 0; i < audio_block_data.size(); i++){
				data_comb[2*i] = (audio_block_data[i] + stereo_filt[i]);
				data_comb[2*i + 1] = (audio_block_data[i] - stereo_filt[i]);
			}


		} else{//mode 2 & 3
			//mono
			//filter_ds_with_us(mono_coeff, fm_demod, audio_state, audio_block_data_beforeallpass, audio_decim, mono_us);

			//allpass(audio_block_data_beforeallpass ,state_mono_allpass, audio_block_data);

			//stereo carrier recovery

			//block_processing_with_state(stereo_carrier_coeff, fm_demod, stereo_carrier_state, carrier_recovery_data);
		}


		audio_data.resize(data_comb.size());
		//audio_data.resize(data_comb.size());
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

		if((std::cin.rdstate()) != 0&& my_queue.empty()){
			break;
		}
    else {
			block_count ++;
		}
 	}





}




void rf_thread(int &mode, std::queue<void *> &my_queue, std::mutex &mutex_lock, std::condition_variable& my_cvar){
	int rf_decim;
	int audio_decim;
	float rf_Fc = 1e5;
	float rf_Fs;
	int num_taps = 501;
	int audio_Fs;
	int audio_Fc = 16e3;
	int mono_us;


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

	unsigned int block_size = 1024 * rf_decim * audio_decim * 2;

	std::vector<float> iq_data,i_data,q_data,i_filter,q_filter,fm_demod;
	std::vector<float> i_state;
	std::vector<float> q_state;

	std::vector<float> my_coeff;

	i_state.resize(num_taps-1,0.0);
	q_state.resize(num_taps-1,0.0);
	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);
	float last_i = 0.0;
	float last_q = 0.0;

	unsigned int block_count = 1;

	static float queue_block[QUEUE_BLOCKS][block_size];
	while(true){
		readStdInBlock(block_size, block_count, iq_data);

		for(unsigned int i = 0; i < (block_size)/2; i++){
			i_data[i] = iq_data[i*2];
			q_data[i] = iq_data[1+i*2];
		}
		impulseResponseLPF(rf_Fs, rf_Fc, num_taps, my_coeff);

		unsigned int queue_entry = block_count % QUEUE_BLOCKS;
		//Demoadulate data
		float * pointer_block = &queue_block[queue_entry][0];
		RF_front_end(my_coeff, i_data, q_data, i_state, q_state, block_size, rf_Fc, rf_Fs, rf_decim, num_taps, i_filter, q_filter);

		demodulator(pointer_block, i_filter, q_filter, last_i, last_q);
		std::unique_lock<std::mutex> my_lock(mutex_lock);
		if(my_queue.size() == QUEUE_BLOCKS-1)
			{
				my_cvar.wait(my_lock);
			}
			my_queue.push((void *)&queue_block[queue_entry][0]);

			//Fills with zeros
			std::fill(i_filter.begin(), i_filter.end(), 0);
			std::fill(q_filter.begin(), q_filter.end(), 0);

			//iterate block id
			block_count ++;

			//Unlock and notify
			my_lock.unlock();
			my_cvar.notify_one();
			if((std::cin.rdstate()) != 0)
			{
				break;
			}

	}

}

int main(){
  int mode = modeSelection(argc, argv);
  std::queue<void *> my_queue;
  std::mutex mutex_lock;
  std::condition_variable my_cvar;


  std::thread rf=std::thread(rf_thread,std::ref(mode),std::ref(my_queue),std::ref(mutex_lock),std::ref(my_cvar));
  std::thread ms=std::thread(ms_thread,std::ref(mode),std::ref(my_queue),std::ref(mutex_lock),std::ref(my_cvar));

  rf.join();
	ms.join();



  return 0;



}
