#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"


void monoPart(std::vector<float> audio_coeff, float Fc,float Fs, int N_taps,std::vector<float> &fm_demod, std::vector<float> &fm_demod_state, int mode_num,std::vector<float> &block_data, int us, int ds){
	//std::vector<float> audio_coeff;
	//std::vector<float> after_us;
	//std::vector<float> block_data;
	//std::vector<float> fm_demod_state;
	//std::vector<float> fm_block_data;

	//block_data.clear();
	//after_us.clear();
	//fm_demod_state.clear();
	//fm_block_data.clear();

	//fm_demod_state.resize(N_taps-1,0.0);
	//fm_block_data.resize(block_size,0.0);

	//fm_block_data.assign((fm_demod.begin()+position),(fm_demod.begin()+position+block_size));

	if(mode_num==1 || mode_num==0){
		//impulseResponseLPF(Fs,Fc, N_taps, audio_coeff);
		filter_ds(audio_coeff,fm_demod,fm_demod_state,block_data,ds);
	}

	// us is the sp sampler factor that to be determined. just put the integer value to replace
  else{
		//impulseResponseLPF(Fs*us,Fc, N_taps*us, audio_coeff);
			//upsampler(fm_block_data,us,after_us);
		filter_ds_with_us(audio_coeff,fm_demod,fm_demod_state,block_data,ds,us);
	}

}
