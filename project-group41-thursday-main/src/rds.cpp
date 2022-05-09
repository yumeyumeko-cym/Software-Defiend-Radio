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

  std::vector<float> info_decode;
  std::vector<float> diff_decode;


	std::vector<float> fm_demod;

	std::vector<float> symbol;
	unsigned int block_count = 0, temp_high, offset;

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

    info_decode.resize((int)(symbol.size()/2)-start_pos,0);

    // Decode the symbol
    for(unsigned int i = 0; i < info_decode.size();i++){
      // break when reaching the array boundary
      if (2*i + 1 + start > symbol.size()-1){
        break;
      }
      if(symbol[2*i+1] > symbol[2*i+1+start]){
        info_decode[i] = 1; // HL
      } else {
        info_decode[i] = 0; // LH
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
