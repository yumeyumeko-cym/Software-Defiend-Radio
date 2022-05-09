#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"


void downsampler(std::vector<float> &data, unsigned int factor){
  for(unsigned int i = 0; i < (data.size()/factor);i++){
    data[i] = data[i*factor];
  }
  data.resize(int(data.size()/factor));
  data.shrink_to_fit();
}

void fmDemodArctan(std::vector<float> &fm_demod, std::vector<float> &I, std::vector<float> &Q, std::vector<float> &prev_phase)
{
    float dQ;
    float dI;
    float scaling;
    float phase;
    // prev_phase = {0, 0};
    fm_demod.resize(I.size());
    for(unsigned int i =0; i<I.size(); i++){
        dQ = Q[i] - prev_phase[1];
        dI = I[i] - prev_phase[0];
        if(pow(I[i],2) + pow(Q[i],2)==0){
            phase = 0;
            // fm_demod.push_back(0);
            fm_demod[i] = 0;
        } else {
            scaling = 1/(pow(I[i], (float)2) + pow(Q[i], (float)2));
            phase = (I[i]*dQ - Q[i]*dI)*scaling;
            fm_demod[i]=phase;
        }

        prev_phase = {I[i], Q[i]};
    }
}


void RF_front_end(std::vector<float> &my_coeff, std::vector<float> &i_data, std::vector<float> &q_data, std::vector<float> &i_state, std::vector<float> &q_state, unsigned int block_size, float Fc, float Fs, unsigned int factor, int num_taps, std::vector<float> &data_filtered_i_block, std::vector<float> &data_filtered_q_block){
	//std::vector<float> my_coeff;

	//impulseResponseLPF(Fs, Fc, num_taps, my_coeff);
	//std::vector<float> i_data_filtered;
  //std::vector<float> q_data_filtered;
	//std::vector<float>	i_state;
	//std::vector<float>	q_state;
	//std::vector<float> data_i_block;
	//std::vector<float> data_q_block;

  //i_data_filtered.clear();
  //i_data_filtered.resize(i_data.size(),0.0);

  //q_data_filtered.clear();
  //q_data_filtered.resize(q_data.size(),0.0);

	//i_state.clear();
	//i_state.resize(num_taps-1,0.0);

	//q_state.clear();
	//q_state.resize(num_taps-1,0.0);

	//data_i_block.clear();
	//data_i_block.resize(block_size,0.0);

	//data_q_block.clear();
	//data_q_block.resize(block_size,0.0);

	data_filtered_i_block.clear();
	data_filtered_i_block.resize(block_size,0.0);

	data_filtered_q_block.clear();
	data_filtered_q_block.resize(block_size,0.0);

	//unsigned int position = 0;
  //unsigned int index_count = 0;

	//while(true){
		//data_i_block.assign((i_data.begin()+position),(i_data.begin()+position+block_size));
		//data_q_block.assign((q_data.begin()+position),(q_data.begin()+position+block_size));

	filter_ds(my_coeff, i_data, i_state, data_filtered_i_block,factor);
	filter_ds(my_coeff, q_data, q_state, data_filtered_q_block,factor);

		//for(unsigned int i=0; i < data_filtered_i_block.size(); i++){
	//		i_data_filtered[index_count+i] = data_filtered_i_block[i];
	//		q_data_filtered[index_count+i] = data_filtered_q_block[i];
	//	}

    //index_count += data_filtered_i_block.size();

		//position += block_size;

		//if(position > i_data.size()){
		//	break;
		//}
	//}
  //i_data_filtered.resize(i_data.size()/factor);
//  q_data_filtered.resize(q_data.size()/factor);
}
