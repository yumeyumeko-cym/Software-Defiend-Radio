/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
  h.clear(); h.resize(num_taps, 0.0);
	float norm_Cutoff = Fc/(Fs/2);

	for(int i=0; i<num_taps; i++){
		if(i == (num_taps-1)/2){
			h[i] = norm_Cutoff;
		} else{
			h[i] = norm_Cutoff * std::sin(PI*norm_Cutoff*(i-(num_taps-1)/2))/(PI*norm_Cutoff*(i-(num_taps-1)/2));
		}
		h[i] = h[i] * pow(std::sin(i*PI/num_taps),2);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	for (unsigned int n=0;n<y.size();n++){
  	for (unsigned int k=0;k<h.size();k++){
   		if ((n-k)>=0 and n-k < x.size()){
    		y[n]=y[n]+h[k]*x[n-k];
			}
		}
	}
}


void block_processing_with_state(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb){
	yb.clear();
	yb.resize(xb.size(), 0.0);
  int count;

	for(unsigned int n = 0;n<xb.size();n++){
    count = 0;
		for(unsigned int k = 0;k<h.size();k++){
			if(n - k >= 0 && (n-k)<xb.size()){
				yb[n] += h[k] * xb[n-k];
			}
			else{
				yb[n] += h[k] * state[state.size()-count-1];
			}
		}
	}
	//state = xb[(xb.size()-state.size()):(xb.size())];
	state = std::vector<float>(xb.begin()+(xb.size()-(h.size()-1)), xb.end());
}

void filter_signal(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb){
	yb.clear();
	yb.resize(xb.size(), 0.0);

	block_processing_with_state(h, xb, state, yb);
}

void lpf_block_processing(std::vector<uint8_t> &data, unsigned int block_size, float Fc, float Fs, int N_taps, std::vector<float> &filtered_data){
	std::vector<float> my_coeff;

	impulseResponseLPF(Fs, Fc, N_taps, my_coeff);

	filtered_data.clear();
	filtered_data.resize(data.size());

	//start at the first block (with relative position zero)
	unsigned int position = 0;

	//intiial filter state - state is the size of the impulse response minus 1
	//we need to channels for the state (one for each audio channel)
	std::vector<float> filtered_state;
	filtered_state.clear();
	filtered_state.resize(my_coeff.size()-1,0.0);

  std::vector<float> data_block;
	data_block.clear();
	data_block.resize(block_size,0.0);

	std::vector<float> my_data;

	my_data.clear();
	my_data.resize(block_size,0.0);

	while(true){
		data_block.assign((data.begin()+position),(data.begin()+position+block_size));

		filter_signal(my_coeff, data_block, filtered_state, my_data);

		for(unsigned int i=0; i < block_size; i++){
			filtered_data[position+i] = my_data[i];

			//std::cout<<filtered_data_left[position+i]<<"\n";
		}

		position += block_size;

		if(position > data.size()){
			break;
		}
	}


}

void filter_ds_rf(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb, const std::vector<float> &xb_1, std::vector<float> &state_1, std::vector<float> &yb_1,const int &factor){
	yb.clear();
	yb.resize(xb.size(), 0.0);
  yb_1.clear();
  yb_1.resize(xb.size(), 0.0);
  int count;

	for(unsigned int n = 0;n<xb.size()/factor;n++){
    count = 0;
		for(unsigned int k = 0;k<h.size();k++){
			if(factor*n - k >= 0 && (factor*n-k < xb.size())){
				yb[n] += h[k] * xb[factor*n-k];
        yb_1[n] += h[k] * xb_1[factor*n-k];

			}
			else{
				yb[n] += h[k] * state[state.size()-count - 1];
        yb_1[n] += h[k] * state_1[state_1.size()-count - 1];
        count += 1;
			}
		}
	}
  yb.resize(yb.size()/factor);
  yb.shrink_to_fit();
  yb_1.resize(yb.size()/factor);
  yb_1.shrink_to_fit();
	//state = xb[(xb.size()-state.size()):(xb.size())];
	//state.assign(xb.begin()+(xb.size()-state.size()), xb.begin()+xb.size());
  //state_1.assign(xb_1.begin()+(xb_1.size()-state_1.size()), xb_1.begin()+xb_1.size());


	state = std::vector<float>(xb.begin()+(xb.size()-(h.size()-1)), xb.end());
	state_1 = std::vector<float>(xb_1.begin()+(xb.size()-(h.size()-1)), xb_1.end());

}

void filter_ds(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb, const int &factor){
	yb.clear();
	yb.resize(xb.size(), 0.0);
  int count;

	for(unsigned int n = 0;n<xb.size()/factor;n++){
    count = 0;
		for(unsigned int k = 0;k<h.size();k++){
			if(factor*n - k >= 0 && (factor*n-k < xb.size())){
				yb[n] += h[k] * xb[factor*n-k];
			}
			else{
				yb[n] += h[k] * state[state.size()-count - 1];
        count += 1;
			}
		}
	}
  yb.resize(yb.size()/factor);
  yb.shrink_to_fit();
	//state = xb[(xb.size()-state.size()):(xb.size())];
	//state.assign(xb.begin()+(xb.size()-state.size()), xb.begin()+xb.size());

  state = std::vector<float>(xb.begin()+(xb.size()-(h.size()-1)), xb.end());

}


void filter_ds_with_us(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb, const int &factor,const int &us){
	yb.clear();
	yb.resize(xb.size()*us/factor, 0.0);
	int phase=0;
  int count;

	for(unsigned int n = 0;n<yb.size();n++){
    count = 0;
    phase=(n*factor)%us;
		for(unsigned int k = 0;k<h.size()/us;k++){
			if((n*factor-phase)/us >= 0 && (n*factor-phase)/us < xb.size()){
				yb[n] += (h[phase+k*us]*xb[(n*factor-phase)/us])*us;
			}else{
				//yb[n] += h[phase+k*us] * state[state.size()+((n*factor-phase)/us)-(phase+k*us)];
        //yb[n] += h[phase+k*us] * state[state.size()+(n-k)];
        yb[n] += h[phase+k*us] * state[state.size()-count -1];
			}
		}
	}
  //yb.resize(xb.size()*us/factor);
  //yb.shrink_to_fit();
	//state = xb[(xb.size()-state.size()):(xb.size())];
	//state.assign(xb.begin()+(xb.size()-state.size()), xb.begin()+xb.size());

  state = std::vector<float>(xb.begin()+(xb.size()-(h.size()-1)), xb.end());

}

void impulseResponseBPF(float Fs, float Fb,float Fe, unsigned short int num_taps, std::vector<float> &h){
	// bring your own functionality
  // allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);
	float norm_center = ((Fe+Fb)/2)/(Fs/2);
  float norm_pass=(Fe-Fb)*2/Fs;

	for(int i=0; i<num_taps; i++){
		if(i == (num_taps-1)/2){
			h[i] = norm_pass;
		} else{
			h[i] = norm_pass * std::sin(PI*(norm_pass/2)*(i-(num_taps-1)/2))/(PI*(norm_pass/2)*(i-(num_taps-1)/2));
		}
    h[i]=h[i]*std::cos(i*PI*norm_center);
		h[i] = h[i] * pow(std::sin(i*PI/num_taps),2);
	}

}

void allpass(const std::vector<float> &input_block, std::vector<float> &state_block, std::vector<float> &output_block){
  state_block.resize(input_block.size()/2);
  output_block.resize(input_block.size(),0.0);

  for(unsigned int i=0; i < state_block.size(); i++){
    output_block[i] = state_block[i];
    output_block[i + state_block.size()] = input_block[i];
    state_block[i] = input_block[i + state_block.size()];
  }

  //output_block.assign(state_block.begin(),input_block.begin()-state_block.size());
  //state_block.assign(input_block.begin()+state_block.size(),input_block.begin()+input_block.size());

//  output_block.clear();
//  output_block = state_block;
//  output_block.reserve(input_block.size());

//  for(unsigned int n = 0; n<(input_block.size()-state_block.size());n++){
//    output_block.push_back(input_block[n]);
//  }

//  std::vector<float> temp(&input_block[input_block.size()-state_block.size()],&input_block[input_block.size()]);
//  state_block = temp;

}
