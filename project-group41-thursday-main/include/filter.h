/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void block_processing_with_state(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb);
void filter_signal(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb);
void filter_ds(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb, const int &factor);
void lpf_block_processing(std::vector<uint8_t> &data, unsigned int block_size, float Fc, float Fs, int N_taps, std::vector<float> &filtered_data);
void filter_ds_with_us(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &state, std::vector<float> &yb, const int &factor,const int &us);
void impulseResponseBPF(float Fs, float Fb,float fe, unsigned short int num_taps, std::vector<float> &h);
void filter_ds_ptr(const std::vector<float> &h, float* &xb, std::vector<float> &state, std::vector<float> &yb, const int &factor,const unsigned int block_size);
void allpass(const std::vector<float> &input_block, std::vector<float> &state_block, std::vector<float> &output_block);
#endif // DY4_FILTER_H
