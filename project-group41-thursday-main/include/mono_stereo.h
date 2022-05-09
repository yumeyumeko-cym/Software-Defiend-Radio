#ifndef mono_stereo_H
#define mono_stereo_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void pll(std::vector<float> &, float, float, float, float, float &, float &, float &, float &, float &, int &, float, std::vector<float> &);
void pllIQ(std::vector<float> &, float, float, float, float, float &, float &, float &, float &, float &, float &, int &, float, std::vector<float> &, std::vector<float> &);
#endif // DY4_FILTER_H
