#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <cmath>

#define pi 3.14159265358979323846

void pll(std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float &state_integrator, float &state_phaseEst, float &state_feedbackI, float &state_feedbackQ, float &state_ncoLast, int &state_trigOffset, float normBandwidth, std::vector<float> &ncoOut){
  float Cp = 2.666;
	float Ci = 3.555;

  float Kp = (normBandwidth)*Cp;
  float Ki = (normBandwidth*normBandwidth)*Ci;

  ncoOut.clear();
  ncoOut.resize(pllIn.size()+1);

  float integrator = state_integrator;
	float phaseEst = state_phaseEst;
	float feedbackI = state_feedbackI;
	float feedbackQ = state_feedbackQ;
  ncoOut[0] = state_ncoLast;
	int trigOffset = state_trigOffset;

  float errorI;
  float errorQ;
  float errorD;
  float trigArg;



  for(unsigned int k=0; k<pllIn.size(); k++){

    //phase detector
		errorI = pllIn[k] * (+feedbackI);  //complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ); //feedback complex exponential

		//four-quadrant arctangent discriminator for phase error detection
		errorD = std::atan2(errorQ, errorI);

		//loop filter
		integrator = integrator + Ki*errorD;

		//update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

    //internal oscillator
		trigOffset += 1;
		trigArg = 2*pi*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);
  }

  state_integrator = integrator;
  state_phaseEst = phaseEst;
  state_feedbackI = feedbackI;
  state_feedbackQ = feedbackQ;
  state_ncoLast = ncoOut[ncoOut.size()-1];
  state_trigOffset = trigOffset + pllIn.size();

  //ncoOut.assign(ncoOut.begin(),ncoOut.end()-1);

  ncoOut = std::vector<float>(ncoOut.begin(), ncoOut.end()-1);

}

void pllIQ(std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float &state_integrator, float &state_phaseEst, float &state_feedbackI, float &state_feedbackQ, float &state_ncoLast, float &state_ncoQLast,int &state_trigOffset, float normBandwidth, std::vector<float> &ncoOut, std::vector<float> &ncoOutQ){
  float Cp = 2.666;
	float Ci = 3.555;

  float Kp = (normBandwidth)*Cp;
  float Ki = (normBandwidth*normBandwidth)*Ci;

  ncoOut.clear();
  ncoOut.resize(pllIn.size()+1);
  ncoOutQ.clear();
  ncoOutQ.resize(pllIn.size()+1);

  float integrator = state_integrator;
	float phaseEst = state_phaseEst;
	float feedbackI = state_feedbackI;
	float feedbackQ = state_feedbackQ;
  ncoOut[0] = state_ncoLast;
  ncoOutQ[0] = state_ncoQLast;
	int trigOffset = state_trigOffset;

  float errorI;
  float errorQ;
  float errorD;
  float trigArg;



  for(unsigned int k=0; k<pllIn.size(); k++){

    //phase detector
		errorI = pllIn[k] * (+feedbackI);  //complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ); //feedback complex exponential

		//four-quadrant arctangent discriminator for phase error detection
		errorD = std::atan2(errorQ, errorI);

		//loop filter
		integrator = integrator + Ki*errorD;

		//update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

    //internal oscillator
		trigOffset += 1;
		trigArg = 2*pi*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);
    ncoOutQ[k+1] = std::sin(trigArg*ncoScale + phaseAdjust);
  }

  state_integrator = integrator;
  state_phaseEst = phaseEst;
  state_feedbackI = feedbackI;
  state_feedbackQ = feedbackQ;
  state_ncoLast = ncoOut[ncoOut.size()-1];
  state_ncoQLast = ncoOutQ[ncoOutQ.size()-1];
  state_trigOffset = trigOffset + pllIn.size();

  //ncoOut.assign(ncoOut.begin(),ncoOut.end()-1);

  ncoOut = std::vector<float>(ncoOut.begin(), ncoOut.end()-1);
  ncoOutQ = std::vector<float>(ncoOutQ.begin(), ncoOutQ.end()-1);
}

void mono_stereo_part(){



}
