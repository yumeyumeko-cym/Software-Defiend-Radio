import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD

from fmMonoBasic import filtering, my_own_coeff
# for take-home add your functions

def block_processing_with_state(h, xb, state):
	yb = np.zeros(len(xb))
	for n in range(len(yb)):
		for k in range(len(h)):
			if n - k >= 0:
				yb[n] += h[k] * xb[n-k]
			else:
				yb[n] += h[k] * state[len(state)+(n-k)]
	state = xb[(len(xb)-len(state)):(len(xb))]
	return yb, state


def filter_signal(coeff, xb, state):
	filtered_signal = np.array([])
	signal, next_state = block_processing_with_state(coeff, xb, state)
	filtered_signal = np.concatenate((filtered_signal,signal))
	#filtered_signal = filtered_signal[:(len(x)+len(coeff)-1)]
	return filtered_signal, next_state

def impulseResponseBPF(Fs, Fb, Fe, num_taps,h):
    h=np.zeros(num_taps)
    norm_center=((Fe+Fb)/2)/(Fs/2)
    norm_pass=((Fe-Fb)*2)/Fs

    for i in range (num_taps):
        if(i == (num_taps-1)/2):
            h[i] = norm_pass
        else:
            h[i] = math.sin(pi*(norm_pass/2)*(i-(num_taps-1)/2))/(pi*(norm_pass/2)*(i-(num_taps-1)/2))

        h[i]=h[i]*math.cos(i*pi*norm_center)
        h[i] = h[i] * (math.sin(i*pi/num_taps))**2
    return h


def filter_ds(h, xb, state,factor):
    yb = np.zeros(int(len(xb)/factor))
    for n in range (len(yb)):
        count=0
        for k in range (len(h)):
            if(factor*n - k >= 0 and (factor*n-k < len(xb))):
                yb[n] += h[k] * xb[factor*n-k]

            else:
                #yb[n] += h[k] * state[len(state)+(n-k)]
                yb[n] += h[k] * state[len(state)-count-1]
                count += 1

    state = xb[(len(xb)-len(state)):(len(xb))]
    return yb, state

def filter_ds_with_us(h, xb, state, factor,us):
    yb = np.zeros(int(us*len(xb)/factor))
    phase=0

    for n in range (len(yb)):

        phase=(n*factor)%us;
        for k in range (int(len(h)/us)):
            if((n*factor-phase)/us >= 0 and (n*factor-phase)/us < len(xb)):
                yb[n] += (h[int(phase+k*us)]*xb[int((n*factor-phase)/us)])*us

            else:
				#yb[n] += h[phase+k*us] * state[state.size()+((n*factor-phase)/us)-(phase+k*us)];
                yb[n] += h[phase+k*us] * state[len(state)+(n-k)]
        state = xb[(len(xb)-len(state)):(len(xb))]
        return yb, state

def update_method(I, Q, last_q = 0.0, last_i = 0.0):
	fm_demod = np.empty(len(I))

	for k in range(len(I)):
		der_q = Q[k] - last_q
		der_i = I[k] - last_i

		last_q = Q[k]
		last_i = I[k]

		fm_demod[k] = (I[k] * der_q - Q[k] * der_i)/(I[k]**2 + Q[k]**2)

	return fm_demod, last_q, last_i

if __name__ == "__main__":


     # read the raw IQ data from the recorded file
    # IQ data is normalized between -1 and +1 and interleaved
    in_fname = "../data/iq_samples.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in float32 format")
    iq_data = (np.float32(raw_data) - 128.0)/128.0
    rf_Fs = 2.4e6
    rf_Fc = 100e3
    rf_taps = 151
    rf_decim = 10

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

    audio_Fc = 16e3

    mode=0
    if(mode == 0):
            rf_Fs = 2.4e6
            rf_decim = 10
            audio_decim = 5
            audio_Fs = 48e3
            mono_us = 1
            num_taps=101

    elif(mode == 1):
            rf_Fs = 1.92e6
            rf_decim = 8
            audio_decim = 5
            audio_Fs = 48e3
            mono_us = 1
            num_taps=101

    elif(mode == 2):
            rf_Fs = 2.4e6
            rf_decim = 10
            audio_decim = 800
            audio_Fs = 441e2
            mono_us = 147
            num_taps=151

    elif(mode == 3):
            rf_Fs = 1.152e6
            rf_decim = 4
            audio_decim = 320
            audio_Fs = 441e2
            mono_us = 49
            num_taps=101

    #block_size = 1024 * rf_decim * audio_decim * 2
    block_size = 307200
    block_count = 0

	# states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    filter_state = np.zeros(rf_taps-1)
    state_phase = 0
    last_q = 0
    last_i = 0
    audio_data = np.array([]) # used to concatenate filtered blocks (audio data)
    audio_coeff = np.array([])

    #set up the subfigures for plotting
    subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
    plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace = .6)

    print(len(iq_data))
    print("blocksize\n")
    print(block_size)

    while(block_count+1)*block_size < len(iq_data):
        print('Processing block ' + str(block_count))
        # filter to extract the FM channel (I samples are even, Q samples are odd)
        i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                       iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
                       zi=state_i_lpf_100k)
        q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                       iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
                       zi=state_q_lpf_100k)

        # downsample the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]


        fm_demod, last_q, last_i = update_method(i_ds, q_ds, last_q, last_i)
        if(mode==0 or mode==1):
            audio_coeff = signal.firwin(num_taps, audio_Fc/(rf_Fs/rf_decim/2), window=('hann'))
            audio_filt, filter_state = filter_ds(audio_coeff, fm_demod, filter_state,audio_decim)
            #audio_filt, filter_state = filter_signal(audio_coeff, fm_demod, filter_state)
            #audio_block = audio_filt[::5]

        else:
            #audio_coeff = signal.firwin(num_taps, audio_Fc/(rf_Fs/rf_decim/2), window=('hann'))
            audio_coeff = my_own_coeff(audio_Fc,audio_Fs*mono_us,num_taps*mono_us)
            audio_filt, filter_state = filter_ds_with_us(audio_coeff, fm_demod, filter_state, audio_decim,mono_us)


        audio_data = np.concatenate((audio_data, audio_filt))
        #
        if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
            ax0.clear()
            fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
            		'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
            fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
            fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as neede
            ax1.clear()
            fmPlotPSD(ax1, audio_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Mono')
			# plot PSD of selected block after downsampling mono audio
			# ... change as needed
            ax2.clear()
            fmPlotPSD(ax2, audio_filt, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')
			# save figure to file
            fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

        block_count += 1

    print('Finished processing all the blocks from the recorded I/Q samples')

    # write audio data to file
    out_fname = "../data/fmMonoBlock.wav"
    wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
    print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

    # uncomment assuming you wish to show some plots
    plt.show()
