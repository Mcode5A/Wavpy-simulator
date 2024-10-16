import matplotlib.pyplot as plt
import numpy as np
import sys
import os
HOME = os.getenv("HOME")
wavpyDIR = HOME + "/gnssr_analysis/trunk/src/waveform_pylib/wavpy/"
sys.path.insert(1,wavpyDIR)
import wavpy

modelWAV = wavpy.ZaVoModel_GNSSR()

#Antenna diagram PIRA flight 2
anglesE = [-90.0, -53.0, -40.0, -25.0, 0.0, 20.0, 35.0, 45.0, 90.0]
patternE = [-40.0, -20.0, -11.0, -4.0, 0.0, -3.0, -9.0, -15.0, -40.0]
anglesH = [-60.0, -52.0, -45.0, -37.0, -30.0, -28.0, -25.0, -20.0, -15.0, 0.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 45.0, 60.0]
patternH = [-40.0, -20.0, -15.0, -13.0, -16.0, -20.0, -40.0, -14.0, -6.0, 0.0, -3.0, -8.0, -15.0, -40.0, -18.0, -15.7, -20.0, -40.0]
modelWAV.receiver_Down.set_antenna_patterns_interp(anglesE, patternE, anglesH, patternH, -40.0)
modelWAV.receiver_Down.set_antenna_orientation_EH([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0])
modelWAV.receiver_Down.set_receiver_params(15.0, 200.0, 3.0, 12000000.0, 0) #Antenna peak gain, ant temp, F, BW, isotropic (1)
modelWAV.receiver_Up.set_antenna_patterns_interp(anglesE, patternE, anglesH, patternH, -40.0)
modelWAV.receiver_Up.set_antenna_orientation_EH([0.0, 1.0, 0.0], [1.0, 0.0, 0.0])
modelWAV.receiver_Up.set_receiver_params(15.0, 10.0, 3.0, 12000000.0, 0)

#CASE: Receiver at H=3km, elevation=90deg, MSS=0.02, sea water
modelWAV.geometry.set_LongLatHeight_Receiver([0.0, 0.0, 3.0])
modelWAV.geometry.set_LongLatHeight_Transmitter([0.0, 0.0, 20000.0])
modelWAV.geometry.compute_specular_point(1)
modelWAV.surface.mss_1 = 0.01
modelWAV.surface.mss_2 = 0.01
maxval = modelWAV.compute_waveform(1)
wav_model = modelWAV.waveform_POW.get_waveform(256)
modelWAV.waveform_POW.compute_delays()
modelWAV.waveform_POW.dump_delays() #get delays from the model
autocorr_length = modelWAV.gnss_signal.lambda_size
[range_autocorr, autocorr] = modelWAV.gnss_signal.get_lambda_func(autocorr_length, autocorr_length)

#plt.plot(10*np.log10(wav_model))
#plt.show()

num_wavs = 100
waveform_array = np.zeros((256, num_wavs))
for index in range(num_wavs):
  noise_exp = np.random.exponential(1.0, 1024)
  noise_corr = np.convolve(noise_exp, (autocorr*autocorr), 'valid')
  waveform_array[:, index] = wav_model*noise_corr[0:256]/np.std(noise_corr)

#Compute average waveform and get delays
wav_mean = np.mean(waveform_array, axis=1)

wav_norm = np.float32(wav_mean/max(wav_mean))
mywaveform = wavpy.Waveform_power()
mywaveform.set_receiver_params(80000000, 12000000)
mywaveform.set_norm_waveform(wav_norm, max(wav_mean))
mywaveform.compute_delays()
mywaveform.dump_delays()



