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
modelWAV.receiver_Down.set_receiver_params(15.0, 200.0, 3.0, 12000000.0, 0)
modelWAV.receiver_Down.set_antenna_patterns_interp(anglesE, patternE, anglesH, patternH, -40.0)
modelWAV.receiver_Down.set_antenna_orientation_BF_EH([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0])
modelWAV.receiver_Up.set_receiver_params(15.0, 10.0, 3.0, 12000000.0, 0)
modelWAV.receiver_Up.set_antenna_patterns_interp(anglesE, patternE, anglesH, patternH, -40.0)
modelWAV.receiver_Up.set_antenna_orientation_BF_EH([0.0, 1.0, 0.0], [1.0, 0.0, 0.0])
#Plot antenna diagram and dump parameters
[pattern_E_out, pattern_H_out] = modelWAV.receiver_Down.get_antenna_patterns()
angles_out = np.concatenate((np.arange(0,181,1), np.arange(-179,0,1)), axis=0)
plt.plot(angles_out, pattern_E_out - pattern_E_out[0])
plt.plot(angles_out, pattern_H_out - pattern_H_out[0])
plt.plot(anglesE, patternE, '.')
plt.plot(anglesH, patternH, '.')
plt.show()
modelWAV.receiver_Down.dump_parameters()
modelWAV.receiver_Up.dump_parameters()

#CASE 1: Receiver at H=3km, elevation=90deg, MSS=0.02, sea water
modelWAV.geometry.set_LongLatHeight_Receiver([0.0, 0.0, 3.0])
modelWAV.geometry.set_LongLatHeight_Transmitter([0.0, 0.0, 20000.0])
modelWAV.geometry.compute_specular_point(1)
modelWAV.surface.mss_1 = 0.01
modelWAV.surface.mss_2 = 0.01
modelWAV.interferometric_flag = True
modelWAV.compute_waveform()
range_wav = modelWAV.waveform_POW.get_range_waveform(256)
wav_Case1_Comp = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.gnss_signal.weight_IM = 0.0
modelWAV.interferometric_flag = False
modelWAV.compute_waveform()
wav_Case1_CA = modelWAV.waveform_POW.get_waveform(256)
plt.plot(range_wav, 10*np.log10(wav_Case1_Comp))
plt.plot(range_wav, 10*np.log10(wav_Case1_CA))
plt.show()

#CASE 2: Receiver at H=3km, elevation=90deg, MSS=0.02, sea water, roll=20deg
modelWAV.gnss_signal.weight_PY = 1.0
modelWAV.gnss_signal.weight_M = 1.0
modelWAV.geometry.set_inertials(20.0, 0.0, 0.0)
modelWAV.interferometric_flag = True
modelWAV.compute_waveform()
range_wav = modelWAV.waveform_POW.get_range_waveform(256)
wav_Case2_Comp = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.interferometric_flag = False
modelWAV.compute_waveform()
wav_Case2_CA = modelWAV.waveform_POW.get_waveform(256)
plt.plot(range_wav, 10*np.log10(wav_Case2_Comp))
plt.plot(range_wav, 10*np.log10(wav_Case2_CA))
plt.show()

#CASE 3: Receiver at H=3km, elevation=90deg, MSS=0.02, sea water, pitch=20deg
modelWAV.gnss_signal.weight_PY = 1.0
modelWAV.gnss_signal.weight_M = 1.0
modelWAV.geometry.set_inertials(0.0, 20.0, 0.0)
modelWAV.interferometric_flag = True
modelWAV.compute_waveform()
range_wav = modelWAV.waveform_POW.get_range_waveform(256)
wav_Case3_Comp = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.interferometric_flag = False
modelWAV.compute_waveform()
wav_Case3_CA = modelWAV.waveform_POW.get_waveform(256)
plt.plot(range_wav, 10*np.log10(wav_Case3_Comp))
plt.plot(range_wav, 10*np.log10(wav_Case3_CA))
plt.show()

#CASE 4: Receiver at H=300km, elevation=90deg, MSS=0.02, sea water
modelWAV.gnss_signal.weight_PY = 1.0
modelWAV.gnss_signal.weight_M = 1.0
modelWAV.geometry.set_inertials(0.0, 0.0, 0.0)
modelWAV.geometry.set_LongLatHeight_Receiver([0.0, 0.0, 300.0])
modelWAV.geometry.compute_specular_point(1)
modelWAV.interferometric_flag = True
modelWAV.curvature_approx_flag = True
modelWAV.compute_waveform()
range_wav = modelWAV.waveform_POW.get_range_waveform(256)
wav_Case4_Comp = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.interferometric_flag = False
modelWAV.compute_waveform()
wav_Case4_CA = modelWAV.waveform_POW.get_waveform(256)
plt.plot(range_wav, 10*np.log10(wav_Case4_Comp))
plt.plot(range_wav, 10*np.log10(wav_Case4_CA))
plt.show()

#CASE 5: Receiver at H=300km, elevation=90deg, MSS=0.01, sea water
modelWAV.gnss_signal.weight_PY = 1.0
modelWAV.gnss_signal.weight_M = 1.0
modelWAV.surface.mss_1 = 0.005
modelWAV.surface.mss_2 = 0.005
modelWAV.interferometric_flag = True
modelWAV.compute_waveform()
range_wav = modelWAV.waveform_POW.get_range_waveform(256)
wav_Case5_Comp = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.interferometric_flag = False
modelWAV.compute_waveform()
wav_Case5_CA = modelWAV.waveform_POW.get_waveform(256)
plt.plot(range_wav, 10*np.log10(wav_Case5_Comp))
plt.plot(range_wav, 10*np.log10(wav_Case5_CA))
plt.show()

#CASE 6: Receiver at H=200m, elevation=90deg, MSS=0.02, sea water
modelWAV.gnss_signal.weight_PY = 1.0
modelWAV.gnss_signal.weight_M = 1.0
modelWAV.surface.mss_1 = 0.01
modelWAV.surface.mss_2 = 0.01
modelWAV.geometry.set_LongLatHeight_Receiver([0.0, 0.0, 0.2])
modelWAV.geometry.compute_specular_point(1)
modelWAV.curvature_approx_flag = False
modelWAV.interferometric_flag = True
modelWAV.compute_waveform()
range_wav = modelWAV.waveform_POW.get_range_waveform(256)
wav_Case6_Comp = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.interferometric_flag = False
modelWAV.compute_waveform()
wav_Case6_CA = modelWAV.waveform_POW.get_waveform(256)
plt.plot(range_wav, 10*np.log10(wav_Case6_Comp))
plt.plot(range_wav, 10*np.log10(wav_Case6_CA))
plt.show()

#CASE 7: Receiver at H=200m, elevation=90deg, MSS=0.01, sea water
modelWAV.gnss_signal.weight_PY = 1.0
modelWAV.gnss_signal.weight_M = 1.0
modelWAV.surface.mss_1 = 0.005
modelWAV.surface.mss_2 = 0.005
modelWAV.interferometric_flag = True
modelWAV.compute_waveform()
range_wav = modelWAV.waveform_POW.get_range_waveform(256)
wav_Case7_Comp = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_PY = 0.0
modelWAV.gnss_signal.weight_M = 0.0
modelWAV.interferometric_flag = False
modelWAV.compute_waveform()
wav_Case7_CA = modelWAV.waveform_POW.get_waveform(256)
modelWAV.gnss_signal.weight_CA = 0.0
modelWAV.gnss_signal.weight_B1I = 1.0
modelWAV.compute_waveform()
wav_Case7_BeiDou = modelWAV.waveform_POW.get_waveform(256)
plt.plot(range_wav, 10*np.log10(wav_Case7_Comp))
plt.plot(range_wav, 10*np.log10(wav_Case7_CA))
plt.plot(range_wav, 10*np.log10(wav_Case7_BeiDou))
plt.show()
