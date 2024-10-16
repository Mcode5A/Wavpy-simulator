import numpy as np
import scipy
from scipy.io import netcdf
import sys
import os
HOME = os.getenv("HOME")
wavpyDIR = HOME + "/gnssr_analysis/trunk/src/waveform_pylib/wavpy/"
sys.path.insert(1,wavpyDIR)
import wavpy

data_path = "../../src/data_share/"
data_file = data_path + "GEOHALO_sample.WAVINTGEO.nc"
f = netcdf.netcdf_file(data_file, 'r')
waveforms_file = f.variables['Waveform']
pow_wav = np.array(waveforms_file[1,:], dtype=np.float32)
mywaveform = wavpy.Waveform_power()
mywaveform.set_float_waveform(pow_wav)
mywaveform.set_sampling_rate(20000000.0)
mywaveform.dump_norm_waveform()
mywaveform.compute_delays()
mywaveform.dump_delays()
print "%f %f %f %f %f %f" % (mywaveform.positionMax, mywaveform.sigma_posMax, mywaveform.powerMax, mywaveform.positionDer, mywaveform.sigma_posDer, mywaveform.power_posDer)
pow_wav_ext = np.zeros([320], dtype=np.float32)
pow_wav_ext[0:64] = pow_wav[0:64]
mywaveform.set_float_waveform(pow_wav_ext)
mywaveform.set_sampling_rate(20000000.0)
mywaveform.compute_delays()
mywaveform.dump_delays()
