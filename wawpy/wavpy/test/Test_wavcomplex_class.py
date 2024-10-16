import numpy as np
import scipy
from scipy.io import netcdf
import sys
import os
HOME = os.getenv("HOME")
wavpyDIR = HOME + "/gnssr_analysis/trunk/src/waveform_pylib/wavpy/"
sys.path.insert(1,wavpyDIR)
import wavpy
import matplotlib.pyplot as plt

data_path = "../../src/data_share/"
data_file = data_path + "Greenland_sample.WAV.nc"
f = netcdf.netcdf_file(data_file, 'r')
waveform_I_comps = f.variables['I_Waveform']
waveform_Q_comps = f.variables['Q_Waveform']
waveform_links = f.variables['Link']
waveform_sow = f.variables['SecondOfWeek']
waveform_msec = f.variables['Millisecond']
size_cluster = 1000
size_waveform = 64
mywavClusterL1 = wavpy.Waveform_complex_cluster()
mywavClusterL1.initialize(size_cluster, size_waveform)
mywavClusterL2 = wavpy.Waveform_complex_cluster()
mywavClusterL2.initialize(size_cluster, size_waveform)

#In addition to its cc equivalent version, integrated waveforms are plotted
mywaveform = wavpy.Waveform_power()
mywaveform.set_sampling_rate(20000000)

for second in xrange(waveform_sow[0], waveform_sow[-1]):
  for msec in xrange(0, size_cluster):
    wavIcomp_link1 = waveform_I_comps[(waveform_sow[:] == second) & (waveform_msec[:] == msec) & (waveform_links[:] == 1) , :]
    wavQcomp_link1 = waveform_Q_comps[(waveform_sow[:] == second) & (waveform_msec[:] == msec) & (waveform_links[:] == 1) , :]
    if ((wavIcomp_link1[:].size == size_waveform) & (wavQcomp_link1[:].size == size_waveform)):
      mywavClusterL1.add_waveform_GOLD(wavIcomp_link1.reshape(size_waveform), wavQcomp_link1.reshape(size_waveform), msec)
      
    wavIcomp_link2 = waveform_I_comps[(waveform_sow[:] == second) & (waveform_msec[:] == msec) & (waveform_links[:] == 2) , :]
    wavQcomp_link2 = waveform_Q_comps[(waveform_sow[:] == second) & (waveform_msec[:] == msec) & (waveform_links[:] == 2) , :]
    if ((wavIcomp_link2[:].size == size_waveform) & (wavQcomp_link2[:].size == size_waveform)):
      mywavClusterL2.add_waveform_GOLD(wavIcomp_link2.reshape(size_waveform), wavQcomp_link2.reshape(size_waveform), msec)
      
  #Plot integrated waveforms
  int_wav = mywavClusterL1.integrate_waveforms(10, 64)
  plt.plot(int_wav)
  int_wav = mywavClusterL2.integrate_waveforms(10, 64)
  plt.plot(int_wav)
      
  mywavClusterL1.store_phasor_wavs((size_waveform/2)+1)
  [phasorI, phasorQ, valid_phasor] = mywavClusterL1.get_phasor(size_cluster, size_cluster, size_cluster)
  phasesL1 = np.arctan2(phasorQ, phasorI)
  mywavClusterL1.initialize(size_cluster, size_waveform)
  mywavClusterL2.store_phasor_wavs((size_waveform/2)+1)
  mywavClusterL2.counterrot_phasor(phasesL1, valid_phasor)
  [phasorI, phasorQ, valid_phasor] = mywavClusterL2.get_phasor(size_cluster, size_cluster, size_cluster)
  phasesL2 = np.arctan2(phasorQ, phasorI)*np.float32(valid_phasor)
  mywavClusterL2.initialize(size_cluster, size_waveform)
  for msec in xrange(0, size_cluster):
    print "PHASE %d %d %f" % (second, msec, phasesL2[msec])

#Show plot with integrated waveforms
plt.show()
