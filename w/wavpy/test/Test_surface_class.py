import numpy as np
import sys
import os
HOME = os.getenv("HOME")
wavpyDIR = HOME + "/gnssr_analysis/trunk/src/waveform_pylib/wavpy/"
sys.path.insert(1,wavpyDIR)
import wavpy

mysurf = wavpy.Reflecting_surface()
mysurf.dump_parameters()
mysurf.c21_coeff = 1.0
mysurf.c03_coeff = 2.0
mysurf.wind_Magnitude = 3.0
mysurf.wind_Azimuth = 4.0
mysurf.epsilon_sea_ice(30.0)
mysurf.dump_parameters()
[rco, rcross] = mysurf.compute_Rfresnel_circular((3.14159/6.0), [1.0, 0.0])
print "Fresnel CIRC: %f %f %f %f" % (rco[0], rco[1], rcross[0], rcross[1])
mysurf.compute_sea_spectrum(4, 100.5, 90.0, 1.0)
for i in range(4):
  for j in range(4):
    mysurf.get_sea_spectrum(i, j)
    


