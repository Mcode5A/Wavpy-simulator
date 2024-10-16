import matplotlib.pyplot as plt
import numpy as np
import sys
import os
HOME = os.getenv("HOME")
wavpyDIR = HOME + "/gnssr_analysis/trunk/src/waveform_pylib/wavpy/"
sys.path.insert(1,wavpyDIR)
import wavpy

angles_out = np.arange(360)
receiver = wavpy.RF_FrontEnd()
receiver.set_receiver_params(0.0, 10.0, 3.0, 12000000.0, 0)

angles = [0.0, 45.0, 50.0, 60.0, 90.0, 130.0, 180.0, 213.2, 250.0, 267.567, 290.0, 310.45, 330.9 ]
pattern = [0.0, -10.0, -12.0, -13.0, -30.0, -20.0, -30.0, -30.0, -20.0, -15.43, -5.0, -20.0, -6.0 ]
receiver.set_antenna_pattern_interp( angles, pattern, -30.0 )
[pattern_E, pattern_H] = receiver.get_antenna_patterns()

receiver.dump_parameters()

plt.plot(angles_out, pattern_E)
plt.plot(angles, pattern, '.-')
plt.show()


