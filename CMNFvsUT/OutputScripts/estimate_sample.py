
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os

#inputfilename = os.path.join(sys.argv[1], sys.argv[2])
#outputfilename = os.path.join(sys.argv[1], sys.argv[3])
inputfilename = sys.argv[1]
outputfilename = sys.argv[2]

#t, err_umf, err_ut = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,3,4), unpack=True, dtype=float)
t, err_umf, err_ut = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,2,3), unpack=True, dtype=float)

max_1 = max(sorted(err_umf)[0 : round(len(t) * 0.98)])
min_1 = min(sorted(err_umf)[round(len(t) * 0.02) : len(t)])

max_2 = max(sorted(err_ut)[0 : round(len(t) * 0.98)])
min_2 = min(sorted(err_ut)[round(len(t) * 0.02) : len(t)])

from pylab import *

f = plt.figure(num=None, figsize=(7, 3.5), dpi=150, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
ax1 = plt.subplot(111)

ls_x = (0, ()) #solid

alpha_UT = 0.4
params_UT = {
            'color' : 'black', 
            'alpha' : alpha_UT,
            'linewidth' : 4.0,
            }
alpha_UMF = 0.7
params_UMF = {
            'color' : 'black', 
            'alpha' : alpha_UMF,
            'linewidth' : 2.5,
            }

ax1.plot(t, err_umf, **params_UMF, linestyle=ls_x)
ax1.plot(t, err_ut, **params_UT, linestyle=ls_x)

ax1.set_ylim(min(min_1, min_1), max(max_2, max_2))

plt.savefig(outputfilename)



