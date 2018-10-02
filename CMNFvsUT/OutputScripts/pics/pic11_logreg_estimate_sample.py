
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os
from multiplypoints import *
inputfilename = "D:/results/logreg-zero/LogisticModelZero_sample_0.txt"

t, x, err_umf, err_ut = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,3,4), unpack=True, dtype=float)
t, err_umf, err_ut = multx(t), multy(err_umf), multy(err_ut)
from pylab import *

f = plt.figure(num=None, figsize=(5,5), dpi=200)
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

ax1.set_ylim(-1000.0, 1000.0)

plt.tight_layout()
plt.show()



