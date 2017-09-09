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

t, mx, Dx = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2), unpack=True, dtype=float)

max_D = max(sorted(Dx)[0 : round(len(t) * 0.98)])
min_m = min(sorted(mx)[round(len(t) * 0.02) : len(t)])

from pylab import *

f = plt.figure(num=None, figsize=(7, 3.5), dpi=150, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=0.06, bottom=0.01, right=0.98, top=0.99, wspace=0.1)
ax1 = plt.subplot(111)

ls_m = (0, ())
ls_D =  (0, (5, 1))

params = {
            'color' : 'black', 
            'linewidth' : 1.5,
            }

ax1.plot(t, mx, **params, linestyle=ls_m)
ax1.plot(t, Dx, **params, linestyle=ls_D)

ax1.set_ylim(min_m, max_D * 1.1)

plt.savefig(outputfilename)


