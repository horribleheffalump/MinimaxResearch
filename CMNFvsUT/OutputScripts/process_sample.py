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

t, x, y = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2), unpack=True, dtype=float)

max_x = max(sorted(x)[0 : round(len(t) * 0.98)])
min_x = min(sorted(x)[round(len(t) * 0.02) : len(t)])

max_y = max(sorted(y)[0 : round(len(t) * 0.98)])
min_y = min(sorted(y)[round(len(t) * 0.02) : len(t)])

from pylab import *

f = plt.figure(num=None, figsize=(7, 3.5), dpi=150, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=0.06, bottom=0.01, right=0.98, top=0.99, wspace=0.1)
ax1 = plt.subplot(111)

ls_x = (0, ()) #solid
ls_y =  (0, (1, 1)) #dots

params = {
            'color' : 'black', 
            'linewidth' : 1.5,
            }

ax1.plot(t, x, **params, linestyle=ls_x)
ax1.plot(t, y, **params, linestyle=ls_y)

#ax1.set_ylim(min(min_x, min_y), max(max_x, max_y))

plt.savefig(outputfilename)



