import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab

filename = u"../output/test1_estimateAvg_0.txt"
#filename = u"../output/test2_estimateAvg.txt"
#t, x, Dx, xhat, err, err2, D, xhatU, errU, errU2, DU, xhatU_, errU_, errU_2, DU_ = np.loadtxt(filename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), unpack=True, dtype=float)
t, x, Dx, xhat, mErr, DErr, DErrTheor, xhatU, mErrU, DErrU, DErrUTheor = np.loadtxt(filename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True, dtype=float)

from pylab import *

f = plt.figure(num=None, figsize=(7, 5), dpi=150, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=0.06, bottom=0.01, right=0.98, top=0.99, wspace=0.1)
#ax1 = plt.subplot(411)
ax1 = plt.subplot(111)
ax1.plot(t, x, '-', color = 'black', linewidth = 1.0)
ax1.plot(t, mErr, '-', color = 'blue', linewidth = 1.0)
ax1.plot(t, DErr, '-', color = 'red', linewidth = 1.0)
ax1.plot(t, DErrTheor, '-', color = 'green', linewidth = 1.0)
ax1.plot(t, mErrU, '--', color = 'blue', linewidth = 1.0)
ax1.plot(t, DErrU, '--', color = 'red', linewidth = 1.0)
ax1.plot(t, DErrUTheor, '--', color = 'green', linewidth = 1.0)
#ax1.plot(t, errU_, ':', color = 'blue', linewidth = 1.0)
#ax1.plot(t, DU_, ':', color = 'red', linewidth = 1.0)
#ax1.plot(t, errU_2, ':', color = 'green', linewidth = 1.0)

ax1.set_ylim(0,6)


#ax1.set_ylim(-0.5,0.5)

plt.show();
