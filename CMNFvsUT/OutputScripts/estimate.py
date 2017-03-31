import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab

filename = u"../output/estimateAverage.txt"
t, x, xhat, err, sigma = np.loadtxt(filename, delimiter = ' ', usecols=(0,1,2,3,4), unpack=True, dtype=float)

from pylab import *

f = plt.figure(num=None, figsize=(7, 5), dpi=150, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=0.06, bottom=0.01, right=0.98, top=0.99, wspace=0.1)
#ax1 = plt.subplot(411)
ax1 = plt.subplot(111)
ax1.plot(t, x, '-', color = 'black', linewidth = 1.0)
ax1.plot(t, err, '-', color = 'blue', linewidth = 1.0)
ax1.plot(t, 3*sigma, '-', color = 'red', linewidth = 1.0)


plt.show();
