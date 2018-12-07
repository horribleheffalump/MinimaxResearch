import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os

#inputfilename = sys.argv[1] 
#outputfilename = sys.argv[2]

for i in range (0,5):
    inputfilename = "D:/results/cont/TargetTracking_average_" + str(i) + ".txt"
    outputfilename = "D:/results/cont/TargetTracking_average_" + str(i) + ".png"




    t, x, Dx, xhat_1, mErr_1, DErr_1, DErrTheor_1, xhat_2, mErr_2, DErr_2, DErrTheor_2 = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True, dtype=float)

    from pylab import *

    f = plt.figure(num=None, figsize=(14, 3.5), dpi=150, facecolor='w', edgecolor='k')
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])     
    gs.update(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1, hspace=0.1)
    ax_1 = plt.subplot(gs[0])
    ax_2 = plt.subplot(gs[1])

    T = 600

    ax_1.plot(t[0:T], mErr_1[0:T], c = 'k')
    ax_1.plot(t[0:T], DErr_1[0:T], c = 'r')

    ax_2.plot(t[0:T], mErr_2[0:T], c = 'k')
    ax_2.plot(t[0:T], DErr_2[0:T], c = 'g')


    plt.savefig(outputfilename)



