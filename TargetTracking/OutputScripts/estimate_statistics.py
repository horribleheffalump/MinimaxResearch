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
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

for i in range (0,5):
    inputfilename = "D:/results/cont_a_1.0_two/TargetTracking_average_" + str(i) + ".txt"
    outputfilename = "D:/results/cont_a_1.0_two/TargetTracking_average_" + str(i) + ".png"




    t, x, Dx, xhat_1, mErr_1, DErr_1, DErrTheor_1, xhat_2, mErr_2, DErr_2, DErrTheor_2, xhat_3, mErr_3, DErr_3, DErrTheor_3 = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), unpack=True, dtype=float)

    #, xhat_4, mErr_4, DErr_4, DErrTheor_4 

    from pylab import *
    
    f = plt.figure(num=None, figsize=cm2inch((12,9)), dpi=200)
    plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
    ax = plt.subplot(111)

    T = 150

    #ax.plot(t[1:T], mErr_1[1:T], c = 'k')
    #ax.plot(t[1:T], DErr_1[1:T], c = 'r')

    #ax.plot(t[1:T], mErr_2[1:T], c = 'k')
    ax.plot(t[1:T], DErr_2[1:T], c = 'g')

    #ax.plot(t[1:T], mErr_3[1:T], c = 'k')
    ax.plot(t[1:T], DErr_3[1:T], c = 'b')

    #ax.plot(t[1:T], mErr_4[1:T], c = 'k')
    #ax.plot(t[1:T], DErr_4[1:T], c = 'orange')

    plt.savefig(outputfilename)



