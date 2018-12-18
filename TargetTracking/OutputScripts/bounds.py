import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os
import pandas as pd
#inputfilename = sys.argv[1]
#outputfilename = sys.argv[2]
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)

for i in range(0,5):
    inputfilename = "D:/results/cont_a_1.0_two/TargetTracking_average_" + str(i) + ".txt"
    outputfilename = "D:/results/cont_a_1.0_two/TargetTracking_bulk_" + str(i) + ".png"

    bulkfilename = "D:/results/cont_a_1.0_two/TargetTracking_bulk_" + str(i) + ".txt"


    t, x, Dx, xhat_1, mErr_1, DErr_1, DErrTheor_1, xhat_2, mErr_2, DErr_2, DErrTheor_2, xhat_3, mErr_3, DErr_3, DErrTheor_3, xhat_4, mErr_4, DErr_4, DErrTheor_4 = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), unpack=True, dtype=float)

    #,

    from pylab import *
    

    f = plt.figure(num=None, figsize=cm2inch((12,9)), dpi=200)
    plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
    ax = plt.subplot(111)

    T = 150



   
    data = pd.read_csv(bulkfilename, delimiter = " ", header=0, dtype=float, engine='python')
    for i in range(0,1000):
        plt.plot(np.arange(1,50), data.iloc[i,0:49], color='black', alpha=0.1)

    ax.plot(t[0:T], x[0:T]-mErr_4[0:T], c = 'k')
    ax.plot(t[0:T], x[0:T]-mErr_4[0:T] - 5 * np.sqrt(DErrTheor_4[0:T]), c = 'orange')
    ax.plot(t[0:T], x[0:T]-mErr_4[0:T] + 5 * np.sqrt(DErrTheor_4[0:T]), c = 'orange')

    #plt.show()
    plt.savefig(outputfilename)



