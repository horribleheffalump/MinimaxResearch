import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif'})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble=r'\usepackage[T2A]{fontenc}')
rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex',preamble=r'\usepackage[russian]{babel}')
import pylab
import sys
import os
#from multiplypoints import *

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

for i in range (0,10):
    inputfilename = "D:/results/cont/TargetTracking_sample_"+str(i)+".txt"

    t, x, err_1, err_2, err_3 = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4), unpack=True, dtype=float)
    from pylab import *

    f = plt.figure(num=None, figsize=cm2inch((12,9)), dpi=200)
    plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
    ax1 = plt.subplot(111)

    ax1.plot(t, x, 'k', label='x')
    ax1.plot(t, x-err_1, 'r', label='CMNF')
    ax1.plot(t, x-err_2, 'g', label='MCMNF')
    ax1.plot(t, x-err_3, 'b', label='UKF')

    ax1.legend()

    plt.tight_layout()
    plt.show()

