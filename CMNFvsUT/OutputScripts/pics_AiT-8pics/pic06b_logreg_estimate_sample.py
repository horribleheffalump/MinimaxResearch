
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
from multiplypoints import *
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)



inputfilename = "D:/results/logreg-zero_all_but_CMNF/LogisticModelZero_sample_0_0.txt"

t, x, err_umf, err_ut = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,3,4), unpack=True, dtype=float)
t, err_umf, err_ut = multx(t), multy(err_umf), multy(err_ut)
from pylab import *

f = plt.figure(num=None, figsize=cm2inch((12,9)), dpi=200)
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
ax1.plot(t, err_umf, **params_UMF, linestyle=ls_x, label='$\hat{x}_t-x_t$ УМНФ')
ax1.plot(t, err_ut, **params_UT, linestyle=ls_x, label='$\hat{x}_t-x_t$ CT-фильтр')

ax1.set_ylim(-1000.0, 500.0)
#ax1.legend()
plt.tight_layout()

#float commas
rc('text.latex',preamble=r'\usepackage{icomma}')
import matplotlib as mpl
import locale
locale.setlocale(locale.LC_ALL, "rus_RUS")
ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useLocale=True, useMathText=True))
ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useLocale=True, useMathText=True))

plt.show()



