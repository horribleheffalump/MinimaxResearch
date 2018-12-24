import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
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

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
#inputfilename = os.path.join(sys.argv[1], sys.argv[2])
#outputfilename = os.path.join(sys.argv[1], sys.argv[3])

inputfilename = "D:/results/cubic_1000/CubicSensor_average_0.txt"
t, x, Dx, xhat_UMF, mErr_UMF, DErr_UMF, DErrTheor_UMF = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4,5,6), unpack=True, dtype=float)

inputfilename = "D:/results/cubic_1000/CubicSensor_average_0.txt"
xhat_MUMF_1000, mErr_MUMF_1000, DErr_MUMF_1000, DErrTheor_MUMF_1000 = np.loadtxt(inputfilename, delimiter = ' ', usecols=(7,8,9,10), unpack=True, dtype=float)

inputfilename = "D:/results/cubic_100/CubicSensor_average_0.txt"
xhat_MUMF_100, mErr_MUMF_100, DErr_MUMF_100, DErrTheor_MUMF_100 = np.loadtxt(inputfilename, delimiter = ' ', usecols=(3,4,5,6), unpack=True, dtype=float)

inputfilename = "D:/results/cubic_25/CubicSensor_average_0.txt"
xhat_MUMF_25, mErr_MUMF_25, DErr_MUMF_25, DErrTheor_MUMF_25 = np.loadtxt(inputfilename, delimiter = ' ', usecols=(3,4,5,6), unpack=True, dtype=float)

from pylab import *

f = plt.figure(num=None, figsize=cm2inch((12,9)), dpi=200)
plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
ax = plt.subplot(111)

#ls_m = (0, ())
ls_D =  (0, ())
ls_D1 =  (0, (1,1))
ls_D2 =  (0, (3,1))
#ls_Dth = (0, (1, 1))

alpha_UT = 0.5
params_UT = {
            'color' : 'black', 
            'alpha' : alpha_UT,
            'linewidth' : 2.5,
            }
alpha_UMF = 0.7
params_UMF = {
            'color' : 'black', 
            'alpha' : alpha_UMF,
            'linewidth' : 2.5,
            }

alpha_MUMF = 1.0
params_MUMF = {
            'color' : 'black', 
            'alpha' : alpha_MUMF,
            'linewidth' : 1.5,
            }



ax.plot(t[1:], DErr_MUMF_25[1:], **params_MUMF, linestyle=ls_D2, label='$D[x_t - \hat{x}_t]$, $D[x_{T} - \hat{x}_T] = ' + "{:.2f}".format(DErr_MUMF_25[-1]) + '$ МУМНФ (25)')
ax.plot(t[1:], DErr_MUMF_100[1:], **params_MUMF, linestyle=ls_D1, label='$D[x_t - \hat{x}_t]$, $D[x_{T} - \hat{x}_T] = ' + "{:.2f}".format(DErr_MUMF_100[-1]) + '$ МУМНФ (100)')

ax.plot(t[1:], DErr_UMF[1:], **params_UMF, linestyle=ls_D, label='$D[x_t - \hat{x}_t]$, $D[x_{T} - \hat{x}_T] = ' + "{:.2f}".format(DErr_UMF[-1]) + '$ УМНФ')

ax.plot(t[1:], DErr_MUMF_1000[1:], **params_MUMF, linestyle=ls_D, label='$D[x_t - \hat{x}_t]$, $D[x_{T} - \hat{x}_T] = ' + "{:.2f}".format(DErr_MUMF_1000[-1]) + '$ МУМНФ (1000)')


ax.set_ylim(-0.01, 0.6)
ax.legend(loc=4)
plt.tight_layout()
plt.show()


