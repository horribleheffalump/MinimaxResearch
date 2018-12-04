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

inputfilename = "D:/results/cont/TargetTracking_average_0.txt"
outputfilename = "D:/results/cont/TargetTracking_average_0.pdf"




t, x, Dx, xhat_UMF, mErr_UMF, DErr_UMF, DErrTheor_UMF, xhat_UT, mErr_UT, DErr_UT, DErrTheor_UT = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True, dtype=float)

max_DErr_UMF = max(sorted(DErr_UMF)[0 : round(len(t) * 0.98)])
max_DErr_UT = max(sorted(DErr_UT)[0 : round(len(t) * 0.98)])

maxDErr = max(max_DErr_UMF, max_DErr_UT)

min_mErr_UMF = min(sorted(mErr_UMF)[round(len(t) * 0.02) : len(t)])
min_mErr_UT = min(sorted(mErr_UT)[round(len(t) * 0.02) : len(t)])

minmErr = min(min_mErr_UMF, min_mErr_UT)

from pylab import *

f = plt.figure(num=None, figsize=(14, 3.5), dpi=150, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])     
gs.update(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1, hspace=0.1)
ax_CMNF = plt.subplot(gs[0])
ax_UKF = plt.subplot(gs[1])

ls_m = (0, ())
ls_D =  (0, (5, 1))
ls_Dth = (0, (1, 1))

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

ax_UKF.plot(t, mErr_UT, **params_UT, linestyle=ls_m)
ax_UKF.plot(t, DErr_UT, **params_UT, linestyle=ls_D)
#ax_UKF.plot(t, DErrTheor_UT, **params_UT, linestyle =ls_Dth)

ax_CMNF.plot(t, mErr_UMF, **params_UMF, linestyle=ls_m)
ax_CMNF.plot(t, DErr_UMF, **params_UMF, linestyle=ls_D)
#ax_CMNF.plot(t, DErrTheor_UMF, color = 'black', alpha = alpha_UMF, linewidth = 1.5, linestyle = ls_Dth)

#ax1.set_ylim(minmErr, maxDErr * 1.1)


#ax_UKF.set_ylim(0,100)
#ax_CMNF.set_ylim(0,100)
#plt.show()
plt.savefig(outputfilename)



