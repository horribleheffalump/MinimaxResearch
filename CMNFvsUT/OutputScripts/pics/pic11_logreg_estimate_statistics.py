import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os
from multiplypoints import *

#inputfilename = os.path.join(sys.argv[1], sys.argv[2])
#outputfilename = os.path.join(sys.argv[1], sys.argv[3])
inputfilename = "D:/results/logreg-zero/LogisticModelZero_average_0.txt"

t, x, Dx, xhat_UMF, mErr_UMF, DErr_UMF, DErrTheor_UMF, xhat_UT, mErr_UT, DErr_UT, DErrTheor_UT = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True, dtype=float)

t, Dx, mErr_UMF, DErr_UMF, DErrTheor_UMF, mErr_UT, DErr_UT, DErrTheor_UT = multx(t), multy(Dx), multy(mErr_UMF), multy(DErr_UMF), multy(DErrTheor_UMF), multy(mErr_UT), multy(DErr_UT), multy(DErrTheor_UT)

from pylab import *

f = plt.figure(num=None, figsize=(5,5), dpi=200)
plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
ax = plt.subplot(111)

ls_m = (0, ())
ls_D =  (0, (5, 1))
ls_Dth = (0, (1, 1))

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

#ax.plot(t, mErr_UT, **params_UT, linestyle=ls_m)
#ax.plot(t, DErr_UT, **params_UT, linestyle=ls_D)
#ax.plot(t, DErrTheor_UT, **params_UT, linestyle =ls_Dth)

ax.plot(t, mErr_UMF, **params_UMF, linestyle=ls_m)
ax.plot(t, DErr_UMF, **params_UMF, linestyle=ls_D)
ax.plot(t, DErrTheor_UMF, color = 'black', alpha = alpha_UMF, linewidth = 1.5, linestyle = ls_Dth)

ax.fill_between(t, np.zeros_like(Dx), Dx, color='black', alpha = 0.2, linewidth=0.0);

ax.set_ylim(-10, 20000)

plt.tight_layout()
plt.show()
