import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os
import pandas as pd

inputfilename = sys.argv[1] 
outputfilename = sys.argv[2]

inputfilename = "D:/results/logreg-zero_mcmnf_revised/LogisticModelZero_average_0.txt"
outputfilename = "D:/results/logreg-zero_mcmnf_revised/LogisticModelZero_estimate_statistics_single_0.pdf"

#xhat_UMF, mErr_UMF, DErr_UMF, DErrTheor_UMF, xhat_UT, mErr_UT, DErr_UT, DErrTheor_UT = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True, dtype=float)

data = pd.read_csv(inputfilename, delimiter = " ", header=None, dtype=float)
n = int((data.shape[1] - 3) / 4)

f = plt.figure(num=None, figsize=(14, 3.5), dpi=150, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
ax = plt.subplot(111)

ls_m = (0, ())
ls_D =  (0, (5, 1))
ls_Dth = (0, (1, 1))

colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow']
for j in range(n):
    #ax.plot(data[[0]], data[[3+j*4+1]], linestyle=ls_m, color=colors[j], linewidth=2.5, alpha=0.7)
    ax.plot(data[[0]], data[[3+j*4+2]], linestyle=ls_D, color=colors[j], linewidth=2.5, alpha=0.7)
    #ax.plot(data[[0]], data[[3+j*4+3]], linestyle=ls_Dth, color=colors[j], linewidth=2.5, alpha=0.7)


plt.savefig(outputfilename)


