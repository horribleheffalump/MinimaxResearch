import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
import sys
import os
#from colormap import *
from matplotlib.gridspec import GridSpec
matplotlib.rc('font',**{'family':'serif'})
matplotlib.rc('text', usetex = True)

colormap = {
    'CMNF': 'red', 
    'BCMNF': 'red', 
    'UKF': 'blue', 
    'EKF': 'green', 
    'Dummy' : 'grey'}

folder = "D:/Наука/_Статьи/__в работе/2019 - CDC - УМНФ/results_CDC/test_DX0_2000_0.01/"
folder_EKF = "D:/Наука/_Статьи/__в работе/2019 - CDC - УМНФ/results_CDC/test_DX0_50_0.01/"
#folder = "D:/results/test_DX0_2000_0.01/"
#folder_EKF = "D:/results/test_DX0_50_0.01/"

fig = plt.figure(num=None, figsize=(9, 6), dpi=300)
#plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.1)
gs = GridSpec(2, 15, figure=fig)

ax_x = fig.add_subplot(gs[0, :6])
ax_y = fig.add_subplot(gs[0, 9:])

ax_legend = fig.add_subplot(gs[0, 6:9])

ax_v = fig.add_subplot(gs[1, :5])
ax_angle = fig.add_subplot(gs[1, 5:10])
ax_acc = fig.add_subplot(gs[1, 10:])

ax = [ax_x, ax_y, ax_v, ax_angle, ax_acc]
#l = ['x', 'y', 'V', 'angle', 'a']
lb = [-5, -5, -5, -0.5, -0.5]
ub = [80, 80, 140, 1.6, 3.5]
for i in range (0,5):
    inputfilename = folder + "TargetTracking_average_" + str(i) + ".txt"
    inputfilename_EKF = folder_EKF + "TargetTracking_average_" + str(i) + ".txt"
    outputfilename = folder + "TargetTracking_average_" + str(i) + ".png"



    
    data = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')
    data_EKF = pd.read_csv(inputfilename_EKF, delimiter = " ", dtype=float, engine='python')

    data[['EKF_mError']] = data_EKF[['EKF_mError']]
    data[['EKF_DError']] = data_EKF[['EKF_DError']]
    
    zero = 2
    T = int(data['t'].max())

    xlabels = ['a) Cartesian coordinate $x_t^1$', 'b) Cartesian coordinate $x_t^2$',  'c) Ground speed $x_t^3$' , 'd) Heading angle $x_t^4$', 'e) Normal acceleration $x_t^5$']

    filters = [s.replace("_mError","") for s in data.columns if "_mError" in s]
    #print(data.columns)
    print(filters)
    #ax.plot(data['t'][zero:T], np.sqrt(data['Dx'][zero:T]), c = 'k', linestyle=':', label = 'Dx')
    maxD = 0
    minM = 0

    filters = ['Dummy', 'EKF', 'UKF', 'BCMNF']
    for p in filters:
        if p == "BCMNF":
            ax[i].plot(data['t'][zero:T], np.sqrt(data[p+'_mKHat'][zero:T]), c=colormap[p],  linestyle='-', linewidth=4.0, alpha=0.3, label= 'theory ' + p)
        ax[i].plot(data['t'][zero:T], np.sqrt(data[p+'_DError'][zero:T] + np.square(data[p+'_mError'][zero:T])), c=colormap[p],  linestyle='-', linewidth=1.5, label= 's ' + p)
        if p != "EKF":
            maxD = max( maxD, np.max(np.sqrt(data[p+'_DError'][zero:T])))
            minM = min( minM, np.min(data[p+'_mError'][zero:T]))
            ax[i].plot(data['t'][zero:T], data[p+'_mError'][zero:T], c=colormap[p],  linestyle=':', linewidth=1.0, label = 'E ' + p)
        else:
            #ax[i].plot(data['t'][zero:T], data[p+'_mError'][zero:T], c=colormap[p],  linestyle=':', linewidth=1.0, label = 'E ' + p)
            ax[i].plot(data['t'][zero:42], data[p+'_mError'][zero:42], c=colormap[p],  linestyle=':', linewidth=1.0, label = 'E ' + p)


        ax[i].set_xlabel(xlabels[i])

    #ax[i].set_title(l[i])
    ax[i].set_ylim(lb[i], ub[i])
    #ax[i].set_ylim(minM * 1.1, maxD* 1.1)

for p in filters:
    pp = p
    if p == "Dummy": pp = "trivial"
    if p == "BCMNF": pp = "CMNF"

    ax_legend.plot([],[], c=colormap[p],  linestyle=':', linewidth=1.0, label = '$m^i$ bias ' + pp)
    ax_legend.plot([],[], c=colormap[p],  linestyle='-', linewidth=1.5, label= '$\sigma^i$ RMSE ' + pp)
    if p == "BCMNF": 
        ax_legend.plot([],[], c=colormap[p],  linestyle='-', linewidth=4.0, alpha=0.3, label= '$\sigma^i$ SD CMNF theory')

ax_legend.set_axis_off()
ax_legend.legend(loc=9, frameon=False)
plt.tight_layout()
plt.savefig("D:/Наука/_Статьи/__в работе/2019 - CDC - УМНФ/pic_stats.pdf")
#plt.show()



