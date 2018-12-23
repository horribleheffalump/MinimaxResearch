import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import sys
import os

#folder = "D:/results/cont_EKF/"
folder = sys.argv[1]

colormap = {'CMNF': 'red', 'UKF': 'blue', 'MCMNF': 'green', 'RCMNF': 'orange', 'EKF': 'yellow'}

outputfilename = folder + "Trajectory.png"

inputfilename = folder + "TargetTracking_sample_" + str(0) + ".txt"
data_stateX = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

inputfilename = folder + "TargetTracking_sample_" + str(1) + ".txt"
data_stateY = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')


f = plt.figure()
ax = plt.subplot(111)

plt.plot(data_stateX['x'], data_stateY['x'], 'k', linewidth = 1.0)

    
protocols = [s.replace("_Error","") for s in data_stateX.columns if "_Error" in s]

for p in protocols:
    ax.plot(data_stateX['x'] - data_stateX[p + '_Error'], data_stateY['x'] - data_stateY[p + '_Error'], c=colormap[p],  linestyle='-', linewidth=1.5, label=p)


ax.legend()


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

#filenamePhi = u'D:/results/cont_million/TargetTracking_sample_obs_0.txt'
#data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", header=None,
#dtype=float, engine='python')

#filenameRho = u'D:/results/cont_million/TargetTracking_sample_obs_1.txt'
#data_stateRho = pd.read_csv(filenameRho, delimiter = " ", header=None,
#dtype=float, engine='python')

#X_R = [20000, 0];
#Xobs, Yobs = pol2cart(data_stateRho[1], data_statePhi[1])
#plt.scatter(Xobs+X_R[0],Yobs+X_R[1], c='k')

#plt.show()

plt.savefig(outputfilename)
