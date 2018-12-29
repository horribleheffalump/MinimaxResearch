import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import sys
import os

folder = "D:/results/cont_EKF_divergance_investigation/"
#folder = sys.argv[1]

colormap = {'CMNF': 'red', 'UKF': 'blue', 'MCMNF': 'green', 'RCMNF': 'orange', 'EKF': 'yellow'}

outputfilename = folder + "Trajectory.png"

inputfilename = folder + "TargetTracking_sample_611_" + str(0) + ".txt"
data_stateX = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

inputfilename = folder + "TargetTracking_sample_611_" + str(1) + ".txt"
data_stateY = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')


f = plt.figure()
ax = plt.subplot(111)

plt.plot(data_stateX['x'], data_stateY['x'], 'k', linewidth = 1.0)

    
protocols = [s.replace("_Error","") for s in data_stateX.columns if "_Error" in s]
#protocols=['UKF']

for p in protocols:
    ax.plot(data_stateX['x'] - data_stateX[p + '_Error'], data_stateY['x'] - data_stateY[p + '_Error'], c=colormap[p],  linestyle='-', linewidth=1.5, label=p)


ax.legend()


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

filenamePhi = u'D:/results/cont_EKF_divergance_investigation/TargetTracking_sample_obs_611_0.txt'
data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", dtype=float, engine='python')

filenameRho = u'D:/results/cont_EKF_divergance_investigation/TargetTracking_sample_obs_611_1.txt'
data_stateRho = pd.read_csv(filenameRho, delimiter = " ", dtype=float, engine='python')

X_R1 = [20000, 0];
Xobs, Yobs = pol2cart(data_stateRho['y'], data_statePhi['y'])
plt.scatter(Xobs+X_R1[0],Yobs+X_R1[1], c='orange')

filenamePhi = u'D:/results/cont_EKF_divergance_investigation/TargetTracking_sample_obs_611_2.txt'
data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", dtype=float, engine='python')

filenameRho = u'D:/results/cont_EKF_divergance_investigation/TargetTracking_sample_obs_611_3.txt'
data_stateRho = pd.read_csv(filenameRho, delimiter = " ", dtype=float, engine='python')

X_R2 = [0, 0];
Xobs, Yobs = pol2cart(data_stateRho['y'], data_statePhi['y'])
plt.scatter(Xobs+X_R2[0],Yobs+X_R2[1], c='cyan')

plt.scatter(X_R1[0], X_R1[1], c='red', s=10)
plt.scatter(X_R2[0], X_R2[1], c='red', s=10)
plt.show()

plt.savefig(outputfilename)
