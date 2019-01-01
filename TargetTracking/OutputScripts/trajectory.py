import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import sys
import os
import glob

if (len(sys.argv)) > 1:
    folder = sys.argv[1]
else:
    folder = "D:/results/cont_EKF/"

colormap = {'CMNF': 'red', 'UKF': 'blue', 'MCMNF': 'green', 'RCMNF': 'orange', 'EKF': 'yellow'}

for file in glob.glob(folder + "TargetTracking_sample_*_0.txt"):
    if file.find('_obs') < 0 :
        print('processing ', file)
        outputfilename = file.replace('TargetTracking_sample', 'Trajectory').replace('_0.txt', '.png')

        data_stateX = pd.read_csv(file, delimiter = " ", dtype=float, engine='python')
    
        file = file.replace('_0.txt', '_1.txt')
        data_stateY = pd.read_csv(file, delimiter = " ", dtype=float, engine='python')


        f = plt.figure()
        ax = plt.subplot(111)

        plt.plot(data_stateX['x'], data_stateY['x'], 'k', linewidth = 1.0)

    
        protocols = [s.replace("_Error","") for s in data_stateX.columns if "_Error" in s]
        #protocols=['UKF']

        for p in protocols:
            ax.plot(data_stateX['x'] - data_stateX[p + '_Error'], data_stateY['x'] - data_stateY[p + '_Error'], c=colormap[p],  linestyle='-', linewidth=1.5, label=p)


        ax.legend()

        #observations
        #def pol2cart(rho, phi):
        #    x = rho * np.cos(phi)
        #    y = rho * np.sin(phi)
        #    return(x, y)

        #filenamePhi = folder + 'TargetTracking_sample_obs_0.txt'
        #data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", dtype=float, engine='python')

        #filenameRho = folder + 'TargetTracking_sample_obs_1.txt'
        #data_stateRho = pd.read_csv(filenameRho, delimiter = " ", dtype=float, engine='python')

        #X_R1 = [20000, 0];
        #Xobs, Yobs = pol2cart(data_stateRho['y'], data_statePhi['y'])
        #plt.scatter(Xobs+X_R1[0],Yobs+X_R1[1], c='orange')

        #filenamePhi = folder + 'TargetTracking_sample_obs_2.txt'
        #data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", dtype=float, engine='python')

        #filenameRho = folder + 'TargetTracking_sample_obs_3.txt'
        #data_stateRho = pd.read_csv(filenameRho, delimiter = " ", dtype=float, engine='python')

        #X_R2 = [0, 0];
        #Xobs, Yobs = pol2cart(data_stateRho['y'], data_statePhi['y'])
        #plt.scatter(Xobs+X_R2[0],Yobs+X_R2[1], c='cyan')

        #radars
        #plt.scatter(X_R1[0], X_R1[1], c='red', s=10)
        #plt.scatter(X_R2[0], X_R2[1], c='red', s=10)

        #plt.show()

        plt.savefig(outputfilename)
