import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import sys
import os
import glob
from colormap import *
matplotlib.rc('font',**{'family':'serif'})
matplotlib.rc('text', usetex = True)


colormap = {
    'CMNF': 'red', 
    'BCMNF': 'magenta', 
    'UKF': 'blue', 
    'EKF': 'green', 
    'Dummy' : 'grey'}


if (len(sys.argv)) > 1:
    folder = sys.argv[1]
else:
    folder = "D:/results/test_DX0_2000_0.01/"

for file in glob.glob(folder + "TargetTracking_sample*_0.txt"):
    if file.find('_obs') < 0 :
        print('processing ', file)
        outputfilename = file.replace('TargetTracking_sample', 'Trajectory').replace('_0.txt', '.pdf')

        data_stateX = pd.read_csv(file, delimiter = " ", dtype=float, engine='python')
    
        file = file.replace('_0.txt', '_1.txt')
        data_stateY = pd.read_csv(file, delimiter = " ", dtype=float, engine='python')


        fig = plt.figure(num=None, figsize=(4, 3), dpi=300)
        ax = plt.subplot(111)

        plt.plot(data_stateX['x']/1000, data_stateY['x']/1000, 'k', linewidth = 1.0, label = 'Sample path')

    
        #protocols = [s.replace("_Error","") for s in data_stateX.columns if "_Error" in s]
        #protocols.remove('EKF')

        #for p in protocols:
        #    ax.plot(data_stateX['x'] - data_stateX[p + '_Error'], data_stateY['x'] - data_stateY[p + '_Error'], c=colormap[p],  linestyle='-', linewidth=1.5, label=p)


        #ax.legend()

        #observations
        #def pol2cart(rho, phi):
        #    x = rho * np.cos(phi)
        #    y = rho * np.sin(phi)
        #    return(x, y)

        #filenamePhi = folder + 'TargetTracking_sample_obs_0.txt'
        #data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", dtype=float, engine='python')

        #filenameRho = folder + 'TargetTracking_sample_obs_1.txt'
        #data_stateRho = pd.read_csv(filenameRho, delimiter = " ", dtype=float, engine='python')

        X_R1 = [-10000, 10000];
        #Xobs, Yobs = pol2cart(data_stateRho['y'], data_statePhi['y'])
        #plt.scatter(Xobs+X_R1[0],Yobs+X_R1[1], c='orange')

        #filenamePhi = folder + 'TargetTracking_sample_obs_2.txt'
        #data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", dtype=float, engine='python')

        #filenameRho = folder + 'TargetTracking_sample_obs_3.txt'
        #data_stateRho = pd.read_csv(filenameRho, delimiter = " ", dtype=float, engine='python')

        X_R2 = [0, 0];
        #Xobs, Yobs = pol2cart(data_stateRho['y'], data_statePhi['y'])
        #plt.scatter(Xobs+X_R2[0],Yobs+X_R2[1], c='cyan')

        #radars
        plt.scatter(X_R1[0]/1000, X_R1[1]/1000, c='red', s=10, label = 'Observers')
        plt.scatter(X_R2[0]/1000, X_R2[1]/1000, c='red', s=10)

        #plt.show()

        start = (0/1000,25000/1000)
        sigma=2000/1000

        plt.plot([],[], color='blue', alpha = 0.1, linewidth=4.0, label ='Starting area')
        c = plt.Circle(start, 3*sigma, color='blue', alpha = 0.1)
        ax.add_artist(c)
        ax.axis('equal')

        ax.set_xlabel('$x_t^1$ [km]')
        ax.set_ylabel('$x_t^2$ [km]')
        
        plt.legend()
        plt.tight_layout()
        plt.savefig(outputfilename)

