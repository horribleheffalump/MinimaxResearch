import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import sys
import os

folder = "D:/results/cont_EKF_divergance_investigation/"
#folder = sys.argv[1]

colormap = {'CMNF': 'red', 'UKF': 'blue', 'MCMNF': 'green', 'RCMNF': 'orange', 'EKF': 'yellow'}
for i in range(0, 1000):

    outputfilename = folder + "Trajectory" + str(i) + ".png"

    inputfilename = folder + "TargetTracking_sample_" + str(i) + "_" + str(0) + ".txt"
    data_stateX = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    inputfilename = folder + "TargetTracking_sample_" + str(i) + "_" + str(1) + ".txt"
    data_stateY = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    mx = data_stateX['x'].min()
    my = data_stateY['x'].min()
    if (mx < 0): 
        print(i, "x", mx)
    if (my < 0): 
        print(i, "y", my)

    inputfilename = folder + "TargetTracking_sample_obs_" + str(i) + "_" + str(0) + ".txt"
    data_obs1 = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    inputfilename = folder + "TargetTracking_sample_obs_" + str(i) + "_" + str(2) + ".txt"
    data_obs2 = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    min_angle1 = data_obs1['y'].min()
    min_angle2 = data_obs2['y'].min()
    if (min_angle1 < 0): 
        print(i, "1angle", min_angle1)
    if (min_angle2 < 0): 
        print(i, "2angle", min_angle2)

    inputfilename = folder + "TargetTracking_sample_obs_" + str(i) + "_" + str(1) + ".txt"
    data_obs1 = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    inputfilename = folder + "TargetTracking_sample_obs_" + str(i) + "_" + str(3) + ".txt"
    data_obs2 = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    min_r1 = data_obs1['y'].min()
    min_r2 = data_obs2['y'].min()
    if (min_r1 < 20000): 
        print(i, "1r", min_r1)
    if (min_r2 < 20000): 
        print(i, "2r", min_r2)


    #f = plt.figure()
    #ax = plt.subplot(111)

    #plt.plot(data_stateX['x'], data_stateY['x'], 'k', linewidth = 1.0)

    
    #protocols = [s.replace("_Error","") for s in data_stateX.columns if "_Error" in s]

    #for p in protocols:
    #    ax.plot(data_stateX['x'] - data_stateX[p + '_Error'], data_stateY['x'] - data_stateY[p + '_Error'], c=colormap[p],  linestyle='-', linewidth=1.5, label=p)


    #ax.legend()

    #plt.savefig(outputfilename)

