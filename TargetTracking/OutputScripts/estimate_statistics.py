import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
import sys
import os
from colormap import *

if (len(sys.argv)) > 1:
    folder = sys.argv[1]
else:
    folder = "D:/results/cubic/"

#inputfilename = sys.argv[1] 
#outputfilename = sys.argv[2]

l = ['x', 'y', 'V', 'angle', 'a']
lb = [-5, -5, -5, -0.5, 0.5]
ub = [100, 100, 5, 0.5, 5.0]
for i in range (0,5):
    inputfilename = folder + "TargetTracking_average_" + str(i) + ".txt"
    outputfilename = folder + "TargetTracking_average_" + str(i) + ".png"



    
    data = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    f = plt.figure()
    ax = plt.subplot(111)

    T = int(data['t'].max())


    protocols = [s.replace("_mError","") for s in data.columns if "_mError" in s]
    #print(data.columns)
    print(protocols)
    #ax.plot(data['t'][1:T], np.sqrt(data['Dx'][1:T]), c = 'k', linestyle=':', label = 'Dx')
    maxD = 0
    minM = 0
    #if "CMNF" in protocols:        
    #    protocols.remove("CMNF")
    #protocols.remove("EKF")
    #protocols.remove("UKFOptRandomStepwise")
    #protocols.remove("BCMNF")
    #protocols.remove("Dummy")
    #protocols.remove("UKF")

    for p in protocols:
        #ax.plot(data['t'][1:T], data[p+'_mError'][1:T], c=colormap[p],  linestyle=':', linewidth=1.0, label = 'E ' + p)
        ax.plot(data['t'][1:T], np.sqrt(data[p+'_DError'][1:T]), c=colormap[p],  linestyle='-', linewidth=1.5, label= 's ' + p)
        #ax.plot(data['t'][1:T], np.sqrt(data[p+'_mKHat'][1:T]), c=colormap[p],  linestyle='--', linewidth=2.5, label= 'theory ' + p)
        if p != "EKF":
            maxD = max( maxD, np.max(np.sqrt(data[p+'_DError'][1:T])))
            minM = min( minM, np.min(data[p+'_mError'][1:T]))

    ax.set_title(l[i])
    #plt.show()
    ax.legend()
    #ax.set_ylim(lb[i], ub[i])
    ax.set_ylim(minM * 1.1, maxD* 1.1)
    plt.savefig(outputfilename)



