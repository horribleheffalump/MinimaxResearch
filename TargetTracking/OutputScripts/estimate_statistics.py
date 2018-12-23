import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
import sys
import os

#folder = "D:/results/cont_EKF/"
folder = sys.argv[1]
#inputfilename = sys.argv[1] 
#outputfilename = sys.argv[2]
colormap = {'CMNF': 'red', 'UKF': 'blue', 'MCMNF': 'green', 'RCMNF': 'orange', 'EKF': 'yellow'}
l = ['x', 'y', 'V', 'angle', 'a']
for i in range (0,5):
    inputfilename = folder + "TargetTracking_average_" + str(i) + ".txt"
    outputfilename = folder + "TargetTracking_average_" + str(i) + ".png"



    
    data = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')

    f = plt.figure()
    ax = plt.subplot(111)

    T = 150

    
    protocols = [s.replace("_mError","") for s in data.columns if "_mError" in s]

    ax.plot(data['t'][1:T], np.sqrt(data['Dx'][1:T]), c = 'k', linestyle=':', label = 'Dx')

    for p in protocols:
        ax.plot(data['t'][1:T], data[p+'_mError'][1:T], c=colormap[p],  linestyle=':', linewidth=1.0, label = 'E ' + p)
        ax.plot(data['t'][1:T], np.sqrt(data[p+'_DError'][1:T]), c=colormap[p],  linestyle='-', linewidth=1.5, label= 's ' + p)

    ax.set_title(l[i])
    #plt.show()
    ax.legend()
    plt.savefig(outputfilename)



