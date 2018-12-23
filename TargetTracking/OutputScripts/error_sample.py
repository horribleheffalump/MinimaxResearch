import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
import sys
import os
#from multiplypoints import *

#folder = "D:/results/cont_EKF/"
folder = sys.argv[1]

colormap = {'CMNF': 'red', 'UKF': 'blue', 'MCMNF': 'green', 'RCMNF': 'orange', 'EKF': 'yellow'}
folder = "D:/results/cont_EKF/"
for i in range (0,5):
    inputfilename = folder + "TargetTracking_sample_"+str(i)+".txt"
    outputfilename = folder + "TargetTracking_errorsample_" + str(i) + ".png"

    data = pd.read_csv(inputfilename, delimiter = " ", dtype=float, engine='python')
    from pylab import *

    f = plt.figure()
    ax1 = plt.subplot(111)
    
    plotcols = [s for s in data.columns if "_Error" in s]

    #ax1.plot(data['t'], data['x'], 'k', label='x')

    for p in plotcols:
        filtname = p.replace("_Error","")
        ax1.plot(data['t'], data[p], c=colormap[filtname], label=filtname)

    ax1.legend()

    #plt.show()    
    plt.savefig(outputfilename)

