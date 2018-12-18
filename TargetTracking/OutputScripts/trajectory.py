import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd


#filename = u'D:/results/cont/state.txt'
#data_state = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')
#filename = u'D:/results/cont/obs.txt'
#data_obs = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')

#plt.plot(data_state[1], data_state[2])
##plt.scatter(data_obs[0], data_obs[1])


#filename = u'D:/results/cont/state1.txt'
#data_state = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')
#plt.plot(data_state[1], data_state[2], 'r')

#filename = u'D:/results/cont/state2.txt'
#data_state = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')
#plt.plot(data_state[1], data_state[2], 'g')

#filename = u'D:/results/cont/state3.txt'
#data_state = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')
#plt.plot(data_state[1], data_state[2], 'b')

#filename = u'D:/results/cont/state4.txt'
#data_state = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')
#plt.plot(data_state[1], data_state[2], 'k')

filenameX = u'D:/results/cont_a_1.0_two/TargetTracking_sample_0.txt'
data_stateX = pd.read_csv(filenameX, delimiter = " ", header=None, dtype=float, engine='python')

filenameY = u'D:/results/cont_a_1.0_two/TargetTracking_sample_1.txt'
data_stateY = pd.read_csv(filenameY, delimiter = " ", header=None, dtype=float, engine='python')

plt.plot(data_stateX[1], data_stateY[1], 'k', linewidth = 3.0)
#plt.plot(data_stateX[1]-data_stateX[2], data_stateY[1]-data_stateY[2], 'r')
#plt.plot(data_stateX[1]-data_stateX[3], data_stateY[1]-data_stateY[3], 'g')
#plt.plot(data_stateX[1]-data_stateX[4], data_stateY[1]-data_stateY[4], 'b')

plt.scatter(data_stateX[1]-data_stateX[2], data_stateY[1]-data_stateY[2], c = 'r', linewidth=2.0)
plt.scatter(data_stateX[1]-data_stateX[3], data_stateY[1]-data_stateY[3], c = 'g', linewidth=2.0)
plt.scatter(data_stateX[1]-data_stateX[4], data_stateY[1]-data_stateY[4], c = 'b')
plt.scatter(data_stateX[1]-data_stateX[5], data_stateY[1]-data_stateY[5], c = 'orange')


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

filenamePhi = u'D:/results/cont_a_1.0_two/TargetTracking_sample_obs_0.txt'
data_statePhi = pd.read_csv(filenamePhi, delimiter = " ", header=None, dtype=float, engine='python')

filenameRho = u'D:/results/cont_a_1.0_two/TargetTracking_sample_obs_1.txt'
data_stateRho = pd.read_csv(filenameRho, delimiter = " ", header=None, dtype=float, engine='python')

#X_R = [20000, 0];
#Xobs, Yobs = pol2cart(data_stateRho[1], data_statePhi[1])
#plt.scatter(Xobs+X_R[0],Yobs+X_R[1], c='k')

#plt.show() 


plt.savefig("D:/results/cont_a_1.0_two/Trajectory.png")