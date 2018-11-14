import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd


filename = u'D:/results/cont/state.txt'
data_state = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')
filename = u'D:/results/cont/obs.txt'
data_obs = pd.read_csv(filename, delimiter = " ", header=None, dtype=float, engine='python')

plt.plot(data_state[0], data_state[1])
plt.scatter(data_obs[0], data_obs[1])
plt.show() 