import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
from scipy.stats import norm


fig = plt.figure(figsize=(5, 5), dpi=200)

filename = u'D:/results/cubic/CubicSensor_bulk_0.txt'
data = pd.read_csv(filename, delimiter = " ", header=0, dtype=float, engine='python')
n, bins, patches = plt.hist(data['49'].values[:], 40, density=1, facecolor='black', alpha=0.5)

m = np.mean(data['49'])
s = np.std(data['49'])
min, max = np.min(data['49']), np.max(data['49'])
print(m,s, min, max)
x = np.arange(min, max, 0.1)
_ = plt.plot(x, norm.pdf(x), color = 'black')

#print(data['49'])
#data[['49']].hist(bins=50)
#plt.hist(data[['49']],  density=True, facecolor='black', alpha=0.75, histtype='step')
plt.tight_layout()
plt.show()

