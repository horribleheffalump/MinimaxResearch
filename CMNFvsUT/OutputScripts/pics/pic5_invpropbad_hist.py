import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd



fig = plt.figure(figsize=(5, 5), dpi=200)

filename = u'D:/results/invprop_bad/InverseProportionBad_bulk_0.txt'
data = pd.read_csv(filename, delimiter = " ", header=0, dtype=float, engine='python')
n, bins, patches = plt.hist(data['49'].values[:], 40, density=1, facecolor='black', alpha=0.5)
#print(data['49'])
#data[['49']].hist(bins=50)
#plt.hist(data[['49']],  density=True, facecolor='black', alpha=0.75, histtype='step')
plt.tight_layout()
plt.show()

