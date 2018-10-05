import numpy as np
import pandas as pd

#many files
for i in range(1,1000):
    filename = "D:/results/invprop_bad/InverseProportionBad_sample_" + str(i) + "_0.txt"
    data = pd.read_csv(filename, delimiter = " ", header=None, usecols=(0,1), dtype=float, names = ["t", "X"])
    max = np.max(data.X)
    if np.abs(max) > 500:
        print(i, max)

#single bulk file
#filename = u"D:/results/invprop_good/InverseProportionGood_bulk_0.txt"
#data = pd.read_csv(filename, delimiter = " ", header=0, dtype=float, engine='python')
#print(data.shape)
#for i in range(0,data.shape[0]):
#    max = np.max(np.abs(data.iloc[i]))
#    if np.abs(max) > 500:
#        print(i, max)
