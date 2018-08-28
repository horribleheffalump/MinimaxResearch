import numpy as np
import pandas as pd


for i in range(1,1000):
    filename = "D:\projects.git\MinimaxResearch\CMNFvsUT\output\InverseProportionGood_sample_0_" + str(i) + ".txt"
    data = pd.read_csv(filename, delimiter = " ", header=None, usecols=(0,1), dtype=float, names = ["t", "X"])
    max = np.max(data.X)
    if np.abs(max) > 10:
        print(i, max)
