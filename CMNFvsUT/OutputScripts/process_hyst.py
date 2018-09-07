import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = "D:\projects.git\MinimaxResearch\CMNFvsUT\output\LogisticModelZero_bulk_0.txt"
data = pd.read_csv(filename, delimiter = " ", header=0, dtype=float)
data[['0', '10', '20', '30', '40', '49']].hist(bins=20)
plt.show()
