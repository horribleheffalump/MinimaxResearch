import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm


filename = "C:\projects.git\MinimaxResearch\CMNFvsUT\output\samplereg-obs\SampledRegression_bulk_0.txt"
data = pd.read_csv(filename, delimiter = " ", header=0, dtype=float)
#data[['0', '10', '20', '30', '40', '49']].hist(bins=20)
#data[['49']].hist(bins=50)

n, bins, patches = plt.hist(data['49'].values[:], 50, density=1, facecolor='black', alpha=0.5)

m = np.mean(data['49'])
s = np.std(data['49'])
min, max = np.min(data['49']), np.max(data['49'])
print(m,s, min, max)




x = np.arange(min, max, 0.1)
_ = plt.plot(x, norm.pdf(x, m, s), color = 'black')


plt.show()
