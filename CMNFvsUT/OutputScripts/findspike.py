import numpy as np
import pandas as pd

#many files
for i in range(1,1000):
    filename = "D:/results/invprop_bad/InverseProportionBad_sample_" + str(i) + "_0.txt"
    data = pd.read_csv(filename, delimiter = " ", header=None, usecols=(0,1), dtype=float, names = ["t", "X"])
    max = np.max(data.X)
    if np.abs(max) > 500:
        print(i, max)

#15 94096.9641330407
#67 3771.8394108797797
#262 13127.529388448402
#358 2424.2835335792
#429 17611.5877732941
#459 1503.4417638628101
#461 558.225180332397
#540 1372.1550713291902
#548 1049.362068331
#738 1268.0231331540801
#741 874.663652261001
#809 702.9439525988389
#960 100018.185615116
#986 99976.4333353638
#990 677.676569179953


#single bulk file
#filename = u"D:/results/invprop_good/InverseProportionGood_bulk_0.txt"
#data = pd.read_csv(filename, delimiter = " ", header=0, dtype=float, engine='python')
#print(data.shape)
#for i in range(0,data.shape[0]):
#    max = np.max(np.abs(data.iloc[i]))
#    if np.abs(max) > 500:
#        print(i, max)
