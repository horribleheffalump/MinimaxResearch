import pandas as pd
import matplotlib.pyplot as plt
folder = 'D:/results/polar_0_400/'
#folder = 'D:/results/polar_0_-400/'
#folder = 'D:/results/polar_400_0/'
#folder = 'D:/results/polar_-400_0/'
filename = folder + 'test_polar_alldata.txt'
df = pd.read_csv(filename, sep=';', decimal=',', header=0)
df.columns = df.columns.str.strip()
dsc = df[['Y_0']].describe()
#hist = df[['Y_0']].hist(bins = 500)
#scatter = plt.scatter(df[['Xinv_0']], df[['Xinv_1']], c='blue')
#scatter = plt.scatter(df[['X_0']], df[['X_1']], c='red')
hist = df[['Xinv_1']].hist(bins = 500)
print(dsc)
plt.savefig(folder+'Yinv_hist.png')
#plt.savefig(folder+'scatter.png')

#c = 0
#for i in range(100,10001):
#    m = 0
#    for k in range(0,10):
#        a = str(i).count(str(k))
#        if a>m:
#            m=a
#    if m == 3:
#        c = c+1
#    #print(i, m)
#print(c)

#(1/2/3)/4/5/(6/7/8/9)
#(1/2/3/4)/(5/6/7/8/9)
#(1/2)/(3/4/5/6/7/8/9)
#1/2/3/(4/5/6/7/8/9)
#1/(2/3/4/5/6/7/8/9)

