import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)


matplotlib.rc('text', usetex = True)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

filename = u"../output/optimize_UKF.txt"
crit, alpha, beta, kappa = np.loadtxt(filename, delimiter = ' ', usecols=(0,1,2,3), unpack=True, dtype=float)


ax.scatter(alpha, beta, kappa, c=crit, cmap=plt.hot())

ax.set_xlabel(r"$\alpha$")
ax.set_ylabel(r"$\beta$")
ax.set_zlabel(r"$\kappa$")

plt.show()