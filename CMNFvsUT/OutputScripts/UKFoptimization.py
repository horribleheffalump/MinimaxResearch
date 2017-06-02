import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)


matplotlib.rc('text', usetex = True)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

filename = u"../output/optimize_UKF_filtered_10best.txt"
crit, alpha, beta, kappa = np.loadtxt(filename, delimiter = ' ', usecols=(0,1,2,3), unpack=True, dtype=float)
L = 2
l = alpha*alpha*(L + kappa) - L
a = pow(l + L, 0.5)
W0m = l / (L + l)
W0c = W0m + 1 - alpha*alpha + beta
Wi = 1 / 2 / (L + l)
#ax.scatter(alpha, beta, kappa, c=crit, cmap=plt.hot())

#ax.set_xlabel(r"$\alpha$")
#ax.set_ylabel(r"$\beta$")
#ax.set_zlabel(r"$\kappa$")
#ax.scatter(W0m, W0c, Wi, c=crit, cmap=plt.hot())


ax.scatter(a, W0m, W0c, c=crit, cmap=plt.hot())
ax.set_xlabel(r"$a$")
ax.set_ylabel(r"$W_0^{(c)}$")
ax.set_zlabel(r"$W_0^{(m)}$")




#ax.scatter(W0m, W0c, Wi, c=crit, cmap=plt.hot())
#ax.set_xlabel(r"$W_0^{(c)}$")
#ax.set_ylabel(r"$W_0^{(m)}$")
#ax.set_zlabel(r"$W_i$")


ax.set_xlim(4,5)
ax.set_ylim(0,1)
ax.set_zlim(0,5)

plt.show()