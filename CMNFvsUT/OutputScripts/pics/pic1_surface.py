import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

#mx, stdx = 30, 30
#my, stdy = 40, 30
mx, stdx = 300, 30
my, stdy = 400, 30

X = np.arange(mx - 2.0 * stdx, mx + 2.0 * stdx, (2.0)/ 10.0)
Y = np.arange(my - 2.0 * stdy, my + 2.0 * stdy, (2.0)/ 10.0)

#X = np.arange(-30.0, 90.0, (2.0)/ 10.0)
#Y = np.arange(-20.0, 100.0, (2.0)/ 10.0)
X, Y = np.meshgrid(X, Y)
Z = np.sqrt(X**2 + Y**2)
#print(X)
#print(Y)
#print(Z)

fig = plt.figure(figsize=(5, 5), dpi=200)
ax = fig.gca(projection='3d')
ax.view_init(15, 30)
surf = ax.plot_surface(X, Y, Z, cmap=cm.gray, linewidth=0, antialiased=False, alpha=0.3)
plt.tight_layout()
plt.show()
        