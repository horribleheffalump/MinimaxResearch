import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

X = np.arange(-30.0, 90.0, (2.0)/ 100.0)
Y = np.arange(-20.0, 100.0, (2.0)/ 100.0)
X, Y = np.meshgrid(X, Y)
Z = np.sqrt(X**2 + Y**2)
#print(X)
#print(Y)
#print(Z)

fig = plt.figure(figsize=(10, 6), dpi=200)
ax = fig.gca(projection='3d')
ax.view_init(50, 30)
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False, alpha=0.3)
plt.show()
        