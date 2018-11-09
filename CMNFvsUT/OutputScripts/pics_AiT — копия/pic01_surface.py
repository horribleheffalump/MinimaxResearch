import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif'})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble=r'\usepackage[T2A]{fontenc}')
rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex',preamble=r'\usepackage[russian]{babel}')

mx, stdx = 30, 30
my, stdy = 40, 30
#mx, stdx = 300, 30
#my, stdy = 400, 30

#X = np.arange(mx - 2.0 * stdx, mx + 2.0 * stdx, (2.0)/ 10.0)
#Y = np.arange(my - 2.0 * stdy, my + 2.0 * stdy, (2.0)/ 10.0)
X = np.arange(-25, 125, (2.0)/ 10.0)
Y = np.arange(-25, 125, (2.0)/ 10.0)

X, Y = np.meshgrid(X, Y)
Z = np.sqrt(X**2 + Y**2)

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


fig = plt.figure(num=None, figsize=cm2inch((9,9)), dpi=200)
ax = fig.gca(projection='3d')
ax.view_init(15, 30)
surf = ax.plot_surface(X, Y, Z, cmap=cm.gray, linewidth=0, antialiased=False, alpha=0.3)
#ax.set_xticks([250,275,300,325,350])
#ax.set_xticklabels(['250','','300','','350'])
#ax.set_yticks([350,375,400,425,450])
#ax.set_yticklabels(['350','','400','','450'])
#ax.set_zticks([450,475,500,525,550])
#ax.set_zticklabels(['450','','500','','550'])

ax.set_xticks([0,25,50,75,100])
ax.set_xticklabels(['0','','50','','100'])
ax.set_yticks([0,25,50,75,100])
ax.set_yticklabels(['0','','50','','100'])
ax.set_zticks([0,25,50,75,100,125,150])
ax.set_zticklabels(['0','','50','','100','','150'])


plt.tight_layout()
plt.show()
        