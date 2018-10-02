import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os
from multiplypoints import *

#inputfilename = os.path.join(sys.argv[1], sys.argv[2])
#outputfilename = os.path.join(sys.argv[1], sys.argv[3])
inputfilename = "D:/results/logreg-zero/LogisticModelZero_sample_0.txt"

t, x, y = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2), unpack=True, dtype=float)

t, x, y = multx(t), multy(x), multy(y)

from pylab import *

f = plt.figure(num=None, figsize=(5,5), dpi=200)
plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
ax = plt.subplot(111)

ls_x = (0, ()) #solid
ls_y =  (0, (1, 1)) #dots

params = {
            'color' : 'black', 
            'linewidth' : 1.5,
            }

ax.plot(t, x, **params, linestyle=ls_x)
ax.plot(t, y, **params, linestyle=ls_y)

ax.set_ylim(-1000.0, 1000.0)

plt.tight_layout()
plt.show()



