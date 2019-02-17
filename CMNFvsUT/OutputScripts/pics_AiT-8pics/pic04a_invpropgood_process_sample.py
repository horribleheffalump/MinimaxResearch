import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc('text', usetex = True)
import pylab
import sys
import os
from multiplypoints import *

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

#inputfilename = os.path.join(sys.argv[1], sys.argv[2])
#outputfilename = os.path.join(sys.argv[1], sys.argv[3])
inputfilename = "D:/results//invprop_good/InverseProportionGood_sample_2982_0.txt"

t, x, y = np.loadtxt(inputfilename, delimiter = ' ', usecols=(0,1,2), unpack=True, dtype=float)

t, x, y = multx(t), multy(x), multy(y)

from pylab import *

f = plt.figure(num=None, figsize=cm2inch((12,9)), dpi=200)
plt.subplots_adjust(left=0.06, bottom=0.07, right=0.98, top=0.95, wspace=0.1)
ax = plt.subplot(111)

ls_x = (0, ()) #solid
ls_y =  (0, (1, 1)) #dots

params = {
            'color' : 'black', 
            'linewidth' : 1.5,
            }

ax.plot(t, x, **params, linestyle=ls_x, label='$x_t$')
ax.plot(t, y, **params, linestyle=ls_y, label='$y_t$')

#ax.set_ylim(-10.0, 10.0)
#ax.legend(loc=1)
plt.tight_layout()

#float commas
rc('text.latex',preamble=r'\usepackage{icomma}')
import matplotlib as mpl
import locale
locale.setlocale(locale.LC_ALL, "rus_RUS")
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useLocale=True, useMathText=True))
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useLocale=True, useMathText=True))

plt.show()



