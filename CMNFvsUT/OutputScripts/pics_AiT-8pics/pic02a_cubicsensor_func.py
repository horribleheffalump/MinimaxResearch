
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif'})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble=r'\usepackage[T2A]{fontenc}')
rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex',preamble=r'\usepackage[russian]{babel}')



X = np.arange(-50.0, 50.0, 0.1)
Y = X / (1 + X*X)

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


f = plt.figure(num=None, figsize=cm2inch((9,9)), dpi=200)
_ = plt.plot(X, Y, linewidth=1.0, color='black')
plt.tight_layout()
ax = plt.gca()

#float commas
rc('text.latex',preamble=r'\usepackage{icomma}')
import matplotlib as mpl
import locale
locale.setlocale(locale.LC_ALL, "rus_RUS")
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useLocale=True, useMathText=True))
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useLocale=True, useMathText=True))

plt.show()
        

