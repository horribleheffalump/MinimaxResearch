import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
from scipy.stats import norm
from matplotlib import rc
rc('font',**{'family':'serif'})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble=r'\usepackage[T2A]{fontenc}')
rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex',preamble=r'\usepackage[russian]{babel}')
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

def textonly(ax, txt, fontsize = 10, loc = 2, *args, **kwargs):
    at = AnchoredText(txt,
                      prop=dict(size=fontsize), 
                      frameon=False,
                      loc=loc)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    return at

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


f = plt.figure(num=None, figsize=cm2inch((9,9)), dpi=200)

filename = u'D:/results/sampled_UMFs_good_predict/SampledRegression_bulk_0.txt'
data = pd.read_csv(filename, delimiter = " ", header=0, dtype=float, engine='python')
n, bins, patches = plt.hist(data['49'].values[:], 100, density=1, facecolor='black', alpha=0.5)


m = np.mean(data['49'])
s = np.std(data['49'])


at = textonly(plt.gca(), '$E[x_T] = ' +  "{:.2f}".format(m) +'$\n$D[x_T] = ' +  "{:.2f}".format(s*s)+'$', loc = 1)


#m = np.mean(data['49'])
#s = np.std(data['49'])
#min, max = np.min(data['49']), np.max(data['49'])
#print(m,s, min, max)
#x = np.arange(min, max, 0.1)
#_ = plt.plot(x, norm.pdf(x, m, s), color = 'black')

#print(data['49'])
#data[['49']].hist(bins=50)
#plt.hist(data[['49']],  density=True, facecolor='black', alpha=0.75, histtype='step')
plt.tight_layout()
plt.show()

