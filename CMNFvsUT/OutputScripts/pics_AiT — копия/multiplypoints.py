import numpy as np
def multx(x):        
    xplot = np.zeros(x.size * 2 - 1)
    for i in range(0, x.size - 1):
        xplot[i * 2] = x[i]
        xplot[i * 2 + 1] = x[i + 1]
    xplot[x.size * 2 - 2] = x[x.size - 1]
    return(xplot)
def multy(x):        
    xplot = np.zeros(x.size * 2 - 1)
    for i in range(0, x.size - 1):
        xplot[i * 2] = x[i]
        xplot[i * 2 + 1] = x[i]
    xplot[x.size * 2 - 2] = x[x.size - 1]
    return(xplot)
