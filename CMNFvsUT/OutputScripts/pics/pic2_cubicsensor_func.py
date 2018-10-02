import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

X = np.arange(-50.0, 50.0, 0.1)
Y = X / (1 + X*X)

fig = plt.figure(figsize=(5, 5), dpi=200)
_ = plt.plot(X, Y, linewidth=1.0, color='black')
plt.tight_layout()
plt.show()
        
