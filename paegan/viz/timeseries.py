import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import pylab

def progressive_vectors(time, U, V, save=None, axes=None, **kwargs):
    # Create figure and axes from matplotlib
    if axes == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        ax = axes
    #print time.shape, U.shape, V.shape, np.zeros_like(U).shape
    handle = ax.quiver(time, np.zeros_like(U), U, V,
                       pivot='tail', **kwargs)
    
    if save == None:
        pass
    else:
        fig.savefig(save, dpi=200)
    return handle
        
    
