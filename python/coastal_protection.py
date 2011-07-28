from invest_models import *
import numpy as np
import matplotlib.pyplot as plt

def plotResult(result, n, m, dimx, dimy, h):
    #print 'generate mesh'
    X, Y = np.meshgrid(np.arange(0, dimx, h), np.arange(0, dimy, h))

#    print 'prepare result plot'
    result = np.array(result).reshape((n, m))
#    print 'plot result'
    plt.pcolormesh(X, Y, result, cmap=plt.cm.gist_earth)
    plt.colorbar()

def test1():
    #coastal_protection test
    coastal_protection_2d(n,m,bathy,E,theta)

if __name__ == "__main__":
    test1()
