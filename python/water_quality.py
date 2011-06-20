from invest_models import water_quality
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    pass #put unit tests here

    #water quality test with all water
    #allow an n*m rectangular grid
    n, m = 3, 3

    #define land
    grid = [True] * n * m

    #cell size
    h = 0.01

    #define constants
    E = map(lambda x: 0.1, grid)
    Ux = map(lambda x: 0.0, grid)
    Uy = map(lambda x: 0.0, grid)
    K = map(lambda x: 0.1, grid)

    #define a source right in the middle
    row = n / 2
    col = m / 2
    s0 = {row * m + col: 1}

    X, Y = np.meshgrid(np.arange(0, n), np.arange(0, m))

    result = water_quality(n, m, grid, E, Ux, Uy, K, s0, h)

    #fast method to re-roll a numpy array into 2D
    def rolling(a, window):
        shape = (a.size - window + 1, window)
        strides = (a.itemsize, a.itemsize)
        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

    Z = rolling(result, m)
    print X.ndim
    print Y.ndim
    print Z.ndim
    plt.figure()
    plt.pcolor(X, Y, Z)
    plt.colorbar()
    plt.show()

    print type(result)

