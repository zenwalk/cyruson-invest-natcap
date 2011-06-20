from invest_models import water_quality
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    pass #put unit tests here

    #water quality test with all water
    #allow an n*m rectangular grid
    n, m = 100, 100

    #define land
    grid = [True] * n * m

    #cell size
    h = 0.01

    #define constants
    E = map(lambda x: 1, grid)
    Ux = map(lambda x:-11.0, grid)
    Uy = map(lambda x: 11.0, grid)
    K = map(lambda x: 1.5, grid)

    #define a source right in the middle
    row = n / 2
    col = m / 2
    s0 = {row * m + col: 100}

    X, Y = np.meshgrid(np.arange(0, m), np.arange(0, n))

    result = water_quality(n, m, grid, E, Ux, Uy, K, s0, h)

    #reroll into a 2D array
    Z = np.zeros(shape=(n, m), dtype=np.float)
    for i in range(n): Z[i, :] = result[i * m:i * m + m]
    print Z

    print X.ndim
    print Y.ndim
    print Z.ndim
    print X
    print Y
    print result
    print Z
    print s0
    fig = plt.figure()
    plt.pcolor(X, Y, Z)
    plt.colorbar()
    #fig.savefig('water.png', format='png')
    plt.show()

    print type(result)

