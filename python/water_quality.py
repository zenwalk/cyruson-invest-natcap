from invest_models import *
import numpy as np
import matplotlib.pyplot as plt

def plotResult(land, result, n, m, dimx, dimy, h):
    print 'generate mesh'
    X, Y = np.meshgrid(np.arange(0, dimx, h), np.arange(0, dimy, h))



    print 'prepare result plot'
    result = np.array(result).reshape((n, m))
    print 'plot result'
    plt.pcolormesh(X, Y, result, cmap=plt.cm.gist_earth)
    plt.colorbar()

    if False in land:
        print 'generate land mask'
        Z = np.zeros(shape=(n, m), dtype=np.float)
        Zmask = np.ma.array(Z, mask=land)
        print 'plot land grid'
        plt.pcolormesh(X, Y, Zmask, cmap=plt.cm.BrBG)

def test2():
    #water quality test with all water

    #define land from file
    print 'load file'
    f = open('cell.inp', 'r')
    n = int(f.readline())
    m = int(f.readline())
    h = float(f.readline())
    gridText = f.readline()
    print n, m, h

    grid = []
    for c in gridText:
        if c == '1':
            grid.append(True)
        else:
            grid.append(False)

    print len(grid)

    dimx, dimy = m * h, n * h

    #define constants
    E = map(lambda x: 4.0, grid)
    Ux = map(lambda x:0.0707, grid)
    Uy = map(lambda x:0.0, grid)
    K = map(lambda x: 0.04, grid)

    #define a source right in the middle
    row = 291
    col = 342
    s0 = {row * m + col: 1}

    print row, col

    result2d = water_quality(n, m, grid, E, Ux, Uy, K, s0, h)
    result1d = water_quality_1d(n, m, E, Ux, Uy, K, s0, h)

    fig = plt.figure()
    f = fig.add_subplot(211, aspect='equal')
    plotResult(grid, result2d, n, m, dimx, dimy, h)
    f = fig.add_subplot(212, aspect='equal')
    plotResult(grid, result1d, n, m, dimx, dimy, h)


def test1():
    #water quality test with all water

    #cell size
    h = 0.5 #0.01km

    #square size
    dimx = 100
    dimy = 100

    n, m = int(dimy / h), int(dimx / h)

    #define land
    grid = [True] * n * m

    print "elements: ", n * m

    #define constants
    E = map(lambda x: 4, grid)
    Ux = map(lambda x:1.0, grid)
    Uy = map(lambda x:0.0, grid)
    K = map(lambda x: 0.1, grid)

    #define a source right in the middle
    row = int(n / 2)
    col = int(m / 2)
    s0 = {row * m + col: 1}

    result2d = water_quality(n, m, grid, E, Ux, Uy, K, s0, h)
    result1d = water_quality_1d(n, m, E, Ux, Uy, K, s0, h)

    fig = plt.figure()
    f = fig.add_subplot(211, aspect='equal')
    plotResult(grid, result2d, n, m, dimx, dimy, h)
    f = fig.add_subplot(212, aspect='equal')
    plotResult(grid, result1d, n, m, dimx, dimy, h)


if __name__ == "__main__":
    test1()
    test2()
    plt.show()
