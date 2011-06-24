from invest_models import water_quality
import numpy as np
import matplotlib.pyplot as plt
import cProfile

def plotResult(land, result, n, m, dimx, dimy, h):
    print 'generate mesh'
    X, Y = np.meshgrid(np.arange(0, dimx, h), np.arange(0, dimy, h))

    fig = plt.figure()

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
    plt.axis('tight')
    fig.savefig('case1.png', format='png')
    plt.show()


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
    E = map(lambda x: 4.0, grid) #5 km/day 
    Ux = map(lambda x:0.0, grid) #8.64 km/day
    Uy = map(lambda x:0, grid)
    K = map(lambda x: 0.1, grid) #0.1%/day

    #define a source right in the middle
    row = 291
    col = 342
    s0 = {row * m + col: 1}

    print row, col

    result = water_quality(n, m, grid, E, Ux, Uy, K, s0, h)

    plotResult(grid, result, n, m, dimx, dimy, h)

def test1():
    #water quality test with all water

    #cell size
    h = 0.01 #0.01km

    #square size
    dimx = dimy = 3

    n, m = int(dimy / h), int(dimx / h)

    #define land
    grid = [True] * n * m

    print "elements: ", n*m

    #define constants
    E = map(lambda x: 0.1, grid) #5 km/day 
    Ux = map(lambda x:0.1, grid) #8.64 km/day
    Uy = map(lambda x:0.1, grid)
    K = map(lambda x: 0.2, grid) #0.1%/day

    #define a source right in the middle
    row = int(n / 2)
    col = int(m / 2)
    s0 = {row * m + col: 1}

    result = water_quality(n, m, grid, E, Ux, Uy, K, s0, h)
    plotResult(grid, result, n, m, dimx, dimy, h)

if __name__ == "__main__":
    test1()
    #cProfile.run('test1()', 'test1prof')
