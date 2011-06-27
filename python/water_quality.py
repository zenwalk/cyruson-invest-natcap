from invest_models import *
import numpy as np
import matplotlib.pyplot as plt

def plotResult(land, result, n, m, dimx, dimy, h):
    #print 'generate mesh'
    X, Y = np.meshgrid(np.arange(0, dimx, h), np.arange(0, dimy, h))



#    print 'prepare result plot'
    result = np.array(result).reshape((n, m))
#    print 'plot result'
    plt.pcolormesh(X, Y, result, cmap=plt.cm.gist_earth)
    plt.colorbar()

    if False in land:
        print 'generate land mask'
        Z = np.zeros(shape=(n, m), dtype=np.float)
        Zmask = np.ma.array(Z, mask=land)
        print 'plot land grid'
        plt.pcolormesh(X, Y, Zmask, cmap=plt.cm.BrBG)

def test4():
    #water quality test with all water

    #cell size
    h = 0.25 #0.01km

    #square size
    dimx = 30
    dimy = 30
    dt = 0.1
    tsteps = 50

    n, m = int(dimy / h), int(dimx / h)

    #define land
    grid = [True] * n * m

    print "elements: ", n * m

    #define constants
    E = map(lambda x: 4, grid)
    Ux = map(lambda x:4.0, grid)
    Uy = map(lambda x:4.0, grid)
    K = map(lambda x: 1.4, grid)

    #define a source right in the middle
    row = int(n / 2)
    col = int(m / 2)
    s0 = {row * m + col: 1}

    result2d_td = water_quality_time(n, m, tsteps, grid, E, Ux, Uy, K, s0, h, dt)

    fig = plt.figure()
    for step in range(tsteps):
        print 'ploting ' + str(step + 1) + ' of ' + str(tsteps)
        plt.clf()
        #f = fig.add_subplot("1"+str(tsteps)+str(step+1), aspect='equal')
        plotResult(grid, result2d_td[step], n, m, dimx, dimy, h)
        plt.suptitle('time: ' + str(dt * step))
        plt.savefig('out' + ('0' * int(2 - np.log(step + 1) / np.log(10))) + str(step) + '.png')

def test3():
    #time domain test

    #define land from file
    print 'load file'
    f = open('cell.inp', 'r')
    n = int(f.readline())
    m = int(f.readline())
    h = float(f.readline())
    gridText = f.readline()
    dt = 0.1
    tsteps = 5
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
    row = 342
    col = 291
    s0 = {row * m + col: 1}

    print row, col

    result2d_td = water_quality_time(n, m, tsteps, grid, E, Ux, Uy, K, s0, h, dt)
    fig = plt.figure()
    for step in range(tsteps):
        #f = fig.add_subplot("1"+str(tsteps)+str(step+1), aspect='equal')
        plt.clf()
        plt.suptitle('time: ' + str(dt * step))
        plotResult(grid, result2d_td[step], n, m, dimx, dimy, h)
        plt.savefig('out' + '0' * int(np.log(step + 1) / np.log(10)) + str(step + 1) + '.png')

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
    row = 342
    col = 291
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
    test4()
