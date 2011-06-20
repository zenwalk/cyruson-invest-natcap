from scipy.sparse import *
from scipy import *
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt

def water_quality(n, m, grid, E, Ux, Uy, K, s0, h):
    """2D Water quality model to track a pollutant in the ocean
    
    Keyword arguments:
    n,m -- the number of rows,coluumns in the 2D grid.  Used to determine 
        indices into list parameters 'grid', 'E', 'Ux', 'Uy', and 'K' via
        sourceIndex (i,j) is position i*m+j in a list
    grid -- 1D list n*m elements long of booleans indicating land/water.  True
            is water, False is land.  
    E -- 1D list n*m elements long of dispersion coefficients
    Ux -- 1D list n*m elements long of x component velocity vectors
    Uy -- 1D list n*m elements long y component velocity vectors
    K -- 1D list n*m elements long of decay coefficients
    s0 -- map of sourceIndex to pollutant density
    h -- scalar describing grid cell size
    
    returns a 2D grid of pollutant densities in the same dimension as  'grid'
    
    """

    #used to abstract the 2D to 1D index calculation below
    def calc_index(i, j):
        if i >= 0 and i < n and j >= 0 and j < m:
            return i * m + j
        else:
            return - 1

    #set up variables to hold the sparse system of equations
    col = []
    row = []
    data = []
    b = []

    #this map is used to quickly test if a neighboring element is a water cell

    #iterate over the non-zero elements in grid to build the linear system
    for i in range(n):
        for j in range(m):
            #diagonal element i,j
            sourceIndex = calc_index(i, j)
            row.append(sourceIndex)
            col.append(sourceIndex)
            data.append(-(8.0 * E[sourceIndex] + 2 * h * h * K[sourceIndex]))
            b.append(0) #initialize source vector

            #formulate the nondiagonal elements as a single array 
            elements = [
             (calc_index(i + 1, j), 2 * E[sourceIndex] - Ux[sourceIndex] * h),
             (calc_index(i - 1, j), 2 * E[sourceIndex] + Ux[sourceIndex] * h),
             (calc_index(i, j + 1), 2 * E[sourceIndex] - Uy[sourceIndex] * h),
             (calc_index(i, j - 1), 2 * E[sourceIndex] + Uy[sourceIndex] * h)]

            #process nondiagonal elements.  might be a source, might not...
            for tmpIndex, term in elements:
                if tmpIndex >= 0:
                    if tmpIndex not in s0:
                        row.append(sourceIndex)
                        col.append(tmpIndex)
                        data.append(term)
                    else:
                        b[sourceIndex] += s0[tmpIndex] * (-term)

    #stamp into numpy formulation to be solved
    row = array(row)
    col = array(col)
    data = array(data)
    b = array(b)
    matrix = csr_matrix((data, (row, col)), shape=(n * m, n * m))
    return spsolve(matrix, b)

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
    E = map(lambda x: 1.0, grid)
    Ux = map(lambda x: 0.0, grid)
    Uy = map(lambda x: 0.0, grid)
    K = map(lambda x: 1.0, grid)

    #define a source right in the middle
    row = n / 2
    col = m / 2
    s0 = {row * m + col: 1}

    X, Y = np.meshgrid(np.arange(0, n) * h, np.arange(0, m) * h)

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

