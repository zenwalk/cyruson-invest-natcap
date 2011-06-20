from scipy.sparse import *
from scipy import *
from scipy.sparse.linalg import spsolve
import numpy as np


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
            b.append(0) #initialize source vector

            #formulate elements as a single array 
            elements = [
             (calc_index(i, j), -8.0 * E[sourceIndex] + 2 * h * h * K[sourceIndex]),
             (calc_index(i + 1, j), 2 * E[sourceIndex] - Ux[sourceIndex] * h),
             (calc_index(i - 1, j), 2 * E[sourceIndex] + Ux[sourceIndex] * h),
             (calc_index(i, j + 1), 2 * E[sourceIndex] - Uy[sourceIndex] * h),
             (calc_index(i, j - 1), 2 * E[sourceIndex] + Uy[sourceIndex] * h)]

            #process elements.  might be a source, might not...
            for tmpIndex, term in elements:
                if tmpIndex >= 0:
                    if tmpIndex not in s0:
                        row.append(sourceIndex)
                        col.append(tmpIndex)
                        data.append(term)
                    else:
                        b[sourceIndex] += s0[tmpIndex] * (-term)
                        if i == j:
                            row.append(sourceIndex)
                            col.append(sourceIndex)
                            data.append(1)

    #stamp into numpy formulation to be solved
    row = array(row)
    col = array(col)
    data = array(data)
    b = array(b)
    matrix = csr_matrix((data, (row, col)), shape=(n * m, n * m))
    return spsolve(matrix, b)

