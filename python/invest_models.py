from scipy.sparse import *
from scipy.sparse.linalg import *
from scipy import *
from scipy.sparse.linalg import spsolve
import numpy as np
import time

def water_quality(n, m, grid, E, Ux, Uy, K, s0, h):
    """2D Water quality model to track a pollutant in the ocean
    
    Keyword arguments:
    n,m -- the number of rows, columns in the 2D grid.  Used to determine 
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

    print 'initialize ...',
    t0 = time.clock()

    #used to abstract the 2D to 1D index calculation below
    def calc_index(i, j):
        if i >= 0 and i < n and j >= 0 and j < m:
            return i * m + j
        else:
            return - 1

    #set up variables to hold the sparse system of equations
    #upper bound  n*m*5 elements
    maxElements = n * m * 5
    col = np.empty(maxElements)
    row = np.empty(maxElements)
    data = np.empty(maxElements)
    b = np.zeros(n * m)
    currentIndex = 0
    print '(' + str(time.clock() - t0) + 's elapsed)'
    t0 = time.clock()
    #this map is used to quickly test if a neighboring element is a water cell

    #iterate over the non-zero elnp.array([])ts in grid to build the linear system
    print 'building system ...',
    t0 = time.clock()
    for i in range(n):
        for j in range(m):
            #diagonal element i,j always in bounds, calculate directly
            rowIndex = i * m + j

            #if not land, don't bother making an entry
            if not grid[rowIndex]:
                continue

            #if source, define value and quit
            if rowIndex in s0:
                row[currentIndex] = col[currentIndex] = rowIndex
                data[currentIndex] = 1
                currentIndex += 1
                b[rowIndex] = s0[rowIndex]
                continue

            #formulate elements as a single array
            termA = 2 * E[rowIndex]
            Uxtmp = Ux[rowIndex] * h
            Uytmp = Uy[rowIndex] * h

            elements = [
             (rowIndex, -4.0 * (termA + h * h * K[rowIndex])),
             (calc_index(i + 1, j), termA - Uxtmp),
             (calc_index(i - 1, j), termA + Uxtmp),
             (calc_index(i, j + 1), termA - Uytmp),
             (calc_index(i, j - 1), termA + Uytmp)]
            #process elements.  might be a source, might not...
            startIndex = currentIndex

            for colIndex, term in elements:
                if colIndex >= 0: #make sure we're in the grid
                    if grid[colIndex]:
                        row[currentIndex] = rowIndex
                        col[currentIndex] = colIndex
                        data[currentIndex] = term
                        currentIndex += 1
                    else:
                        #handle the land boundary case s_ij' = s_ij
                        data[startIndex] += term

    print '(' + str(time.clock() - t0) + 's elapsed)'
    #truncate the unused columns off of the sparse matrix arrays
    print 'building sparse matrix ...',
    t0 = time.clock()

    #truncate the elements we don't need
    col = np.delete(col, np.s_[currentIndex::])
    row = np.delete(row, np.s_[currentIndex::])
    data = np.delete(data, np.s_[currentIndex::])

    #create sparse matrix
    matrix = csc_matrix((data, (row, col)), shape=(n * m, n * m))

    print '(' + str(time.clock() - t0) + 's elapsed)'
    t0 = time.clock()

    print 'solving ...',
    result = spsolve(matrix, b)
    print '(' + str(time.clock() - t0) + 's elapsed)'
    return result
