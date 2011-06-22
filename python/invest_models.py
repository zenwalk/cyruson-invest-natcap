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
    #upper bound  n*m*5 elements minus 4*the number of land cells
    maxElements = n * m * 5 - grid.count(True) * 4
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
    case1 = 0
    case2 = 0
    case3 = 0
    case4 = 0
    for i in range(n):
        for j in range(m):
            #diagonal element i,j
            rowIndex = calc_index(i, j)

            #check for boundary condition. if i,j on source or land
            if not grid[rowIndex] or rowIndex in s0:
                #not water, just set to 0 and quit
                row[currentIndex] = col[currentIndex] = rowIndex
                data[currentIndex] = 1
                currentIndex += 1
                if rowIndex in s0:
                    b[rowIndex] = s0[rowIndex]
                    case3 += 1
                case1 += 1
                continue

            #formulate elements as a single array 
            elements = [
             (calc_index(i, j), -(8.0 * E[rowIndex] + 2 * h * h * K[rowIndex])),
             (calc_index(i + 1, j), 2 * E[rowIndex] - Ux[rowIndex] * h),
             (calc_index(i - 1, j), 2 * E[rowIndex] + Ux[rowIndex] * h),
             (calc_index(i, j + 1), 2 * E[rowIndex] - Uy[rowIndex] * h),
             (calc_index(i, j - 1), 2 * E[rowIndex] + Uy[rowIndex] * h)]
            #process elements.  might be a source, might not...
            startIndex = currentIndex

            for colIndex, term in elements:
                if colIndex >= 0: #make sure we're in the grid
                    if not grid[colIndex]:
                        #handle the land boundary case s_ij' = s_ij
                        data[startIndex] += term
                        case2 += 1
                    else:
                        row[currentIndex] = rowIndex
                        col[currentIndex] = colIndex
                        data[currentIndex] = term
                        currentIndex += 1
                        case4 += 1

    print case1, case2, case3, case4
    print '(' + str(time.clock() - t0) + 's elapsed)'
    #truncate the unused columns off of the sparse matrix arrays
    print 'building sparse matrix ...',
    t0 = time.clock()
    col = np.delete(col, np.s_[currentIndex::])
    row = np.delete(row, np.s_[currentIndex::])
    data = np.delete(data, np.s_[currentIndex::])
    #stamp into numpy formulation to be solved
    matrix = csc_matrix((data, (row, col)), shape=(n * m, n * m))

    print '(' + str(time.clock() - t0) + 's elapsed)'
    print 'solving ...',
    t0 = time.clock()
    result = spsolve(matrix, b)
    print '(' + str(time.clock() - t0) + 's elapsed)'
    return result

