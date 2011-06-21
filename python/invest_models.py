from scipy.sparse import *
from scipy.sparse.linalg import *
from scipy import *
from scipy.sparse.linalg import spsolve
import numpy as np


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
    print 'building system'
    for i in range(n):
        for j in range(m):
            #diagonal element i,j
            rowIndex = calc_index(i, j)
            b.append(0) #initialize source vector

            #formulate elements as a single array 
            elements = [
             (calc_index(i, j), -(8.0 * E[rowIndex] + 2 * h * h * K[rowIndex])),
             (calc_index(i + 1, j), 2 * E[rowIndex] - Ux[rowIndex] * h),
             (calc_index(i - 1, j), 2 * E[rowIndex] + Ux[rowIndex] * h),
             (calc_index(i, j + 1), 2 * E[rowIndex] - Uy[rowIndex] * h),
             (calc_index(i, j - 1), 2 * E[rowIndex] + Uy[rowIndex] * h)]
            elementNum = 0
            #process elements.  might be a source, might not...
            for colIndex, term in elements:
                if colIndex >= 0:
                    if grid[colIndex]: #check if water
                        if colIndex not in s0: #check if source
                            row.append(rowIndex)
                            col.append(colIndex)
                            data.append(term)
                            elementNum += 1
                        else:
                            #test to see if current row is a source
                            if rowIndex == colIndex:
                                #we encounter the diagonal first because of 
                                #the way we defined elements.  Set that row
                                #to a 1 on the diagonal and put the source
                                #on the forcing function side
                                row.append(rowIndex)
                                col.append(colIndex)
                                data.append(1)
                                elementNum += 1
                                b[rowIndex] = s0[colIndex]
                                #no need to process the rest of the elements
                                break
                            else :
                                #we know what that variable is, so just push
                                #it's value to the b vector
                                b[rowIndex] += s0[colIndex] * (-term)
                    else:
                        #we're on a land point, define s as 0
                        if rowIndex == colIndex: #check if we're on a land point
                            row.append(rowIndex)
                            col.append(colIndex)
                            data.append(1)
                            elementNum += 1
                            #b is implicitly zero, no need to process the rest
                            #of the elements
                            break
                        else:
                            #handle the land boundary case s_ij' = s_ij
                            data[len(data) - elementNum] += term



    print 'building sparse matrix'
    #stamp into numpy formulation to be solved
    row = array(row)
    col = array(col)
    data = array(data)
    b = array(b)
    matrix = csc_matrix((data, (row, col)), shape=(n * m, n * m))

    print 'solving'
    return spsolve(matrix, b)

