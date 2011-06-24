from scipy.sparse import *
from scipy import *
from scipy.sparse.linalg import spsolve
import numpy as np
import time
import matplotlib
import scipy.linalg
import math

def water_quality_1d(n, m, E, Ux, Uy, K, s0, h):
    """1D analytical solution to water quality model
    
    Keyword arguments:
    n,m -- the number of rows, columns in the 2D grid.  Used to determine 
        indices into list watermeters 'E', 'Ux', 'Uy', and 'K' via
        sourceIndex (i,j) is position i*m+j in a list
    E -- 1D list n*m elements long of dispersion coefficients
    Ux -- 1D list n*m elements long of x component velocity vectors
    Uy -- 1D list n*m elements long y component velocity vectors
    K -- 1D list n*m elements long of decay coefficients
    s0 -- map of sourceIndex to pollutant density
    h -- scalar describing grid cell size
    
    returns a rough approximation to 2D via forcing 1d rays all over the grid
    
    """

    t0 = time.clock()
    print 'doing 1d solution ...',

    x0 = np.zeros(n * m)
    for i in range(n):
        for j in range(m):
            index = i * m + j
            if Ux[index] == 0:
                continue
            alpha = Ux[index] / (2 * E[index]) * (1 - math.sqrt(1 + 4 * K[index] * E[index] / (Ux[index] * Ux[index])))
            for p in s0:
                i0, j0 = p % m, p / m
                val = s0[p]
                d = math.sqrt(math.pow(h * (i - i0), 2) + math.pow(h * (j - j0), 2))
                x0[index] += val * math.exp(d * alpha)
    return x0

    print '(' + str(time.clock() - t0) + 's elapsed)'

def water_quality(n, m, inWater, E, Ux, Uy, K, s0, h, directSolve=False):
    """2D Water quality model to track a pollutant in the ocean
    
    Keyword arguments:
    n,m -- the number of rows, columns in the 2D grid.  Used to determine 
        indices into list parameters 'water', 'E', 'Ux', 'Uy', and 'K' i*m+j in
        a list
    water -- 1D list n*m elements long of booleans indicating land/water.  True
            is water, False is land.  
    E -- 1D list n*m elements long of dispersion coefficients
    Ux -- 1D list n*m elements long of x component velocity vectors
    Uy -- 1D list n*m elements long y component velocity vectors
    K -- 1D list n*m elements long of decay coefficients
    s0 -- map of sourceIndex to pollutant density
    h -- scalar describing grid cell size
    directSolve -- if True uses a direct solver that may be faster, but use
        more memory.  May crash in cases where memory is fragmented or low
        Default False.
    
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
    b = np.zeros(n * m)
    #holds the columns for diagonal sparse matrix creation later
    A = np.zeros((5, n * m))

    print '(' + str(time.clock() - t0) + 's elapsed)'
    t0 = time.clock()

    #iterate over the non-zero elments in grid to build the linear system
    print 'building system A...',
    t0 = time.clock()
    for i in range(n):
        for j in range(m):
            #diagonal element i,j always in bounds, calculate directly
            rowIndex = calc_index(i, j)

            #if land then s = 0 and quit
            if not inWater[rowIndex]:
                A[2, rowIndex] = 1
                continue

            #formulate elements as a single array
            termA = 2 * E[rowIndex]
            Uxtmp = Ux[rowIndex] * h
            Uytmp = Uy[rowIndex] * h

            elements = [
             (2, 0, rowIndex, -4.0 * (termA + h * h * K[rowIndex])),
             (4, m, calc_index(i + 1, j), termA - Uytmp),
             (0, -m, calc_index(i - 1, j), termA + Uytmp),
             (3, 1, calc_index(i, j + 1), termA - Uxtmp),
             (1, -1, calc_index(i, j - 1), termA + Uxtmp)]

            for k, offset, colIndex, term in elements:
                if colIndex >= 0: #make sure we're in the grid
                    if inWater[colIndex]: #if water
                        A[k, rowIndex + offset] += term
                    else:
                        #handle the land boundary case s_ij' = s_ij
                        A[2, rowIndex] += term

    #define sources by erasing the rows in the matrix that have already been set
    for rowIndex in s0:
        #the magic numbers are the diagonals and their offsets due to gridsize
        for i, offset in [(4, m), (0, -m), (3, 1), (1, -1)]:
            #zero out that row
            A[i, rowIndex + offset] = 0
        A[2, rowIndex] = 1
        b[rowIndex] = s0[rowIndex]
    print '(' + str(time.clock() - t0) + 's elapsed)'

    print 'building sparse matrix ...',
    t0 = time.clock()
    matrix = spdiags(A, [-m, -1, 0, 1, m], n * m, n * m, "csc")
    print '(' + str(time.clock() - t0) + 's elapsed)'

    if directSolve:
        t0 = time.clock()
        print 'direct solving ...',
        result = spsolve(matrix, b)
    else:
        print 'generating preconditioner via sparse ilu ',
        #normally factor will use m*(n*m) extra space, we restrict to 
        #\sqrt{m}*(n*m) extra space
        P = scipy.sparse.linalg.spilu(matrix, fill_factor=int(math.sqrt(m)))
        print '(' + str(time.clock() - t0) + 's elapsed)'
        t0 = time.clock()
        print 'gmres iteration starting ',
        #create linear operator for precondioner
        M_x = lambda x: P.solve(x)
        M = scipy.sparse.linalg.LinearOperator((n * m, n * m), M_x)
        result = scipy.sparse.linalg.lgmres(matrix, b, tol=1e-5, M=M)[0]
    print '(' + str(time.clock() - t0) + 's elapsed)'
    return result
