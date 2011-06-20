from scipy.sparse import *
from scipy import *
from scipy.sparse.linalg import spsolve

def water_quality(n, m, grid, E, Ux, Uy, K, s0, h):
    """2D Water quality model to track a pollutant in the ocean
    
    Keyword arguments:
    n,m -- the number of rows,coluumns in the 2D grid.  Used to determine 
        indices into list parameters 'grid', 'E', 'Ux', 'Uy', and 'K' via
        index (i,j) is position i*m+j in a list
    grid -- 1D list n*m elements long of booleans indicating land/water.  True
            is water, False is land.  
    E -- 1D list n*m elements long of dispersion coefficients
    Ux -- 1D list n*m elements long of x component velocity vectors
    Uy -- 1D list n*m elements long y component velocity vectors
    K -- 1D list n*m elements long of decay coefficients
    s0 -- map of index to pollutant density
    h -- scalar describing grid cell size
    
    returns a 2D grid of pollutant densities in the same dimension as  'grid'
    
    """

    #set up variables to hold the sparse system of equations
    col = []
    row = []
    data = []
    b = []

    #this map is used to quickly test if a neighboring element is a water cell

    #iterate over the non-zero elements in grid to build the linear system
    for row in range(n):
        col, val = grid[row].nonzero()
        for i in range(len(col)):
            print col[i], val[i]

    matrix = csr_matrix((data, (row, col)), shape=(n, m))

    return True


if __name__ == "__main__":
    pass #put unit tests here

    #water quality test with all water
    #allow an n*m rectangular grid
    n, m = 10, 10

    #define land
    grid = [True] * n * m

    #cell size
    h = 0.01

    #define constants
    E = map(lambda x: 1.0, grid)
    Ux = map(lambda x: 1.0, grid)
    Uy = map(lambda x: 1.0, grid)
    K = map(lambda x: 1.0, grid)

    #define a source right in the middle
    row = array([n / 2])
    col = array([m / 2])
    so = {(row, col): 1}

    water_quality(grid, E, Ux, Uy, K, s0, h)

