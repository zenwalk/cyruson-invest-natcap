from scipy.sparse import *
from scipy import *
from scipy.sparse.linalg import spsolve

def water_quality(grid, E, Ux, Uy, K, s0, h):
    """2D Water quality model to track a pollutant in the ocean
    
    Keyword arguments:
    grid -- 2D grid of booleans indicating land/water.  True is water, False
            is land.
    E -- 2D grid of dispersion coefficients, must be same dimensions as 'grid'
    Ux -- 2D grid of x component velocity vectors, must be same dimensions as 'grid'
    Uy -- 2D grid of y component velocity vectors, must be same dimensions as 'grid'
    K -- 2D grid of decay coefficients, must be same dimensions as 'grid'
    s0 -- 2D array of pollutant density, must be same dimensions as 'grid'
    h -- scalar describing grid cell size
    
    returns a 2D grid of pollutant densities in the same dimension as  'grid'
    
    """

    return True


if __name__ == "__main__":
    pass #put unit tests here

    #water quality test with all water
    n = 10

    #define land
    row = array([x / n for x in range(n * n)])
    col = array([x % n for x in range(n * n)])
    data = array([True] * n * n)
    grid = csr_matrix((data, (row, col)), shape=(n, n))

    #cell size
    h = 0.01

    #define constants
    E = csr_matrix((map(lambda x: 1.0, data), (row, col)), shape=(n, n))
    Ux = csr_matrix((map(lambda x: 1.0, data), (row, col)), shape=(n, n))
    Uy = csr_matrix((map(lambda x: 1.0, data), (row, col)), shape=(n, n))
    K = csr_matrix((map(lambda x: 1.0, data), (row, col)), shape=(n, n))

    #define source
    row = array([n / 2])
    col = array([n / 2])
    data = array([1])
    s0 = csr_matrix((data, (row, col)), shape=(n, n))



    water_quality(grid, E, Ux, Uy, K, s0, h)

