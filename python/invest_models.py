import scipy
def water_quality(grid, E, U, K, s0, h):
    """2D Water quality model to track a pollutant in the ocean
    
    Keyword arguments:
    grid -- 2D grid of booleans indicating land/water.  True is water, False
            is land.
    E -- 2D grid of dispersion coefficients, must be same dimensions as 'grid'
    U -- 2D grid of velocity vectors, must be same dimensions as 'grid'
    K -- 2D grid of decay coefficients, must be same dimensions as 'grid'
    s0 -- 2D array of pollutant density, must be same dimensions as 'grid'
    h -- scalar describing grid cell size
    
    returns a 2D grid of pollutant densities in the same dimension as  'grid'
    
    """

    return True


if __name__ == "__main__":
    pass #put unit tests here
