import numpy
import matplolib.pyplot as plt

def cutoff(E, N0, index, Ecut, E0=300.):
    '''Power law with energy cutoff
    '''
    return (N0 * (E/E0)**-index)*numpy.exp(-E/Ecut)

def grid_to_E (grid_point, extent, E_step):
    ''' Convert a point on an energy grid to physical units
    '''
    return (extent[0] + grid_point * E_step)

def E_to_grid (E, extent, E_step):
    ''' Convert physical unit to grid points
    '''
    return (E - extent[0]) / E_step
