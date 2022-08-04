from pickle import BINSTRING
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

def expcutoff_cdf(E, N0, index, divs, Ecut):
    '''Returns the cumulative of an exponential cutoff function with
    divs bins. 

    ***Warning: this has changed substantially wrt to the old CDF funcion
    defined in NewMontecarlo.py, and needs to be called differently.
    '''
    expcutoff = cutoff(E, N0, index, Ecut)
    binned = numpy.histogram (expcutoff, divs)
    CDF = numpy.cumsum (binned) / numpy.cumsum(binnes)[-1]
    return CDF

