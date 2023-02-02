from pickle import BINSTRING
import numpy
import matplotlib.pyplot as plt

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
    CDF = numpy.cumsum (binned) / numpy.cumsum(binned)[-1]
    return CDF

def xyoffset_to_r (x,y):
    ''' From cartesian offset x,y to radial distance from the centre \
        (assumed to be 0 by default) 
    '''
    return (numpy.sqrt(x**2 + y**2))

def radians_to_degree (radians):
    ''' Correctly scaled, but input of pi/2 means 1.5707963267948966, not 0.5.
    '''
    return radians * 180/numpy.pi

def degree_to_radians (degrees):
    ''' Already multiplied by pi
    '''
    return degrees * numpy.pi/180

