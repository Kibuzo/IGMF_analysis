from astropy.io import fits
from astropy import units
from astropy.table import Table
import matplotlib.pyplot as plt
from gammapy.modeling.models import EBLAbsorptionNormSpectralModel as absorb
import numpy
# pylint: disable=invalid-name
# pylint: disable=no-member

from Constants.J1943 import REDSHIFT

#The old process worked like this:
#1) Open the fit scirpt (FittaCherenkiSomething) and execute it, mostly.
#2) Within this file you had the x,y,xerr,yerr vectors defined somehow by
#   hand
#3) You decided which one to use, provided that it didn't change much, by
#   memory
#4) Something in this file united all the vectors (Fermi, VERITAS, HESS)
#   and created xcombnine, ycombined, etc.
#
#We are in the middle of rethinking the process entirely.
#
#- Now The spectrum class is capable of reading the data and plotting, and
#is ready to host a fit method, and accomodate for parameters so that
#the user can read them.
#- Veritas/HESS columns can now be written to a file and the generic class
#from_file sort of works (still needs the correct scaling for the energies
#for compatibility with fermi files)
#- More methods to come.

file_name = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/\
            Backup HESS/Fermipy/FromScratch/SED.fits'
dominguez = absorb.read_builtin('dominguez', redshift=REDSHIFT)

class spectrum:
    ''' Spectrum superclass containing a common interface and basic
    functions
    '''
    def __init__(self, flux, err_high, err_low, ulims, energy, ts = None):
        '''Constructor
        '''
        self.flux = flux
        self.err_high = err_high
        self.err_low = err_low
        self.ts = ts
        self.ulims = ulims
        self.energy = energy

    def _has_ulims(self):
        '''Search for upper limits column
        '''
        return numpy.any (self.ts < 9)

    def plot(self, plot_ulims = False, **kwargs):
        '''Utility to plot spectra in the same format
        '''
        if self._has_ulims():
            mask=self.ts > 9
            umask = ~mask
            ulims = self.ulims[umask]
            energy_ulims = self.energy[umask]
            if plot_ulims:
                plt.errorbar(energy_ulims, ulims, yerr = ulims/2, fmt ='',
                      linestyle = 'None', uplims = True)
        else:
            mask = numpy.ones(len(self.flux), dtype = bool)
        energy = self.energy[mask]
        flux = self.flux[mask]
        err_hi = self.err_high[mask]/2
        err_lo = self.err_low[mask]/2
        errorbars = (err_hi, err_lo)
        plt.errorbar(energy, flux, yerr=errorbars, markersize=4,
                     elinewidth = 1, linestyle ='None', **kwargs)
        plt.grid(linestyle = '--')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('[Energy [GeV]]')
        plt.ylabel('Flux [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
        plt.legend()

    def deabsorb(self):
        '''De-absorb the spectrum with a dominguez EBL model.
        '''
        self.flux = self.flux*numpy.exp(dominguez.spectral_index(self.energy*units.GeV))

    def absorb (self):
        '''Absorb the spectrum with a dominguez EBL model.
        '''
        self.flux = self.flux*numpy.exp(-dominguez.spectral_index(self.energy*units.GeV))

    def write (self, name, overwrite = 'False'):
        ''' Write a unified fits file from each possible class. This is
        designed to work in the same fashion for each possible spectrum,
        so that the units, column names etc are unified and a single class
        can be used to load them all by justy defining another path.
        The ts column will be a placeholder for most extracted data, and is
        designed merely to work with the implemented filters: 100 for data
        points and 0 for upper limits.

        Please, leave the old functions for reference and debugging issues.
        '''
        t = Table ([self.energy, self.flux, self.err_high, self.err_low,\
            self.ulims, self.ts], names = ('energy','flux','err_high','err_low',\
            'ulims','ts'))
        t.write(name, format = 'fits', overwrite = overwrite)

class from_fits(spectrum):
    '''Generic utility to load SED from self-processed fits file
    '''
    def __init__(self, file_path):
        '''Constructor
        '''
        hdul = fits.open(file_path)
        hdul.info()
        flux = hdul[1].data['flux']
        err_low = hdul[1].data['err_low']
        err_high = hdul[1].data['err_high']
        ts = hdul[1].data['ts']
        ulims = hdul[1].data['ulims']
        energy = hdul[1].data['energy']
        spectrum.__init__(self, flux, err_high, err_low,
                          ulims, energy, ts)

class fermi(spectrum):
    '''Read/output fermipy file relevant data, converted in GeVs
    '''
    def __init__(self, file_path):
        '''Constructor
        '''
        hdul = fits.open(file_path)
        hdul.info()
        scaling = 1/units.MeV.to(units.GeV)#From 1/MeV to 1/GeV
        flux = hdul[1].data['dnde']*scaling
        err_low = hdul[1].data['dnde_err']*scaling
        err_high = err_low
        ts = hdul[1].data['ts']
        ulims = hdul[1].data['dnde_ul']*scaling
        energy = hdul[1].data['e_ref']/scaling
        spectrum.__init__(self, flux, err_high, err_low,
                          ulims, energy, ts)
        self.energy_max = hdul[1].data['e_max']/scaling
        self.energy_min = hdul[1].data['e_min']/scaling

class veritas(spectrum):
    ''' This is sadly based on vectors instead of files.
    The flux is intrinsic, and is obtained from
    https://iopscience.iop.org/article/10.3847/1538-4357/aacbd0/pdf
    and then de-absorbed (see the hess class for method).

    Sadness apart, it works.
    '''
    def __init__(self):
        '''Constructor
        '''
        scaling = 1/units.MeV.to(units.GeV)#From 1/MeV to 1/GeV
        energy = numpy.array([0.21, 0.299, 0.421, 0.592, 0.836,
                             1.18, 1.673])*scaling
        flux = numpy.array([7.22e-11, 4.12e-11, 1.76e-11, 1.26e-11,
                             9.87e-12, 4e-12, 3.49e-12])/scaling
        err_high = numpy.array([7.22e-11, 4.12e-11, 1.76e-11,
                               1.26e-11, 9.87e-12, 4e-12, 3.49e-12])/scaling
        err_low = numpy.array([1.2e-11, 6.13e-12, 3.04e-12, 2.15e-12,
                              1.49e-12, 1.28e-12, 1.15e-12])/scaling
        ulims = numpy.zeros(len(err_low))
        ts = 100. + numpy.zeros(len(err_low))
        spectrum.__init__(self, flux, err_high, err_low,
                          ulims, energy, ts)

class hess(spectrum):
    ''' This is sadly based on vectors instead of files.
    The flux is intrinsic, obtained from
    https://iopscience.iop.org/article/10.3847/1538-4357/aacbd0/pdf
    Achtung (this took me a while to debug): the points that you see here
    are not taken directly from the paper, rather they have been de-absorbed
    with an EBL model from dominguez (such as the one we use in this very
    script). We shall provide a unit test for this. The original vector of
    fluxes in units of cm^-2·s^-2·Tev^-1 is provided here

    flux = [3.94e-12, 1.44e-12, 7.86e-13, 2.5e-13, 
                            1e-13, 5.8e-14, 1.3e-14, 2.1e-14, 2.9e-15, 
                            2.3e-16, 2.4e-15, 4.8e-16]

    and one can see that, apart from a scaling factor, it trends exactly
    as the flux column from the HESS fits file, as opposed to de-absorbing
    or absorbing it 

    Despite being written in a despicable way, it works.
    '''
    def __init__(self):
        '''Constructor
        '''
        scaling = 1/units.MeV.to(units.GeV)#From 1/MeV to 1/GeV
        energy = numpy.array([0.538, 0.702, 0.934, 1.248, 1.67,
                              2.26, 3.018, 4.058, 5.4, 7.3, 9.8, 13])*scaling
        flux = numpy.array([1.30311e-11, 7.44939e-12, 6.42965e-12,
          2.95444e-12, 1.78543e-12, 1.35665e-12, 4.30136e-13, 1.32585e-12,
          2.9e-14, 2.4e-15,2.4e-14,4.8e-15])/scaling
        err_high = numpy.array([2.64e-12, 2.52e-12, 1.46e-12,
          1.09e-12, 8.09e-13, 7.42e-13, 5.06e-13, 6.96e-13,0,0,0,0])/scaling
        err_low = numpy.array([5.87e-12, 2.35e-12, 1.69e-12, 9.54e-13,
          7.11e-13, 6.01e-13, 3.24e-13, 6.13e-13, 0, 0, 0, 0])/scaling
        ulims = flux
        ts = 100. + numpy.zeros(len(err_low))
        ts[-4:] = 0
        spectrum.__init__(self, flux, err_high, err_low,
                          ulims, energy, ts)

VeritasSpec = from_fits('data/Veritas.fits')
HessSpec = from_fits('data/Hess.fits')
FermiSpec = from_fits('data/Fermi.fits')

FermiSpec.deabsorb()
Hspec = hess()

plt.figure('Multifrequency spectrum')
VeritasSpec.plot(label = 'Veritas', marker = '.')
HessSpec.plot(label = 'Hess', marker = 'x')
FermiSpec.plot(label = 'Fermi', marker = 'v')

plt.title('Intrinsic multifrequency spectrum of HESS J1943+213')
plt.legend()

plt.show()
