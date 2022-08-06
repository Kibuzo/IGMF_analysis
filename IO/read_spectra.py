from astropy.io import fits
from astropy import units
import matplotlib.pyplot as plt
from gammapy.modeling.models import EBLAbsorptionNormSpectralModel as absorb
import numpy

from Constants.J1943 import REDSHIFT

'''The old process worked like this:
1) Open the fit scirpt (FittaCherenkiSomething) and execute it, mostly.
2) Within this file you had the x,y,xerr,yerr vectors defined somehow by 
   hand
3) You decided which one to use, provided that it didn't change much, by 
   memory
4) Something in this file united all the vectors (Fermi, VERITAS, HESS) 
   and created xcombnine, ycombined, etc.

We are in the middle of rethinking the process entirely. 

- Now The spectrum class is capable of reading the data and is ready to 
host a plot method, a fit method, and accomodate for parameters so that 
the user can read it.
- The VERITAS/HESS columns will be soon written to a fits file and the 
reading will be unified so that there will be only one spectrum class and
will accept the file_path as argument
- More methods to come.
'''

file_name = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/Backup HESS/Fermipy/FromScratch/SED.fits'
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
   
   def _has_ts(self):
      '''Search for upper limits column
      '''
      return (self.ts is not None)
   
   def plot(self, **kwargs):
      '''Utility to plot spectra in the same format
      '''
      if self._has_ts():
         mask=self.ts > 9 
         umask = ~mask
         ulims = self.ulims[umask]
         energy_ulims = self.energy[umask]
         plt.errorbar(energy_ulims, ulims, yerr = ulims/2, fmt ='',
                      linestyle = 'None', uplims = True)
      else:
         mask = numpy.ones(len(self.flux), dtype = bool)
      #mask = numpy.ones(len(self.flux), dtype = bool)
      #umask = ~mask
      energy = self.energy[mask]
      flux = self.flux[mask]
      err_hi = self.err_high[mask]/2
      err_lo = self.err_low[mask]/2
      errorbars = (err_hi, err_lo)
      plt.errorbar(energy, flux, yerr=errorbars, markersize=4,\
                 elinewidth = 1, linestyle ='None', **kwargs)
      plt.grid(linestyle = '--')
      plt.xscale('log')
      plt.yscale('log')
      plt.legend()

   def deabsorb(self):
      self.flux = self.flux*numpy.exp(dominguez.spectral_index(self.energy*units.GeV))


   def absorb (self):
      self.flux = self.flux*numpy.exp(-dominguez.spectral_index(self.energy*units.GeV))
    

   


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
      spectrum.__init__(self, flux, err_high, err_low,\
         ulims, energy, ts)
      self.energy_max = hdul[1].data['e_max']/scaling
      self.energy_min = hdul[1].data['e_min']/scaling


class veritas(spectrum):
   ''' This is sadly based on vectors instead of files.
   The flux is absorbed.

   Sadness apart, it works. See if you can save everything in a file.
   '''
   def __init__(self):
      '''Constructor
      '''
      scaling = 1/units.MeV.to(units.GeV)#From 1/MeV to 1/GeV
      energy = numpy.array([0.21, 0.299, 0.421, 0.592, 0.836,\
          1.18, 1.673])*scaling
      flux = numpy.array([7.22e-11, 4.12e-11, 1.76e-11, 1.26e-11,\
          9.87e-12, 4e-12, 3.49e-12])/scaling
      err_high = numpy.array([7.22e-11, 4.12e-11, 1.76e-11,\
          1.26e-11, 9.87e-12, 4e-12, 3.49e-12])/scaling
      err_low = numpy.array([1.2e-11, 6.13e-12, 3.04e-12, 2.15e-12,\
          1.49e-12, 1.28e-12, 1.15e-12])/scaling
      ulims = []
      spectrum.__init__(self, flux, err_high, err_low,\
         ulims, energy)
      

class hess(spectrum):
   ''' This is sadly based on vectors instead of files.
   The flux is absorbed.

   Sadness apart, it works. See if you can save everything in a file.
   '''
   def __init__(self):
      '''Constructor
      '''
      scaling = 1/units.MeV.to(units.GeV)#From 1/MeV to 1/GeV
      energy = numpy.array([0.538, 0.702, 0.934, 1.248, 1.67,\
         2.26, 3.018, 4.058, 5.4, 7.3, 9.8, 13])*scaling
      flux = numpy.array([1.30311e-11, 7.44939e-12, 6.42965e-12,\
          2.95444e-12, 1.78543e-12, 1.35665e-12, 4.30136e-13, 1.32585e-12,\
          2.9e-14, 2.4e-15,2.4e-14,4.8e-15])/scaling
      err_high = numpy.array([2.64e-12, 2.52e-12, 1.46e-12,\
          1.09e-12, 8.09e-13, 7.42e-13, 5.06e-13, 6.96e-13,0,0,0,0])\
          /scaling
      err_low = numpy.array([5.87e-12, 2.35e-12, 1.69e-12, 9.54e-13,\
          7.11e-13, 6.01e-13, 3.24e-13, 6.13e-13, 0, 0, 0, 0])/scaling
      ulims = []
      spectrum.__init__(self, flux, err_high, err_low,\
         ulims, energy)

FermiSpec = fermi(file_name)
#FermiSpec.plot(label = 'Fermi', marker = 'v')
FermiSpec.deabsorb()
FermiSpec.plot(label = 'Fermi', marker = 'v')
VeritasSpec = veritas()
VeritasSpec.plot(label = 'Veritas', marker = '.')
HessSpec = hess()
HessSpec.plot(label = 'Hess', marker = 'x')
plt.show()