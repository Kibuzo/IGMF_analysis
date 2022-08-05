from astropy.io import fits
from astropy import units
import numpy

'''The old process worked like this:
1) Open the fit scirpt (FittaCherenkiSomething) and execute it, mostly.
2) Within this file you had the x,y,xerr,yerr vectors defined somehow by 
   hand
3) You decided which one to use, provided that it didn't change much, by 
   memory
4) Something in this file united all the vectors (Fermi, VERITAS, HESS) 
   and created xcombnine, ycombined, etc.

We should really aim to rely on reading fits files for the Fermipy analysis 
and use hand-made vectors only for the salvaged Cherenkov data so that a 
new analysis can update the results on the fly. This might include also 
extending the pipeline to newer data since >2 years have passed.
This file exists for that specific reason.
'''

file_name = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/Backup HESS/Fermipy/FromScratch/SED.fits'

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
      err_low = hdul[1].data['dnde_err']*scaling*2
      err_high = err_low
      ts = hdul[1].data['ts']
      ulims = hdul[1].data['dnde_ul']*scaling*2
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

