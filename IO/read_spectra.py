from astropy.io import fits

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

#file_name = '/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/Backup HESS/Fermipy/FromScratch/SED.fits'

class fermipy:
   '''Read/output fermipy file relevant data, converted in GeVs
   '''
   def __init__(self, file_path):
      '''Constructor
      '''
      self.hdul = fits.open(file_path)
      self.hdul.info()
      self.scaling = 1e3

   
   def dnde(self):
      '''Number of photons in the energy bin
      '''
      return self.hdul[1].data['dnde']*self.scaling
   
   def dnde_err(self):
      '''number of photons in the energy bin
      '''
      return self.hdul[1].data['dnde_err']*self.scaling
   
   def ts(self):
      '''self-explanatory
      '''
      return self.hdul[1].data['ts']
   
   def dnde_ulims(self):
      '''self-explanatory
      '''
      return self.hdul[1].data['dnde_ul']*self.scaling
   
   def energy(self):
      '''self-explanatory
      '''
      return self.hdul[1].data['e_ref']/self.scaling

   def energy_max(self):
      '''self-explanatory
      '''
      return self.hdul[1].data['e_max']/self.scaling

   def energy_min(self):
      '''self-explanatory
      '''
      return self.hdul[1].data['e_min']/self.scaling