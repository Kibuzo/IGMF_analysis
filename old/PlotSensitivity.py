import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from gammapy.spectrum.models import Absorption
import astropy.units as u
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy import stats
from scipy import integrate
from gammapy.spectrum import FluxPoints
from astropy.table import Table
from gammapy.spectrum.models import PowerLaw
from gammapy.utils.fitting import Fit
from iminuit import Minuit
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import math
Ezero=300.
Nzero=3.83e-14
gamma=1.5
E=np.logspace(2,5,num=100)

redshift=0.2
dominguez = Absorption.read_builtin('dominguez').table_model(redshift)
franceschini = Absorption.read_builtin('franceschini').table_model(redshift)
finke = Absorption.read_builtin('finke').table_model(redshift)

def assorbi (energy,flux,model):
    a=[]
    for j in range (0,len(energy)):
        a.append(flux[j]*model.evaluate(energy[j]*u.GeV,1))
    return(a)


savedir='/home/kibuzo/TesiMagistrale/Python/Data/'

def ErgToTeV(erg):
    return 0.624151*erg

CTA=np.loadtxt(savedir+'CTA_50h.txt', unpack=True)
HESS=np.loadtxt(savedir+'HESS_25h.txt', unpack=True)

HESS50=HESS[1]*0.5*(HESS[0]**2)*1e3
ctacentre=(CTA[1]+CTA[0])/2
CTA50=ErgToTeV(CTA[2])*1e3#/ctacentre**2)

plt.plot(HESS[0]*1e3,2*HESS50, label='HESS 25h')
plt.plot(ctacentre*1e3,2*CTA50,label='CTA south 25h')
plt.plot(E,assorbi(E,(E**2)*Nzero*(E/300)**-gamma,franceschini), label='Hypotesized SED (No cutoff)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy [GeV]')
plt.ylabel('E$^2$ flux sensitiity (GeVcm$^{-2}$s$^{-1}$)')
plt.legend()
plt.show()