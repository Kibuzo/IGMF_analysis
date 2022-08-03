from gammapy.spectrum.models import Absorption
import astropy.units as u

redshift=0.2

dominguez = Absorption.read_builtin('dominguez').table_model(redshift)
franceschini = Absorption.read_builtin('franceschini').table_model(redshift)
finke = Absorption.read_builtin('finke').table_model(redshift)

def deassorbi(energy, flux, model):
    a=[]
    for j in range (0,len(flux)):
        a.append(flux[j]/model.evaluate(energy[j]*u.GeV,1))
    return(a)

def assorbi (energy,flux,model):
    a=[]
    for j in range (0,len(energy)):
        a.append(flux[j]*model.evaluate(energy[j]*u.GeV,1))
    return(a)
    
def degtorad(deg):
    return deg*np.pi/180.
    
def radtodeg(rad):
    return rad*180/np.pi