import numpy as np
import matplotlib.pyplot as plt

#Parametri della logparabola del catalogo
alpha=1.3095159530639648
beta=0.351797878742218
norm=0.5439949672095181*1e-14
Eb=15794.2490234375
def propagaspettro(spettro,energia,dspettro,denergia):
    return (np.sqrt((energia**4)*(dspettro**2)+4*(spettro**2)*(energia**2)*(denergia**2)+4*(dspettro*denergia*spettro*(energia)**3)))

#$FGL psc_v18 Extended sources P8R3IRF
nphot=[2.07e-11,5.69e-12,2.11e-12,7.3e-13,2.688e-13,2.345e-14,9.9e-15]
Fermix=np.array([6.94,13.41,25.9,50,96.54,186.38,359.84]) #in GeV
Fermixmin=np.array([5,9.65,18.6,35.98,69.47,134.13,258.97])
Fermixmax=np.array([9.65,18.6,35.98,69.47,134.13,258.97,500])
dNphot=np.array([4.5e-12,1.52e-12,5.9e-13,2.31e-13,9.95e-14,2.26e-14/2,1.01e-14/2])
Fermidx=(Fermixmax-Fermixmin).tolist()
fermiuplims=np.array([0,0,0,0,0,1,1])
Nphot=np.array(nphot)
deltafermi=propagaspettro(Nphot,Fermix,dNphot/2,np.array(Fermidx)/2) #Corretto



e=np.logspace(0.6,2.7) #in GeV
DnDE=[]
eMeV=e*1e3
for E in eMeV:
    x=E/Eb
    DnDE.append(norm*x**(-(alpha+beta*np.log10(x))))
DnDE=np.array(DnDE)
plt.errorbar(Fermix,nphot*Fermix**2,yerr=deltafermi,marker='s', markersize=2,uplims=fermiuplims,linestyle='',label='Fermi')
    
plt.loglog(e,DnDE*1e3*e**2)
plt.show()
    