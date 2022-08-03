from pylab import *
import os
import numpy as np
import time
from scipy.optimize import curve_fit
savedir='./BPL/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
savedirtmp='./BPL/TEMP/'
if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)

redshift=0.21
dominguez = Absorption.read_builtin('dominguez').table_model(redshift)
z1es=0.14
dominguez1es = Absorption.read_builtin('dominguez').table_model(z1es)

def deassorbi(energy, flux, model):
    a=[]
    for j in range (0,len(flux)):
        #a.append(flux[j]*np.exp(model.spectral_index(energy[j]*u.GeV,epsilon=1e-05)))
        a.append(flux[j]/model.evaluate(energy[j]*u.GeV,1))
    return(a)

def assorbi (energy,flux,model):
    a=[]
    for j in range (0,len(energy)):
        a.append(flux[j]*model.evaluate(energy[j]*u.GeV,1))
    return(a)

puntidatogliere=6


#ESx=np.array([291,427,628,922,1360,1990,2920,4290,7640])
#ESy=np.array([1.2e-10,6.1e-11,2.2e-11,4.3e-12,2.2e-12,1e-12,3.5e-13,2.4e-13,2.8e-14])*1e-4 #Per trasformare in cm^-2
#ESdy=np.array([2.8e-11,1e-11,4.4e-12,1.7e-12,8.1e-13,4.2e-13,2.1e-13,1.2e-13,2.5e-14])*1e-4 #idem
ESx=([599,875,1284,1886,2762,4065,6464,11499])
ESy=np.array([3.315e-12, 9.193e-13, 2.343e-13, 1.622e-13, 7.372e-14, 2.444e-14, 3.841e-15, 1.415e-15])*1e-4
ulimsEsy=np.array([1.196e-12,3.350e-13, 2.088e-13, 9.700e-14, 3.532e-14, 6.994e-15, 2.470e-15, 4.54e-12])*1e-4
lowlimsEsy=np.array([6.430e-13, 1.300e-13, 1.158e-13, 5.049e-14, 1.214e-14, 6.037e-16, 4.06e-16, 1.87e-12])*1e-4
ESdy=(ulimsEsy-lowlimsEsy)
ES=[]
ES.append(ESx)
ES.append(ESy)
ES.append(ESdy)
 #np.savetxt(savedir+'1es0229_data.txt',ES)
Mio=np.loadtxt(savedir+'DatiBlacchino.txt')
ymio=deassorbi(Mio[0][puntidatogliere:],Mio[1][puntidatogliere:],dominguez)
dy=Mio[2][puntidatogliere:]
dyMiointrinsic=dy*deassorbi(Mio[0][puntidatogliere:],Mio[1][puntidatogliere:],dominguez)/Mio[1][puntidatogliere:]
ESyntrinsic=deassorbi(ESx,ESy,dominguez1es)
ESdyntrinsic=ESdy*deassorbi(ESx,ESy,dominguez1es)/ESy

parsmio,covm=curve_fit(PowerFit2,Mio[0][puntidatogliere:],ymio,sigma=dyMiointrinsic)
parsuo,comsuo=curve_fit(PowerFit2,ESx,ESyntrinsic,sigma=ESdyntrinsic)

plt.figure(figsize=(10,5))
plt.errorbar(Mio[0][puntidatogliere:],ymio,yerr=dyMiointrinsic,marker='o',linestyle='',label='Hess J1943+213 deabsorbed (z=0.21)' ,markersize=3,color='r',alpha=0.5)
plt.errorbar(ESx[:-1],ESyntrinsic[:-1],yerr=ESdyntrinsic[:-1],marker='v',linestyle='',label='1ES0229+200 deabsorbed (z=0.14)',markersize=3,color='b',alpha=0.5)
xrange=np.logspace(2,4, num=50)
#plt.plot(Mio[0][puntidatogliere:],PowerFit2(Mio[0][puntidatogliere:],parsmio[0],parsmio[1]),color='r')
plt.plot(xrange,PowerFit2(xrange,parsmio[0],parsmio[1]), color='r')
#plt.plot(ESx,PowerFit2(ESx[:-1],parsuo[0],parsuo[1]),color='b')
plt.plot(xrange,PowerFit2(xrange,parsuo[0],parsuo[1]), color='b')
plt.ylabel('Flux [GeV$^{-1}$cm$^{-2}$s$^{-1}$]')
plt.xlabel('Energy [GeV]')
plt.title('1ES0229+200 and HESS J1943+213 high energy deabsorbed flux comparison')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
plt.savefig(savedir+'FluxComparison1ES.eps', format='eps')