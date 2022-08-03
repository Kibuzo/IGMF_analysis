from scipy.special import expn
from scipy.special import erf
from scipy import integrate
import mpmath
import numpy as np
import matplotlib.pyplot as plt
import os

Cut=2080 #In GeV
#Cut=20000
Nzero=3.83e-14
nsample=int(1e6)                                            #Montecarlo samples
extent=[100,2e4]                                        #estremi in energia della PDF
gridpoints=int(extent[1]-extent[0])*10                  #numero di punti della griglia per la discretizzazione della CDF
estep=(extent[1]-extent[0])/gridpoints                  #passo della griglia
egrid=np.linspace(extent[0],extent[1],num=gridpoints)   #definizione della griglia

def Cutoff(E,N0,gamma,Ecut):
    return ((N0*(E/300.)**-gamma)*np.exp(-E/Ecut))
    
# partCDF restituisce l'integrale definito della pdf definita da N0, gamma, ecut tra le energie zero e uno. Lo faccio numerico che la versione analitica fa schifo
def partCDF(zero, uno, ecut, gamma, N0):
    return integrate.quad(Cutoff, zero, uno, args=(N0, gamma, ecut))[0]

def CDF(Emin, Emax, divs, N0, gamma, ecut):
    a=[]
    a.append(0)
    estep=float(Emax)/float(divs)
    for j in range (1,divs):
        a.append(a[-1]+partCDF(j*estep,(j+1)*estep,ecut,gamma,N0))
    return (a/np.max(a))
    
def CDFLOG(Emin, Emax, divs, N0, gamma, ecut):
    x=np.logspace(np.log10(Emin), np.log10(Emax), divs)
    a=[]
    a.append(0)
    for j in range (1, len(x)):
        a.append(a[-1]+(partCDF(x[j-1], x[j], ecut, gamma, N0)))
    return (a/np.max(a),x)


def searchlog (xi, CDF):
    idx=len(CDF[0][CDF[0]<xi])
    elzero=CDF[1][idx-1]
    eluno=CDF[1][idx]
    x=elzero+((xi-CDF[0][idx-1])/(CDF[0][idx] - CDF[0][idx-1])) * (eluno - elzero)
    return (x)
    
def randexpcutoff (a, cdf):
    return (search (np.random.random(), CDF))

def randexpcutofflog (a, cdf):
    return (searchlog (np.random.random(), cdf))
    
cutofflog=[]
cutoff=[]
cdflog=CDFLOG (100, 20000 , int(20000-1)/10 , Nzero, 1.5, 2080)
#cdf=CDF (5, 20000 , int(20000-5)*10 , 1, 1.5, 2080)
for j in range (0,nsample):
    cutofflog.append(searchlog(np.random.random(), cdflog))
    #cutoff.append(search(np.random.random(), cdf))

#plt.hist(cutoff, bins=50)
p=plt.hist(cutofflog,bins=50)#np.logspace(np.log10(100),np.log10(20000),50))

#p=[]
#for j in range (0,10000):
#    p.append(len(cdflog[cdflog<np.random.random()]))


# Parte di test: verifca che una pl con expcutoff vada d'accordo col rispettivo MC
bincentres=p[1][:-1]+(p[1][1:]-p[1][:-1])/2
rescaling=partCDF(100,20000,Cut,1.5,Nzero)/np.sum(p[0]*((p[1][1:]-p[1][:-1])))
plt.plot(bincentres,Cutoff(bincentres,Nzero,1.5,Cut)/rescaling)
plt.xscale('log')
plt.yscale('log')
plt.show()
    