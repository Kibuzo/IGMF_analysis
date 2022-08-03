##Istruzioni per l'uso##
#Esegui Il montecarlo fino alla definizione della savedir
#Esegui FittaCherenki_new_mapchi_paolos fino a riga 444 e successivamente il blocco 466-486 dello stesso file
#Ora puoi eseguire PlotBPL
from pylab import *
import os
import numpy as np
import time
savedir='./BPL/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
savedirtmp='./BPL/TEMP/'
if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)
nphot=1e5-1e2
CutoffTeV=str(4)

d=887.6799277667532 #Distanza in megaparsec
redshift=0.21
#Brms=1e-09
Brms=int(0)
nfile=100
dominguez = Absorption.read_builtin('dominguez').table_model(redshift)

def Gettau(energy,model):
    return(-np.log(model.evaluate(energy*u.GeV,1)))

def countinside(radius, a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)

def profilo(radius,a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)
    
def calcolachi2de(o,e,s,x):
        return np.sum(((o-e)**2)/(s*deassorbi(x,e,dominguez)/e)**2)
    
def Ritcascade(flusso,Tdelay):
    return(flusso*(10/(10+Tdelay)))

def TtoB (T,z,tau,E100):
    T0=22.08e-12
    B18sq=(T/T0)*((1+z)**5)*(E100**(5./2.))/(1-tau**-1)
    Bng=np.sqrt(B18sq)*1e-09

def EtoT (E,z,tau,BnG):
    kappa=1
    B18=BnG*1e09
    #T0=22.08e+12
    T0=7e05/3.154e07
    E100=E/100
    T=T0*kappa*(1-tau**-1)*((1+z)**-5)*(E100**-(5./2.))*B18**2
    return(T)

XSORT=[XCOMBINED[0],XCOMBINED[1],XCOMBINED[2],XCOMBINED[3],XCOMBINED[4],XCOMBINED[5],XCOMBINED[6],XCOMBINED[7],XCOMBINED[12],XCOMBINED[8],XCOMBINED[13],XCOMBINED[9],XCOMBINED[14],XCOMBINED[10],XCOMBINED[15],XCOMBINED[16],XCOMBINED[11],XCOMBINED[17],XCOMBINED[18],XCOMBINED[19]]
YSORT=[YCOMBINED[0],YCOMBINED[1],YCOMBINED[2],YCOMBINED[3],YCOMBINED[4],YCOMBINED[5],YCOMBINED[6],YCOMBINED[7],YCOMBINED[12],YCOMBINED[8],YCOMBINED[13],YCOMBINED[9],YCOMBINED[14],YCOMBINED[10],YCOMBINED[15],YCOMBINED[16],YCOMBINED[11],YCOMBINED[17],YCOMBINED[18],YCOMBINED[19]]
DYSORT=[DYCOMBINED[0],DYCOMBINED[1],DYCOMBINED[2],DYCOMBINED[3],DYCOMBINED[4],DYCOMBINED[5],DYCOMBINED[6],DYCOMBINED[7],DYCOMBINED[12],DYCOMBINED[8],DYCOMBINED[13],DYCOMBINED[9],DYCOMBINED[14],DYCOMBINED[10],DYCOMBINED[15],DYCOMBINED[16],DYCOMBINED[11],DYCOMBINED[17],DYCOMBINED[18],DYCOMBINED[19]]
YSORT=np.array(YSORT)
DYSORT=np.array(DYSORT)
XSORT=np.array(XSORT)

YSORT=deassorbi(XSORT,YSORT,dominguez)
DYSORT=DYSORT*deassorbi(XSORT,YSORT,dominguez)/YSORT
def Cutoff(F0,E,Gamma,Ecut):
    return (F0*((E/Ezero)**-Gamma)*np.exp(-E/Ecut))

# #Crea la lista a coi dati partendo dal primo
# filename0=savedirtmp+'photon_electron_output'+str(Brms)+'0.txt'
# if not (os.stat(filename0).st_size==0):
#     a=np.loadtxt(filename0,unpack=True)
# 
# a=[]
# for i in range (0,nphot):
#     filename=savedirtmp+'photon_electron_output'+str(Brms)+str(i)+'.txt' #Usare per BRMS uguali
#     if not (os.stat(filename).st_size==0):
#         b=loadtxt(filename,unpack=True)
#         if (len(shape(b))<2):##Istruzioni per l'uso##
#Esegui Il montecarlo fino alla definizione della savedir
#Esegui FittaCherenki_new_mapchi_paolos fino a riga 444 e successivamente il blocco 466-486 dello stesso file
#Ora puoi eseguire PlotBPL
from pylab import *
import os
import numpy as np
import time
savedir='./BPL/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
savedirtmp='./BPL/TEMP/'
if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)
#nphot=1e5-1e3

d=887.6799277667532 #Distanza in megaparsec
Brms=int(0)
nfile=100

def countinside(radius, a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)

def profilo(radius,a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)
    
def calcolachi2de(o,e,s,x):
        return np.sum(((o-e)**2)/(s*deassorbi(x,e,dominguez)/e)**2)

# YSORT=[YCOMBINED[0],YCOMBINED[1],YCOMBINED[2],YCOMBINED[3],YCOMBINED[4],YCOMBINED[5],YCOMBINED[6],YCOMBINED[7],YCOMBINED[12],YCOMBINED[8],YCOMBINED[13],YCOMBINED[9],YCOMBINED[14],YCOMBINED[11],YCOMBINED[10],YCOMBINED[15],YCOMBINED[16],YCOMBINED[17],YCOMBINED[18],YCOMBINED[19]]
# DYSORT=[DYCOMBINED[0],DYCOMBINED[1],DYCOMBINED[2],DYCOMBINED[3],DYCOMBINED[4],DYCOMBINED[5],DYCOMBINED[6],DYCOMBINED[7],DYCOMBINED[12],DYCOMBINED[8],DYCOMBINED[13],DYCOMBINED[9],DYCOMBINED[14],DYCOMBINED[11],DYCOMBINED[10],DYCOMBINED[15],DYCOMBINED[16],DYCOMBINED[17],DYCOMBINED[18],DYCOMBINED[19]]
# YSORT=np.array(YSORT)
# DYSORT=np.array(DYSORT)

XSORT=[XCOMBINED[0],XCOMBINED[1],XCOMBINED[2],XCOMBINED[3],XCOMBINED[4],XCOMBINED[5],XCOMBINED[6],XCOMBINED[7],XCOMBINED[12],XCOMBINED[8],XCOMBINED[13],XCOMBINED[9],XCOMBINED[14],XCOMBINED[10],XCOMBINED[15],XCOMBINED[16],XCOMBINED[11],XCOMBINED[17],XCOMBINED[18],XCOMBINED[19]]
YSORT=[YCOMBINED[0],YCOMBINED[1],YCOMBINED[2],YCOMBINED[3],YCOMBINED[4],YCOMBINED[5],YCOMBINED[6],YCOMBINED[7],YCOMBINED[12],YCOMBINED[8],YCOMBINED[13],YCOMBINED[9],YCOMBINED[14],YCOMBINED[10],YCOMBINED[15],YCOMBINED[16],YCOMBINED[11],YCOMBINED[17],YCOMBINED[18],YCOMBINED[19]]
DYSORT=[DYCOMBINED[0],DYCOMBINED[1],DYCOMBINED[2],DYCOMBINED[3],DYCOMBINED[4],DYCOMBINED[5],DYCOMBINED[6],DYCOMBINED[7],DYCOMBINED[12],DYCOMBINED[8],DYCOMBINED[13],DYCOMBINED[9],DYCOMBINED[14],DYCOMBINED[10],DYCOMBINED[15],DYCOMBINED[16],DYCOMBINED[11],DYCOMBINED[17],DYCOMBINED[18],DYCOMBINED[19]]
YSORT=np.array(YSORT)
DYSORT=np.array(DYSORT)
XSORT=np.array(XSORT)

YSORT=deassorbi(XSORT,YSORT,dominguez)
DYSORT=DYSORT*deassorbi(XSORT,YSORT,dominguez)/YSORT
def Cutoff(F0,E,Gamma,Ecut):
    return (F0*((E/Ezero)**-Gamma)*np.exp(-E/Ecut))

#Lettura file: ricordati di rimuovere il file finale vecchio se usi parziali nuovi.
if os.path.exists(savedir+CutoffTeV+'final'+str(Brms)+'.txt'):
    a=np.loadtxt(savedir+CutoffTeV+'final'+str(Brms)+'.txt')
else:
    a=[]
    a=np.loadtxt(savedirtmp+CutoffTeV+'data_final0'+str(Brms)+'.txt')
    for i in range (1,int(nfile-1)):
        filename=savedirtmp+CutoffTeV+'data_final'+str(i)+str(Brms)+'.txt'
        if os.path.exists(filename):
            s=np.loadtxt(savedirtmp+CutoffTeV+'data_final'+str(i)+str(Brms)+'.txt')
            a=np.concatenate((a,s),axis=1)
        else:
            print ('Warning: '+ CutoffTeV+filename +' does not exist')


#Creazione dei vettori di cascade scalati per il numero di fotoni generati dal montecarlo
XSORT=np.sort(XCOMBINED) #Ricordati di prendere xcombined dal file FittaCherenki_new_mapchi_paolos.py
cascadeidx=~(a[2]==a[11])
idx=np.abs(a[3])==22
MCset=np.loadtxt('./mc_complete.txt')
linbin=np.logspace(2,4.3,num=20)
linbin=np.copy(XBIN)
intflusso=np.sum(estep*Cutoff(Nzero,egrid,1.5,2080))
scaling=nphot/intflusso
cascade,bins=np.histogram(a[4][idx*cascadeidx]*1e9,bins=linbin)
cascade=cascade/scaling

# for j in range (0,len(cascade)-1):
#     cascaderit=Ritcascade(cascade[j],delay)

ExtremeDelay=[]
BestDelay=[]
for i in range (0,6):
    cascaderit=[]
    Chi=[]
    for delay in range (0,int(1e3)):
        cascaderit.append(Ritcascade(cascade[i],float(delay)))
        Chi.append(((cascaderit[-1]-YSORT[i])**2)/DYSORT[i]**2)
    BestDelay.append(np.where(Chi==np.min(Chi))[0][0])
    ExtremeDelay.append(np.where(np.array(Chi)>=BestDelay[-1]+2.7)[0][0])

def plottacurva(BnG):
    E=np.linspace(0,1000,num=1000)
    eprim=np.sqrt(20/0.32)*E*1e3
    T=EtoT(np.array(E),redshift,Gettau(eprim,dominguez),BnG)
    plt.plot(E,T,label='B= '+str(BnG)+'nG')


Extremeulims=np.array(ExtremeDelay)/1.5
lowlims=[1,1,1,1,1,1]
loXERR=[XSORT[1]-XSORT[0],XSORT[1]-XSORT[0],XSORT[2]-XSORT[1],XSORT[3]-XSORT[2],XSORT[4]-XSORT[3],XSORT[5]-XSORT[4]]
upXERR=[XSORT[1]-XSORT[0],XSORT[2]-XSORT[1],XSORT[3]-XSORT[2],XSORT[4]-XSORT[3],XSORT[5]-XSORT[4],XSORT[6]-XSORT[5]]
loXERR=np.array(loXERR)/2
upXERR=np.array(upXERR)/2
plottacurva(1e-08)
plottacurva(1e-07)
plottacurva(1e-06)
#plt.plot(XSORT[:6],ExtremeDelay,linestyle='none',marker='+', label='Lower limits')
plt.errorbar(XSORT[:6],ExtremeDelay,yerr=Extremeulims,xerr=(loXERR,upXERR),capsize=2,linestyle='none',marker='.', lolims=Extremeulims, label='Lower limits')

plt.xlabel('Energy [GeV]')
plt.ylabel('Time delay [years]')
#plt.ylim(np.log(1e-1),5)
#plt.xlim=(1,3)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()
            
    
