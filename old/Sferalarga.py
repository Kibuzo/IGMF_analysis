##Istruzioni per l'uso##
#Esegui Il montecarlo fino alla definizione della savedir
#Esegui FittaCherenki_new_mapchi_paolos fino a riga 444 e successivamente il blocco 466-486 dello stesso file
#Ora puoi eseguire PlotBPL
from pylab import *
import os
import numpy as np
import time
import pandas as pd
savedir='./BPL/Cono/filtered/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
savedirtmp='./BPL/Cono/Data/'
if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)
BinTeV=0.1 #Dimensione del bin della misura al TeV
Opening=5 #Apertura del getto
nphot=1e6*(BinTeV/Opening) #Numero di fotoni inizialmente "on target" su cui fare la normalizzazione
CutoffTeV='500Mpc'
CutoffTeVe='2'

def loadfile(l):
    l=str(l)
    z=[]
    print(l)
    z=pd.read_csv(savedirtmp+l+'_data_final_filtered0'+str(Brms)+'.txt', header= None)
    for i in range (1,int(nfile-1)):
        filename=savedirtmp+l+'_data_final_filtered'+str(i)+str(Brms)+'.txt'
        if os.path.exists(filename):
            s_z=pd.read_csv(savedirtmp+l+'_data_final_filtered'+str(i)+str(Brms)+'.txt', header=None)
            z=pd.concat([z,s_z], axis=0)
        else:
            print ('Warning: '+ filename +' does not exist')
    #np.savetxt(savedir+l+'_filtered_final'+str(Brms)+'.txt',z)
    return(z)


#CutoffTeV=str(4)


d=887.6799277667532 #Distanza in megaparsec
Brms=1e-07
#Brms=int(0)
nfile=10000

def countinside(radius, a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)

def profilo(radius,a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)
    
def calcolachi2de(o,e,s,x):
        return np.sum(((o-e)**2)/(s*deassorbi(x,e,dominguez)/e)**2)

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

if os.path.exists(savedir+'pd_data_filtered_final'+str(Brms)+'.csv'):
    data=pd.read_csv('output_list.txt', sep=" ", header=None)
    data.columns=["IDX", "ID0", "ID1", "x", "y", "z", "px", "py", "pz", "E"]
else:
    data=loadfile('idx')
    data['ID0']=loadfile('id0')
    data['ID1']=loadfile('id1')
    data['x']=loadfile('x')
    data['y']=loadfile('y')
    data['z']=loadfile('z')
    data['px']=loadfile('px')
    data['py']=loadfile('py')
    data['pz']=loadfile('pz')
    data['E']=loadfile('E')
    data.columns=["IDX", "ID0", "ID1", "x", "y", "z", "px", "py", "pz", "E"]



detectorsizedeg=0.1
hitx = np.abs(data['x'])<detectorsizedeg*(867*np.tan(2*np.pi/180))
hity = np.abs(data['y'])<detectorsizedeg*(867*np.tan(2*np.pi/180))
hitz = np.abs(data['z'])<detectorsizedeg*(867*np.tan(2*np.pi/180))
photon=np.abs(data['IDX'])==22
halo=~(data['ID0']==data['ID1'])
source=(data['ID0']==data['ID1'])
hit=np.sqrt(np.abs(data['x'])**2+np.abs(data['y'])**2+np.abs(data['z'])**2)<detectorsizedeg*(867*np.tan(2*np.pi/180))
binwidthdeg=0.1
data_hit=data[photon*hit]
data_halo=data[photon*hit*halo]
data_src=data[photon*hit*source]

#plot halo+source
xarr=data_hit['px']/np.sqrt(data_hit['py']**2+data_hit['pz']**2)
yarr=data_hit['py']/np.sqrt(data_hit['px']**2+data_hit['pz']**2)
zarr=data_hit['pz']/np.sqrt(data_hit['px']**2+data_hit['py']**2)
theta=np.arctan(data_hit['px']/np.sqrt(data_hit['py']**2+data_hit['pz']**2))#*(np.pi/180)
phi=np.arctan(data_hit['pz']/data_hit['py'])#*(np.pi/180)
rng=10/binwidthdeg
plt.hist2d(np.arctan(xarr)*180/np.pi,np.arctan(yarr)*180/np.pi,bins=2*rng,range=[[-9.9,9.9],[-9.9,9.9]],cmin=1, norm=LogNorm())#,norm=LogNorm())
#plt.scatter(np.arctan(xarr)*180/np.pi, np.arctan(yarr)*180/np.pi,s=0.1)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction, 0.1 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.show()
plt.savefig(savedir+'PSF_large_all'+str(Brms)+'.png')
plt.close()

#plot halo
xarr=data_halo['px']/np.sqrt(data_halo['py']**2+data_halo['pz']**2)
yarr=data_halo['py']/np.sqrt(data_halo['px']**2+data_halo['pz']**2)
zarr=data_halo['pz']/np.sqrt(data_halo['px']**2+data_halo['py']**2)
theta=np.arctan(data_halo['px']/np.sqrt(data_halo['py']**2+data_halo['pz']**2))#*(np.pi/180)
phi=np.arctan(data_halo['pz']/data_halo['py'])#*(np.pi/180)
rng=10/binwidthdeg
plt.hist2d(np.arctan(xarr)*180/np.pi,np.arctan(yarr)*180/np.pi,bins=2*rng,range=[[-9.9,9.9],[-9.9,9.9]],cmin=1, norm=LogNorm())#,norm=LogNorm())
#plt.scatter(np.arctan(xarr)*180/np.pi, np.arctan(yarr)*180/np.pi,s=0.1)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction, 0.1 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.show()
plt.savefig(savedir+'PSF_large_halo'+str(Brms)+'.png')
plt.close()

#plot source
xarr=data_src['px']/np.sqrt(data_src['py']**2+data_src['pz']**2)
yarr=data_src['py']/np.sqrt(data_src['px']**2+data_src['pz']**2)
zarr=data_src['pz']/np.sqrt(data_src['px']**2+data_src['py']**2)
theta=np.arctan(data_src['px']/np.sqrt(data_src['py']**2+data_src['pz']**2))#*(np.pi/180)
phi=np.arctan(data_src['pz']/data_src['py'])#*(np.pi/180)
rng=10/binwidthdeg
plt.hist2d(np.arctan(xarr)*180/np.pi,np.arctan(yarr)*180/np.pi,bins=2*rng,range=[[-9.9,9.9],[-9.9,9.9]],cmin=1, norm=LogNorm())#,norm=LogNorm())
#plt.scatter(np.arctan(xarr)*180/np.pi, np.arctan(yarr)*180/np.pi,s=0.1)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction, 0.1 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.show()
plt.savefig(savedir+'PSF_large_source'+str(Brms)+'.png')
plt.close()


# Salva il file completo a meno che non esista già
if os.path.exists(savedir+CutoffTeV+'final'+str(Brms)+'.txt'):
    print ('File already exists, doing nothing')
else:
    data.to_csv(savedir+'pd_data_filtered_final'+str(Brms)+'.csv')
    #np.savetxt(savedir+CutoffTeV+'final'+str(Brms)+'.txt',a)
