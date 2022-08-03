from pylab import *
import os
import numpy as np
import time
import matplotlib.patches as mplarr
from mpl_toolkits.mplot3d import Axes3D
savedir='./BPL/Cone_benchmark/new_expcutoff/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
savedirtmp='./BPL/Cone_benchmark/new_expcutoff/TEMP/'

if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)
    
#definizione della funzione dello spettro intrinseco
def Cutoff(F0,E,Gamma,Ecut):
    return (F0*((E/Ezero)**-Gamma)*np.exp(-E/Ecut))

#Parametri della simulazione e dei file
nphot=1e4
d=887.6799277667532
Brms=3e-5
nfile=100
Ezero=300.
Nzero=3.83e-14
extent=[100,2e4]                                        #estremi in energia della PDF
gridpoints=int(extent[1]-extent[0])*10                  #numero di punti della griglia per la discretizzazione della CDF
estep=(extent[1]-extent[0])/gridpoints                  #passo della griglia
egrid=np.linspace(extent[0],extent[1],num=gridpoints)   #definizione della griglia

#Carica la sequenza di file: alla fine a è un vettore senza nessuna label ma con le colonne con indici identici a come escono da crpropa 
if os.path.exists(savedir+'final_largecone'+str(Brms)+'.txt'):
    a=np.loadtxt(savedir+'final_largecone'+str(Brms)+'.txt')
else:
    a=[]
    a=np.loadtxt(savedirtmp+'2data_final_largecone0'+str(Brms)+'.txt')
    length=len(np.unique(a[11]))
    for i in range (1,int(nfile-1)):
        filename=savedirtmp+'2data_final_largecone'+str(i)+str(Brms)+'.txt'
        if os.path.exists(filename):
            s=np.loadtxt(savedirtmp+'2data_final_largecone'+str(i)+str(Brms)+'.txt')
            length+=len(np.unique(s[11]))
            if s.shape[0]>0:
                if (a.ndim==1):
                    a=a[:, np.newaxis]
                if (s.ndim==1):
                    s=s[:, np.newaxis]
                a=np.concatenate((a,s),axis=1)
            else:
                print (filename + 'is empty, skipping...')
        else:
            print ('Warning: '+ filename +' does not exist')
print('A total of '+ str(length)+ ' primary photons generated the cascade photons') #così eviti muri del pianto


#Primissimo filtro: idx sono i fotoni, uguali e diversi sono fotoni che non hanno interagito o che
#uguali=(a[2]==a[11])
#diversi=~(a[2]==a[11])
cascadeidx=~(a[2]==a[11])
idx = np.abs(a[3]) == 22

#Rinomino un po' di cose
x0=a[14][idx*cascadeidx]
y0=a[15][idx*cascadeidx]
z0=a[16][idx*cascadeidx]
x1=a[23][idx*cascadeidx]
y1=a[24][idx*cascadeidx]
z1=a[25][idx*cascadeidx]
X=a[5][idx*cascadeidx]
Y=a[6][idx*cascadeidx]
Z=a[7][idx*cascadeidx]
Px=a[8][idx*cascadeidx]
Py=a[9][idx*cascadeidx]
Pz=a[10][idx*cascadeidx]
Px0=a[17][idx*cascadeidx]
Py0=a[18][idx*cascadeidx]
Pz0=a[19][idx*cascadeidx]

#Creo i vettori P tra cui fare il prodotto scalare
Deltavec=np.transpose((Px, Py, Pz))
Deltavec0=np.transpose((Px0, Py0, Pz0))
#arcoseno del prodotto scalare, somma dei prodotti elemento a elemento
Dot_product=np.arccos(np.clip(np.multiply(Deltavec, Deltavec0).sum(1), -1,1))
#normalizzazione della componente y
phi=y1/np.sqrt(y1**2+x1**2)
phi=np.array(np.arccos(phi))
#coseno e seno daranno x e y
cphi=np.sign(x1)*np.cos(np.array(phi))
sphi=np.sign(y1)*np.sqrt(1-cphi**2)
#Calcolo di lambdea e theta
lambdaxx=np.sqrt((x1-x0)**2+(z1-z0)**2+(y1-y0)**2)
Theta=np.arcsin((lambdaxx/d)*np.cos(np.pi/2-Dot_product))
Thetax=Theta*cphi
Thetay=Theta*sphi
#selezione degli eventi in theta
Hitpaolo=(np.sqrt((Theta*(180/np.pi))**2))<0.15

#Inizia il plot e il calcolo del chi²
intflusso=np.sum(estep*Cutoff(Nzero,egrid,1.5,2080))
scaling=nphot/intflusso

# --- Da qui in poi serve di avere i vettori dell'analisi Fermi e di HESS/VERITAS, però la sintassi dovrebbe essere autoesplicativa
cascade,bins=np.histogram(a[4][idx*cascadeidx][Hitpaolo]*1e9,bins=linbin)
cascade=cascade/scaling
bincenter = 0.5*(bins[1:]+bins[:-1])
intrinsic=Cutoff(Nzero,XSORT,1.5,2080)
cascade[11:]=0  #Taglia i fotoni fluttuanti a energie troppo alte
plt.figure(figsize=(10,5))
plt.plot(XSORT,intrinsic,label='Intrinsic emission model',linestyle='--')
plt.plot(XSORT,cascade,label='Cascade emission',linestyle='-.')
plt.plot(XSORT,cascade+intrinsic,label='Cascade +intrinsic emission')
plt.ylim(2e-17,3e-10)
plottapuntide()
plt.title('Observed vs simulated flux, 5 degree jet ''Brms='+str(Brms)+'nG')
plt.savefig(savedir+'fit_paolo'+str(Brms)+'.eps', format='eps')
plt.savefig(savedir+'fit_paolo'+str(Brms)+'.png')
print ('saved '+'fit_paolo'+str(Brms)+'.png in' + savedir )
plt.close()
residuals=cascade+intrinsic-YSORT
chifermi=calcolachi2de((cascade+intrinsic)[:-13],YSORT[:-13],DYSORT[:-13],XSORT[:-13])
print ('Chi quadro nella banda 0-300 GeV (metodo di Paolo)= '+str(chifermi))

