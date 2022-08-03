##Istruzioni per l'uso##
#Entra nella cartella Python di TesiMagistrale
#Esegui Il montecarlo fino alla definizione della savedir 
#Esegui FittaCherenki_new_mapchi_paolos fino a riga 444 e successivamente il blocco 466-486 dello stesso file -> outdated
#Ora puoi eseguire PlotBPL
from pylab import *
import os
import numpy as np
import time
import pandas as pd
from funct import *   #Achtung: va definito il redshift
fluxdir='/home/kibuzo/TesiMagistrale/Python/Flux_data/'
# # old directories
# savedir='/media/kibuzo/6317E4611EC4DD5F/TesiMagistrale/Python/BPL/'
# savedirtmp='/home/kibuzo/TesiMagistrale/Python/BPL/TEMP/'
# savedfig='/home/kibuzo/TesiMagistrale/Python/BPL/Figs/'

# new dir
savedir='/home/kibuzo/TesiMagistrale/Python/newdata/PointSource/BPL/'
savedfig='/home/kibuzo/TesiMagistrale/Python/newdata/PointSource/Figs/'
savedirtmp='/home/kibuzo/TesiMagistrale/Python/newdata/PointSource/Temp/'

if not os.path.exists(savedir):
    os.mkdir(savedir)
if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)
if not os.path.exists(savedfig):
    os.mkdir(savedfig)
nphot=1e5
CutoffTeV='2'
#CutoffTeV=str(4)



d=887.6799277667532 #Distanza in megaparsec
Brms=7e-6
#rms=6e-5
nfile=1000
activity=8 #tempo di attività della sorgente

def MpctoLy(dmpc):
    return (dmpc*3261564.78)

def countinside(radius, a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)

def profilo(radius,a):
    idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
    return np.sum(np.sqrt(a[6][idx]**2+a[5][idx]**2)<radius)
    
def calcolachi2de(o,e,s,x):
        return np.sum(((o-e)**2)/(s*deassorbi(x,e,dominguez)/e)**2)

combinati=pd.read_csv(fluxdir+'combined.csv')
XCOMBINED=np.transpose(combinati.as_matrix(columns=['XCOMBINED']))[0]
YCOMBINED=np.transpose(combinati.as_matrix(columns=['YCOMBINED']))[0]
DYCOMBINED=np.transpose(combinati.as_matrix(columns=['DYCOMBINED']))[0]

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
#         if (len(shape(b))<2):
#             b=b[:,newaxis]
#         if (len(a)==0):
#             a=b.copy()
#         else:
#             a=np.concatenate((a,b),axis=1)

#Lettura file: ricordati di rimuovere il file finale vecchio se usi parziali nuovi.
if os.path.exists(savedir+CutoffTeV+'final'+str(Brms)+'.txt'):
    a=np.loadtxt(savedir+CutoffTeV+'final'+str(Brms)+'.txt', dtype = np.float64)
else:
    a=[]
    index=0
    # check=True
    # while check:
    #     if os.path.exists(savedirtmp+CutoffTeV+'data_final_large'+str(index)+str(Brms)+'.txt'):
    #         a=np.loadtxt(savedirtmp+CutoffTeV+'data_final_large'+str(index)+str(Brms)+'.txt')
    #         if (len(shape(a))<2):
    #             a=a[:,newaxis]
    #         check=False
    #     else:
    #         index+=1
    a=np.loadtxt(savedir+CutoffTeV+'data_final_large0'+str(Brms)+'.txt', dtype = np.float64)
    for i in range (1,int(nfile-1)):
        filename=savedir+CutoffTeV+'data_final_large'+str(i)+str(Brms)+'.txt'
        if os.path.exists(filename):
            s=np.loadtxt(savedir+CutoffTeV+'data_final_large'+str(i)+str(Brms)+'.txt', dtype = np.float64)
            if (len(shape(s))<2):
                s=s[:,newaxis]
            if len(s)>0:
                a=np.concatenate((a,s),axis=1)
            #print (filename)
        else:
            print ('Warning: '+ filename +' does not exist or is an empty file')
            
# # The old way
# if os.path.exists(savedir+CutoffTeV+'final'+str(Brms)+'.txt'):
#     a=np.loadtxt(savedir+CutoffTeV+'final'+str(Brms)+'.txt')
# else:
#     a=[]
#     index=0
#     # check=True
#     # while check:
#     #     if os.path.exists(savedirtmp+CutoffTeV+'data_final_large'+str(index)+str(Brms)+'.txt'):
#     #         a=np.loadtxt(savedirtmp+CutoffTeV+'data_final_large'+str(index)+str(Brms)+'.txt')
#     #         if (len(shape(a))<2):
#     #             a=a[:,newaxis]
#     #         check=False
#     #     else:
#     #         index+=1
#     a=np.loadtxt(savedirtmp+CutoffTeV+'data_final'+str(Brms)+'.txt')
#     for i in range (1,int(nfile-1)):
#         filename=savedirtmp+CutoffTeV+'data_final'+str(i)+str(Brms)+'.txt'
#         if os.path.exists(filename):
#             s=np.loadtxt(savedirtmp+CutoffTeV+'data_final'+str(i)+str(Brms)+'.txt')
#             if (len(shape(s))<2):
#                 s=s[:,newaxis]
#             if len(s)>0:
#                 a=np.concatenate((a,s),axis=1)
#             #print (filename)
#         else:
#             print ('Warning: '+ filename +' does not exist or is an empty file')


#a=np.array(a)
vecs=[]
vettore=np.linspace(0,180,num=90)
for j in vettore:
    vecs.append(countinside(j,a))
plot(vecs,label='B='+str(Brms)+' ng')                     #Usare per BRMS diversi
print 'angolo massimo di emissione= ' +str(np.round(np.arcsin(np.max(np.sqrt(a[17]**2+a[18]**2)/np.sqrt(a[19]**2+a[17]**2+a[18]**2)))*180/np.pi,decimals=2))+ ' gradi'
theta=(np.sqrt(a[17]**2+a[18]**2)/np.sqrt(a[19]**2+a[17]**2+a[18]**2))*180/np.pi

px0=(a[17])
py0=(a[18])
pz0=(a[19])
x=px0/np.sqrt(px0**2+py0**2+pz0**2)*180/np.pi
y=py0/np.sqrt(px0**2+py0**2+pz0**2)*180/np.pi
z=pz0/np.sqrt(px0**2+py0**2+pz0**2)*180/np.pi
#plt.hist(1/(2*np.pi*tan(theta)))
plt.hist2d(x,y,bins=50)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Distribuzione dell angolo di emissione dei fotoni della sorgente')
plt.colorbar()
plt.savefig(savedfig+CutoffTeV+'photon_initial'+str(Brms)+'nG.png')
plt.close()

    
xlabel('Radius(Mpc)')
ylabel('Counts')
title('Number counts of photons as a fuction of radius')
legend(loc="upper left")
savefig(savedfig+CutoffTeV+'countsrad_all'+str(Brms)+'.png')
close()


idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni
theta=np.arcsin((a[5][idx])/np.sqrt((d-a[5][idx])**2+(d-a[6][idx])**2+(d-a[7][idx])**2))
phi=np.arcsin((a[6][idx])/np.sqrt((d-a[5][idx])**2+(d-a[6][idx])**2+(d-a[7][idx])**2))
theta=theta*360/np.pi
phi=phi*360/np.pi
plot(a[5][idx],a[6][idx],marker='.',markersize=1,linestyle='',label=str(Brms)+'ng')
xlabel('Theta [degrees]')
ylabel('Phi [degrees]')
title('Events scatter plot')
legend(loc="lower right")
savefig(savedfig+CutoffTeV+'PSF_'+str(Brms)+'ng'+'.png')              #Usare per BRMS diversi
close()
#hist2d(a[5][idx],a[6][idx],bins=100,cmap='jet')
hist2d(theta,phi,bins=100,norm=LogNorm())
colorbar()
xlabel('Theta [Degrees]')
ylabel('Phi [Degrees]')
title('Binned event map, B='+str(Brms)+'nG')
savefig(savedfig+CutoffTeV+'PSF_'+str(Brms)+'ng'+'_binned.png')       #Usare per BRMS diversi
close()

#a=np.loadtxt('./data_final.txt')

XSORT=np.sort(XCOMBINED) #Ricordati di prendere xcombined dal file FittaCherenki_new_mapchi_paolos.py
#scaling=np.loadtxt('./scaling.txt')
cascadeidx=~(a[2]==a[11])
MCset=np.loadtxt('./mc_complete.txt')
linbin=np.logspace(2,4.3,num=20)
linbin=np.copy(XBIN)
#binwidth=(np.array(XBIN[1:])-np.array(XBIN[:-1])) #Non serve!
intflusso=np.sum(estep*Cutoff(Nzero,egrid,1.5,2080))
influssothresh=np.sum(estep*Cutoff(Nzero,egrid[egrid>np.min(a[4][idx*~cascadeidx]*1e9)],1.5,2080))
scaling=nphot/intflusso #old, quando i fotoni sorgente erano sempre un numero fisso
#scaling=len(np.unique(a[20]))/intflusso #a[20] è la colonna degli id unici delle particelle in partenza
# #1) istogrammo l'intrinseco assorbito simulato da crpropa
# histo=plt.hist(a[4][idx*~cascadeidx]*1e9,bins=100)
# #2) trovo i bincenter
# histocenter=(histo[1][1:]+histo[1][:-1])/2
# #3) deassorbo il tutto usando come riferimento il bincenter
# deassorbitosim=np.sum(deassorbi(histocenter,histo[0],dominguez))
# #La somma dei contenuti dei bin è il flusso intrinseco ad alta energia della sorgente e lo posso usare per lo scaling.
# scaling=deassorbitosim/influssothresh
#MCset=MCset/scaling
px=(a[8])
py=(a[9])
pz=(a[10])
xarr=px/np.sqrt(px**2+py**2+pz**2)
yarr=py/np.sqrt(px**2+py**2+pz**2)
zarr=pz/np.sqrt(px**2+py**2+pz**2)
centbin=np.sqrt((np.arcsin(xarr)*180/np.pi+np.arcsin(yarr)*180/np.pi)**2)<0.15 #considero solo il bin centrale
cascade,bins=np.histogram(a[4][idx*cascadeidx*centbin]*1e9,bins=linbin)
#considero tutti i dati
#cascade,bins=np.histogram(a[4][idx*cascadeidx]*1e9,bins=linbin)
cascade=cascade/scaling
#cascade*=1e3
bincenter = 0.5*(bins[1:]+bins[:-1])
#intrinsic,bins=np.histogram(MCset,bins=linbin)
intrinsic=Cutoff(Nzero,XSORT,1.5,float(CutoffTeV)*1e3)
cascade[13:]=0  #Taglia i fotoni fluttuanti a energie troppo alte
plt.figure(figsize=(10,5))
plt.plot(XSORT,intrinsic,label='Intrinsic emission model',linestyle='--')
plt.plot(XSORT,cascade,label='Cascade emission',linestyle='-.')
plt.plot(XSORT,cascade+intrinsic,label='Cascade +intrinsic emission')
plt.ylim(2e-17,3e-10)
# plt.loglog(bincenter,intrinsic/scaling,label='Simulated spectrum')
# plt.loglog(bincenter,intrinsic/scaling+cascade,label='Simulated observed spectrum')
# plt.legend(loc='top right')
plottapuntide()
plt.title('Observed vs simulated flux, ''Brms='+str(Brms)+'nG')
plt.savefig(savedfig+CutoffTeV+'fit'+str(Brms)+'.eps', format='eps')
plt.savefig(savedfig+CutoffTeV+'fit'+str(Brms)+'.png')
plt.show()

plt.close()
contot=[]
conticerchio=[]
renorm=[]
sfera=0
xcount=np.linspace(1,30,30)
xalpha=np.arccos(xcount/30)
xbeta=np.arccos((xcount-1)/30)
h=30*(1-np.sin(xbeta))
xtheta=xbeta-xalpha
gusciosferico=0
for j in range (1,30):
    sfera+=gusciosferico#xtheta[j]*1800
    contaraggio=(countinside(j,a)-countinside(j-1,a))
    #gusciosferico=2*xtheta[j]*900
    gusciosferico=2*np.pi*30*(h[j]-h[j-1])
    gusciocilindrico=(np.pi*j**2-np.pi*(j-1)**2)
    renorm.append((gusciocilindrico/gusciosferico)*2)
    contot.append(contaraggio/gusciocilindrico)
    conticerchio.append(contaraggio/gusciosferico)
contot=np.array(contot)
plt.plot(xcount[:-1],conticerchio)
#plt.plot(xcount,contot, label=str(Brms))
plt.title('Event profile'+str(Brms))
plt.xlabel('radius [Mpc]')
plt.ylabel('Average counts')
plt.savefig(savedfig+CutoffTeV+'profile'+str(Brms)+'.png')
plt.close()

plt.hist2d(np.arcsin(xarr)*180/np.pi,np.arcsin(yarr)*180/np.pi,bins=1800,range=[[-90,90],[-90,90]],norm=LogNorm())
plt.clim(1,2*1e5)
plt.xlim(-4,4)
plt.ylim(-4,4)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction, 1 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.show()
plt.savefig(savedfig+CutoffTeV+'Impulsi'+str(Brms)+'.png')
plt.close()


# plt.hist2d(np.arcsin(xarr)*180/np.pi,np.arcsin(yarr)*180/np.pi,bins=23,range=[[-90,90],[-90,90]],norm=LogNorm())
# plt.clim(1,4*1e6)
# plt.xlim(-45,45)
# plt.ylim(-45,45)
# plt.xlabel('Theta [Degrees]')
# plt.ylabel('Phi [Degrees]')
# plt.title('Arrival direction, 4 degree bins, Brms='+str(Brms)+'nG')
# plt.colorbar()
# plt.savefig(savedir+CutoffTeV+'Impulsi_binlarghi'+str(Brms)+'.png')
# plt.close()

residuals=cascade+intrinsic-YSORT
chi=calcolachi2de(cascade+intrinsic,YSORT,DYSORT,XSORT)
chifermi=calcolachi2de((cascade+intrinsic)[:-13],YSORT[:-13],DYSORT[:-13],XSORT[:-13])
print('Risultati del fit per un campo magnetico di '+str(Brms)+' nG')
print('Chi quadro completo= '+str(chi))
print ('Chi quadro nella banda 5-300 GeV= '+str(chifermi))
print('\n')

delayed=MpctoLy(a[0]-np.min(a[0]))<activity
cascadelay,bins=np.histogram(a[4][idx*cascadeidx*centbin*delayed]*1e9,bins=linbin)
cascadelay=cascadelay/scaling
plt.figure(figsize=(10,5))
plt.plot(XSORT,intrinsic,label='Intrinsic emission model',linestyle='--')
plt.plot(XSORT,cascadelay,label='Cascade emission',linestyle='-.')
plt.plot(XSORT,cascadelay+intrinsic,label='Cascade +intrinsic emission')
plt.ylim(2e-17,3e-10)
plottapuntide()
plt.title('Observed vs simulated flux [Delay<'+"8"+'yr], ''Brms='+str(Brms)+'nG')
plt.savefig(savedfig+CutoffTeV+'fit_delay'+str(Brms)+'.eps', format='eps')
plt.savefig(savedfig+CutoffTeV+'fit_delay'+str(Brms)+'.png')
plt.show()
plt.close()

residuals=cascadelay+intrinsic-YSORT
chi=calcolachi2de(cascadelay+intrinsic,YSORT,DYSORT,XSORT)
chifermi=calcolachi2de((cascadelay+intrinsic)[:-13],YSORT[:-13],DYSORT[:-13],XSORT[:-13])
print('Risultati del fit per un campo magnetico di '+str(Brms)+' nG e tempo di ritardo di '+str(activity)+' anni')
print('Chi quadro completo= '+str(chi))
print ('Chi quadro nella banda 5-300 GeV= '+str(chifermi))


# binspace=np.logspace(np.log10(10),np.log10(2000),num=20)
# plt.hist(a[4][idx*cascadeidx]*1e9,bins=binspace)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# Salva il file completo a meno che non esista già
if os.path.exists(savedir+CutoffTeV+'final'+str(Brms)+'.txt'):
    print ('File already exists, doing nothing')
else:
    np.savetxt(savedir+CutoffTeV+'final'+str(Brms)+'.txt',a)



