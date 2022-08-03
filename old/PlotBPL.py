##Istruzioni per l'uso##
#Esegui Il montecarlo fino alla definizione della savedir
#Esegui FittaCherenki_new_mapchi_paolos fino a riga 467 
#Ora puoi eseguire PlotBPL
from pylab import *
import os
import numpy as np
import time
import matplotlib.patches as mplarr
from mpl_toolkits.mplot3d import Axes3D
#savedir='./BPL/Cone_final/'
savedir='./BPL/Cone_benchmark/200GeV/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
#savedirtmp='./BPL/Cone_final/TEMP/'
savedirtmp='./BPL/Cone_benchmark/200GeV/TEMP/'

if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)
BinTeV=0.1**2 #Dimensione del bin della misura al TeV
Opening=5**2 #Apertura del getto
# Un detaglio sul numero di fotoni: ogni simulazione del cnaf lancia tutti i fotoni che trova nel file, ogni file è composto da 100 fotoni, pertanto se lanci 1000 simulazioni otterrai 1000*100 fotoni ovvero 1e5. Il totale di fotoni disponibili nel montecarlo è di 1e6.
#nphot=1e5*(BinTeV/Opening) #Numero di fotoni inizialmente "on target" su cui fare la normalizzazione
nphot=1e4

#Returna l'angolo solido da moltiplicare per pi in funzione dell'apertura (non semiapertura) del cono in gradi
def solid (alphadeg):
    alphamezzi=alphadeg/2
    return (2*(1-np.cos(alphamezzi*(np.pi/180))))
    
Cut=2080 #In GeV
#Cut=20000
Nzero=3.83e-14
extent=[200,2e4]                                        #estremi in energia della PDF
gridpoints=int(extent[1]-extent[0])*10                  #numero di punti della griglia per la discretizzazione della CDF
estep=(extent[1]-extent[0])/gridpoints                  #passo della griglia
egrid=np.linspace(extent[0],extent[1],num=gridpoints)   #definizione della griglia

d=887.6799277667532 #Distanza in megaparsec
# dcm=d*3.086e+24 #distanza in cm

# solidarea=solid(5)*np.pi*dcm**2
# nphot=nphot/solidarea
#Lv=solid(5)*np.pi*intflusso*dcm**2
#Brms=3e-5
Brms=int(0)
#Achtung: al momento mancano 3 file nel campione 1e-5, sono 100 negli altri
nfile=100
Ezero=300. #in GeV, energia di decorrelazione.

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
#if os.path.exists(savedir+'final_largecone'+str(Brms)+'.txt'):
#    a=np.loadtxt(savedir+'final_largecone'+str(Brms)+'.txt')
if True:
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
print('A total of '+ str(length)+ ' primary photons generated the cascade that hit the detector')

#Indici dei fotoni primari e dei fotoni solo cascade
uguali=(a[2]==a[11])
diversi=~(a[2]==a[11])
idx = np.abs(a[3]) == 22 #Seleziona gli indici che contengono fotoni



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
plt.savefig(savedir+'photon_initial'+str(Brms)+'nG.png')
plt.close()

    
xlabel('Radius(Mpc)')
ylabel('Counts')
title('Number counts of photons as a fuction of radius')
legend(loc="upper left")
savefig(savedir+'countsrad_all'+str(Brms)+'.png')
close()


# theta=np.arcsin((a[5][idx])/np.sqrt((d-a[5][idx])**2+(d-a[6][idx])**2+(d-a[7][idx])**2))
# phi=np.arcsin((a[6][idx])/np.sqrt((d-a[5][idx])**2+(d-a[6][idx])**2+(d-a[7][idx])**2))
# theta=theta*360/np.pi
# phi=phi*360/np.pi
# plot(a[5][idx],a[6][idx],marker='.',markersize=1,linestyle='',label=str(Brms)+'ng')
# xlabel('Theta [degrees]')
# ylabel('Phi [degrees]')
# title('Events scatter plot')
# legend(loc="lower right")
# savefig(savedir+'PSF_'+str(Brms)+'ng'+'.png')              #Usare per BRMS diversi
# close()
# #hist2d(a[5][idx],a[6][idx],bins=100,cmap='jet')
# hist2d(theta,phi,bins=100,cmap='jet',norm=LogNorm())
# colorbar()
# xlabel('Theta [Degrees]')
# ylabel('Phi [Degrees]')
# title('Binned event map, B='+str(Brms)+'nG')
# savefig(savedir+'PSF_'+str(Brms)+'ng'+'_binned.png')       #Usare per BRMS diversi
# close()

#a=np.loadtxt('./data_final.txt')

acceptance=3.
halfopening=degrees(np.arctan(np.max(np.sqrt(px0**2+py0**2))))

px=(a[8])
py=(a[9])
pz=(a[10])
xarr=px/np.sqrt(px**2+py**2+pz**2)
yarr=py/np.sqrt(px**2+py**2+pz**2)
zarr=pz/np.sqrt(px**2+py**2+pz**2)
dentrofermi=np.sqrt((np.arcsin(xarr)*180/np.pi)**2+(np.arcsin(yarr)*180/np.pi)**2)<0.3
XSORT=np.sort(XCOMBINED) #Ricordati di prendere xcombined dal file FittaCherenki_new_mapchi_paolos.py
#scaling=np.loadtxt('./scaling.txt')
cascadeidx=~(a[2]==a[11])
#MCset=np.loadtxt('./mc_complete.txt')
linbin=np.logspace(2,4.3,num=20)
linbin=np.copy(XBIN)
#binwidth=(np.array(XBIN[1:])-np.array(XBIN[:-1])) #Non serve!
intflusso=np.sum(estep*Cutoff(Nzero,egrid,1.5,2080))
scaling=length/intflusso
if (halfopening>acceptance/2):
    smallarea=np.pi*(acceptance/2.)**2
    largearea=np.pi*halfopening**2
    scaling*=(smallarea/largearea) #The effective simulated photons (hitting the acceptance when undeflected) are less than the simulated.
#MCset=MCset/scaling
cascade,bins=np.histogram(a[4][idx*cascadeidx*dentrofermi]*1e9,bins=linbin)
cascade=cascade/scaling
#cascade*=1e3
bincenter = 0.5*(bins[1:]+bins[:-1])
#intrinsic,bins=np.histogram(MCset,bins=linbin)
intrinsic=Cutoff(Nzero,XSORT,1.5,2080)
cascade[11:]=0  #Taglia i fotoni fluttuanti a energie troppo alte
plt.figure(figsize=(10,5))
plt.plot(XSORT,intrinsic,label='Intrinsic emission model',linestyle='--')
plt.plot(XSORT,cascade,label='Cascade emission',linestyle='-.')
plt.plot(XSORT,cascade+intrinsic,label='Cascade +intrinsic emission')
plt.ylim(2e-17,3e-10)
# plt.loglog(bincenter,intrinsic/scaling,label='Simulated spectrum')
# plt.loglog(bincenter,intrinsic/scaling+cascade,label='Simulated observed spectrum')
# plt.legend(loc='top right')
plottapuntide()
plt.title('Observed vs simulated flux, 5 degree jet ''Brms='+str(Brms)+'nG')
plt.savefig(savedir+'fit'+str(Brms)+'.eps', format='eps')
plt.savefig(savedir+'fit'+str(Brms)+'.png')
print ('saved '+'fit'+str(Brms)+'.png in' + savedir )
plt.show()
plt.close()

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
plt.savefig(savedir+'profile'+str(Brms)+'.png')
print ('saved '+'profile'+str(Brms)+'.png in' + savedir )
plt.close()


plt.hist2d(np.arctan(xarr)*180/np.pi,np.arctan(yarr)*180/np.pi,bins=1800,range=[[-90,90],[-90,90]],cmin=1)#,norm=LogNorm())
#plt.clim(1,2*1e5)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction (source+cascade), 0.1 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.show()
plt.savefig(savedir+'Impulsi_tot'+str(Brms)+'.png')
print ('saved '+'Impulsi_tot'+str(Brms)+'.png in' + savedir )

plt.close()

px=(a[8][idx*uguali])
py=(a[9][idx*uguali])
pz=(a[10][idx*uguali])
xarr=px/np.sqrt(px**2+py**2+pz**2)
yarr=py/np.sqrt(px**2+py**2+pz**2)
zarr=pz/np.sqrt(px**2+py**2+pz**2)
plt.hist2d(np.arcsin(xarr)*180/np.pi,np.arcsin(yarr)*180/np.pi,bins=1800,range=[[-90,90],[-90,90]],cmin=1)#,norm=LogNorm())
#plt.clim(1,2*1e5)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction (source only), 0.1 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.show()
plt.savefig(savedir+'Impulsi_source'+str(Brms)+'.png')
print ('saved '+'Impulsi_source'+str(Brms)+'.png in' + savedir )

plt.close()

px=(a[8][idx*cascadeidx])
py=(a[9][idx*cascadeidx])
pz=(a[10][idx*cascadeidx])
xarr=px/np.sqrt(px**2+py**2+pz**2)
yarr=py/np.sqrt(px**2+py**2+pz**2)
zarr=pz/np.sqrt(px**2+py**2+pz**2)
dentrozerotre=np.sqrt(np.arcsin(xarr)**2+np.arcsin(yarr)**2)*180/np.pi<0.15
plt.hist2d(np.arcsin(xarr)*180/np.pi,np.arcsin(yarr)*180/np.pi,bins=1800,range=[[-90,90],[-90,90]],cmin=1)#,norm=LogNorm())
#plt.clim(1,2*1e5)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction (cascade), 0.1 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.show()
plt.savefig(savedir+'Impulsi_halo'+str(Brms)+'.png')
print ('saved '+'Impulsi_halo'+str(Brms)+'.png in' + savedir )
plt.close()


plt.hist2d(np.arcsin(xarr)*180/np.pi,np.arcsin(yarr)*180/np.pi,bins=23,range=[[-90,90],[-90,90]],norm=LogNorm())
plt.clim(1,4*1e6)
plt.xlim(-45,45)
plt.ylim(-45,45)
plt.xlabel('Theta [Degrees]')
plt.ylabel('Phi [Degrees]')
plt.title('Arrival direction, 4 degree bins, Brms='+str(Brms)+'nG')
plt.colorbar()
plt.savefig(savedir+'Impulsi_binlarghi'+str(Brms)+'.png')
plt.close()

residuals=cascade+intrinsic-YSORT
chi=calcolachi2de(cascade+intrinsic,YSORT,DYSORT,XSORT)
chifermi=calcolachi2de((cascade+intrinsic)[:-13],YSORT[:-13],DYSORT[:-13],XSORT[:-13])
print('Risultati del fit per un campo magnetico di '+str(Brms)+' nG')
print('Chi quadro completo= '+str(chi))
print ('Chi quadro nella banda 0-300 GeV= '+str(chifermi))

print ('\n Paolo fit begins')
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


DeltaY=np.transpose((Py, Pz))
DeltaX=np.transpose((Px, Pz))
DeltaY0=np.transpose((Py0, Pz0))
DeltaX0=np.transpose((Px0, Pz0))
Deltavec=np.transpose((Px, Py, Pz))
Deltavec0=np.transpose((Px0, Py0, Pz0))
# # old and slow for theta and phi, left for debugging
# Dot_product=[]
# phi=[]
# for j in range (0,len(DeltaX0)):
#     Dot_product.append(np.arccos(np.clip(np.dot(Deltavec0[j], Deltavec[j]), -1.0,1.0)))
#     phi.append((np.clip(np.dot(((x1[j],y1[j],0)/np.linalg.norm([x1[j],y1[j],0])), (0,1,0)), -1.0,1.0)))

# the new fast way for theta and phi
# Pnorm=np.sqrt(np.multiply(Deltavec,Deltavec).sum(1))
# P0norm=np.sqrt(np.multiply(Deltavec0,Deltavec0).sum(1))
#In realtà sono già normalizzati in crpropa, quindi skippo
Dot_product=np.arccos(np.clip(np.multiply(Deltavec, Deltavec0).sum(1), -1,1))
phi=y1/np.sqrt(y1**2+x1**2)

phi=np.array(np.arccos(phi))
cphi=np.sign(x1)*np.cos(np.array(phi))
sphi=np.sign(y1)*np.sqrt(1-cphi**2)
Dot_product=np.array(Dot_product)
# Delta_X=np.array(Dot_productX)
# Delta_Y=np.array(Dot_productY)
lambdaxx=np.sqrt((x1-x0)**2+(z1-z0)**2+(y1-y0)**2)
Theta=np.arcsin((lambdaxx/d)*np.cos(np.pi/2-Dot_product))
Thetax=Theta*cphi
Thetay=Theta*sphi
# plt.arrow(x1[0],z1[0],Px[0],Pz[0], width=0.000001, head_length=0.07, head_width=0.0025)

#Definisce i tagli sulle coordinate di arrivo
Deltacut=Dot_product<np.pi/4
tancut=(np.sqrt(X**2+Y**2)/np.max(np.abs(Z-z0)))<np.tan(2.5*np.pi/180)
zcut=Z<d

#Risultato del fit fatto con le condizioni di Hit di Paolo
Hitpaolo=(np.sqrt((Theta*(180/np.pi))**2))<0.15
#scaling con taglio sulle coordinate di arrivo, al momento è disattivato perché non ha nessun senso
#ratio=(np.float(np.sum(tancut*zcut))/len(tancut))
scaling=length/intflusso
# #scaling senza taglio
# scaling=nphot/(intflusso*(BinTeV/Opening))
#MCset=MCset/scaling
cascade,bins=np.histogram(a[4][idx*cascadeidx][Hitpaolo]*1e9,bins=linbin)
cascade=cascade/scaling
#cascade*=1e3
bincenter = 0.5*(bins[1:]+bins[:-1])
#intrinsic,bins=np.histogram(MCset,bins=linbin)
intrinsic=Cutoff(Nzero,XSORT,1.5,2080)
cascade[11:]=0  #Taglia i fotoni fluttuanti a energie troppo alte
plt.figure(figsize=(10,5))
plt.plot(XSORT,intrinsic,label='Intrinsic emission model',linestyle='--')
plt.plot(XSORT,cascade,label='Cascade emission',linestyle='-.')
plt.plot(XSORT,cascade+intrinsic,label='Cascade +intrinsic emission')
plt.ylim(2e-17,3e-10)
# plt.loglog(bincenter,intrinsic/scaling,label='Simulated spectrum')
# plt.loglog(bincenter,intrinsic/scaling+cascade,label='Simulated observed spectrum')
# plt.legend(loc='top right')
plottapuntide()
plt.title('Observed vs simulated flux, 5 degree jet ''Brms='+str(Brms)+'nG')
plt.savefig(savedir+'fit_paolo'+str(Brms)+'.eps', format='eps')
plt.savefig(savedir+'fit_paolo'+str(Brms)+'.png')
print ('saved '+'fit_paolo'+str(Brms)+'.png in' + savedir )
plt.show()
plt.close()

residuals=cascade+intrinsic-YSORT
chi=calcolachi2de(cascade+intrinsic,YSORT,DYSORT,XSORT)
chifermi=calcolachi2de((cascade+intrinsic)[:-13],YSORT[:-13],DYSORT[:-13],XSORT[:-13])
print('Risultati del fit per un campo magnetico di '+str(Brms)+' nG')
#print('Chi quadro completo= '+str(chi))
print ('Chi quadro nella banda 0-300 GeV (paolo)= '+str(chifermi))

def plottamappapaolo(range):
    nbins=int(range/0.001)
    plt.hist2d(Thetax[zcut]*180/np.pi, Thetay[zcut]*180/np.pi, bins=nbins, range=[[-range,range],[-range,range]], norm=LogNorm())
    plt.xlabel('$\Theta$ cos$\phi$')
    plt.ylabel('$\Theta$ sin$\phi$')
    plt.title('Event map')
    plt.colorbar()
    
def plottazytutti():
    plt.scatter(z1, y1,s=0.1)
    
def electrondensity(depth):
    lower=y1>depth-0.5
    upper=y1<depth+0.5
    plt.hist2d(z1[lower*upper],x1[lower*upper],bins=1000, norm=LogNorm())
    plt.xlabel('z1[Mpc]')
    plt.ylabel('x1[Mpc]')
    plt.colorbar()
    plt.show()

def scatter3d():
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(z1, y1, x1, s=1)
    plt.show()


# lambdaxx=np.sqrt((x1-x0)**2+(z1-z0)**2+(y1-y0)**2)
# lambdayy=np.sqrt((y1-y0)**2+(z1-z0)**2+(x1-x0)**2)
# 
# thetapaoloX=np.arcsin((-lambdaxx/np.max(np.abs(Z-z0)))*np.sin(Delta_X))
# thetapaoloY=np.arcsin((-lambdayy/np.max(np.abs(Z-z0)))*np.sin(Delta_Y))

# plt.scatter(z1[:1000],x1[:1000],s=9)
# for j in range (0,1000):
#    plt.arrow(z1[j],x1[j],10*Pz0[j],10*Px0[j], width=0.000001, head_length=10*0.07, head_width=10*0.0025,color='r')
#    plt.arrow(z1[j],x1[j],10*Pz[j],10*Px[j], width=0.000001, head_length=10*0.07, head_width=10*0.0025)
#    plt.arrow(z1[j],x1[j],10*Pz[j],10*Pz[j]*np.tan(thetapaoloX[j]), width=0.000001, head_length=10*0.07, head_width=10*0.0025,color='g')
# plt.xlim(0,-np.max(Z-z0))
# plt.xlabel('Z[Mpc]')
# plt.ylabel('X[Mpc]')

# thetapaolox[np.isnan(thetapaolox)] = 0
# thetapaoloy[np.isnan(thetapaoloy)] = 0
# thetapaolo=[]
# thetapaolo.append(thetapaolox)
# thetapaolo.append(thetapaoloy)
# thetapaolo=np.array(thetapaolo)
# delta=[]
# delta.append(deltax)
# delta.append(deltay)
# delta=np.array(delta)
#Restituisce la colonna richiesta in ordine di arrivo delle particelle, filtrate sulla 0 che dà la distanza percorsa
def ordinarrivo (b, colonna):
    indices=np.argsort(b[0])
    return (np.take_along_axis(b[colonna], indices, axis=0))
  
    
#Achtung: filtrare prima de vuoi plottare i fotoni, devi mandargli a[:,idx] anziché a in input
def plottascatterfreccinexz(b, intervallo):
    phot = np.abs(b[3]) == 22
    zordered=ordinarrivo(b, 25)
    xordered=ordinarrivo(b, 23)
    zorderedp=ordinarrivo(b, 10)
    xorderedp=ordinarrivo(b, 8)
    zorderedhit=ordinarrivo(b, 7)
    xorderedhit=ordinarrivo(b, 5)
    dordered=ordinarrivo(b, 0)
    thetapaoloxordered=ordinarrivo(thetapaolo,0)
    thetapaoloyordered=ordinarrivo(thetapaolo,1)
    
    # plt.scatter(zordered[intervallo[0]:intervallo[1]],xordered[intervallo[0]:intervallo[1]], s=1, label='Point of creation')
    # plt.scatter(zorderedhit[intervallo[0]:intervallo[1]],xorderedhit[intervallo[0]:intervallo[1]], s=16, color='r', marker='x', label='Hit position')
    plt.scatter(zordered[0:intervallo[1]],xordered[0:intervallo[1]], s=1, label='Point of creation')
    plt.scatter(zorderedhit[0:intervallo[1]],xorderedhit[0:intervallo[1]], s=16, color='r', marker='x', label='Hit position')
    plt.rcParams["figure.figsize"] = (17.5,6.5)
    #plt.xlim(min(b[25]),max(b[25]))
    #plt.ylim(min(b[23]),max(b[23]))
    # #My way
    # for j in range (intervallo[0],intervallo[1]):
    #     plt.arrow(x=zordered[j],y=xordered[j],dx=zorderedp[j]*5, dy=xorderedp[j]*5, alpha=0.2)
    #Paolo's way
    for j in range (intervallo[0],intervallo[1]):
        plt.arrow(x=zordered[j],y=xordered[j],dx=-5, dy=-5*np.tan(thetapaoloxordered[j]), alpha=0.2)
    plt.xlim([min(b[25]),max(b[25])])
    plt.ylim([min(b[23]),max(b[23])])
    plt.xlabel('Z [Mpc]')
    plt.ylabel('X [Mpc]')
    plt.title('Dtot[Mpc]='+str(dordered[intervallo[0]])+ ' - ' + str (dordered[intervallo[1]]))
    plt.legend()

def plottascatterfreccineyz(b, intervallo):
    phot = np.abs(b[3]) == 22
    zordered=ordinarrivo(b, 25)
    yordered=ordinarrivo(b, 24)
    zorderedp=ordinarrivo(b, 10)
    yorderedp=ordinarrivo(b, 9)
    zorderedhit=ordinarrivo(b, 7)
    yorderedhit=ordinarrivo(b, 6)
    dordered=ordinarrivo(b, 0)
    # plt.scatter(zordered[intervallo[0]:intervallo[1]],yordered[intervallo[0]:intervallo[1]], s=1, label='Point of creation')
    plt.scatter(zorderedhit[intervallo[0]:intervallo[1]],yorderedhit[intervallo[0]:intervallo[1]], s=16, color='r', marker='x', label='Hit position')
    plt.scatter(zordered[0:intervallo[1]],yordered[0:intervallo[1]], s=1, label='Point of creation')
    # plt.scatter(zorderedhit[0:intervallo[1]],yorderedhit[0:intervallo[1]], s=16, color='r', marker='x', label='Hit position')
    plt.rcParams["figure.figsize"] = (10,10) #era 24x13
    #plt.xlim(min(b[25]),max(b[25]))
    #plt.ylim(min(b[23]),max(b[23]))
    for j in range (intervallo[0],intervallo[1]):
        plt.arrow(x=zordered[j],y=yordered[j],dx=zorderedp[j]*20, dy=yorderedp[j]*20, alpha=0.2)
    plt.xlim([0,max(b[25])])
    plt.ylim([-11,11])
    plt.xlabel('Z [Mpc]')
    plt.ylabel('Y [Mpc]')
    plt.title('Dtot[Mpc]='+str(dordered[intervallo[0]])+ ' - ' + str (dordered[intervallo[1]]))
    plt.legend()

def plottarrivoxysmall(b, intervallo):
    phot = np.abs(b[3]) == 22
    yordered=ordinarrivo(b, 24)
    xordered=ordinarrivo(b, 23)
    yorderedp=ordinarrivo(b, 9)
    xorderedp=ordinarrivo(b, 8)
    yorderedhit=ordinarrivo(b, 6)
    xorderedhit=ordinarrivo(b, 5)
    dordered=ordinarrivo(b, 0)
    # plt.scatter(yordered[intervallo[0]:intervallo[1]],xordered[intervallo[0]:intervallo[1]], s=1, label='Point of creation')
    #plt.hist2d(yorderedhit[0:intervallo[1]],xorderedhit[0:intervallo[1]], s=16, color='r', marker='x', bins=100, range=[[min(b[24]),max(b[24])],[min(b[23]),max(b[23])]])
    plt.scatter(yorderedhit[0:intervallo[1]],xorderedhit[0:intervallo[1]], s=1, color='b', marker='x', label='All hit position')
    plt.scatter(yorderedhit[intervallo[0]:intervallo[1]],xorderedhit[intervallo[0]:intervallo[1]], s=16, color='r', marker='x', label='Current hit position')
    # plt.scatter(zordered[0:intervallo[1]],xordered[0:intervallo[1]], s=1, label='Point of creation')
    # plt.scatter(zorderedhit[0:intervallo[1]],xorderedhit[0:intervallo[1]], s=16, color='r', marker='x', label='Hit position')
    plt.rcParams["figure.figsize"] = (6.5,6.5)
    #plt.xlim(min(b[25]),max(b[25]))
    #plt.ylim(min(b[23]),max(b[23]))
    # for j in range (intervallo[0],intervallo[1]):
    #     plt.arrow(x=zordered[j],y=xordered[j],dx=zorderedp[j]*5, dy=xorderedp[j]*5, alpha=0.2)
    plt.xlim([min(b[23]),max(b[23])])
    plt.ylim([min(b[24]),max(b[24])])
    plt.xlabel('Y [Mpc]')
    plt.ylabel('X [Mpc]')
    plt.title('Dtot[Mpc]='+str(dordered[intervallo[0]])+ ' - ' + str (dordered[intervallo[1]]))
    plt.legend()
    
def plottarrivoxy(b, intervallo):
    phot = np.abs(b[3]) == 22
    yordered=ordinarrivo(b, 24)
    xordered=ordinarrivo(b, 23)
    yorderedp=ordinarrivo(b, 9)
    xorderedp=ordinarrivo(b, 8)
    yorderedhit=ordinarrivo(b, 6)
    xorderedhit=ordinarrivo(b, 5)
    dordered=ordinarrivo(b, 0)
    # plt.scatter(yordered[intervallo[0]:intervallo[1]],xordered[intervallo[0]:intervallo[1]], s=1, label='Point of creation')
    #plt.hist2d(yorderedhit[0:intervallo[1]],xorderedhit[0:intervallo[1]], s=16, color='r', marker='x', bins=100, range=[[min(b[24]),max(b[24])],[min(b[23]),max(b[23])]])
    plt.scatter(yorderedhit[0:intervallo[1]],xorderedhit[0:intervallo[1]], s=1, color='b', marker='x', label='All hit position')
    plt.scatter(yorderedhit[intervallo[0]:intervallo[1]],xorderedhit[intervallo[0]:intervallo[1]], s=16, color='r', marker='x', label='Current hit position')
    # plt.scatter(zordered[0:intervallo[1]],xordered[0:intervallo[1]], s=1, label='Point of creation')
    # plt.scatter(zorderedhit[0:intervallo[1]],xorderedhit[0:intervallo[1]], s=16, color='r', marker='x', label='Hit position')
    plt.rcParams["figure.figsize"] = (6.5,6.5)
    #plt.xlim(min(b[25]),max(b[25]))
    #plt.ylim(min(b[23]),max(b[23]))
    # for j in range (intervallo[0],intervallo[1]):
    #     plt.arrow(x=zordered[j],y=xordered[j],dx=zorderedp[j]*5, dy=xorderedp[j]*5, alpha=0.2)
    plt.xlim([min(b[5]),max(b[5])])
    plt.ylim([min(b[6]),max(b[6])])
    plt.xlabel('Y [Mpc]')
    plt.ylabel('X [Mpc]')
    plt.title('Dtot[Mpc]='+str(dordered[intervallo[0]])+ ' - ' + str (dordered[intervallo[1]]))
    plt.legend()

def plottanimazione (fotogrammi, b):
    phot = np.abs(b[3]) == 22
    c=b[:,idx]
    size=np.int(math.floor(len(c[23])/fotogrammi))
    for j in range (0,fotogrammi):
        print ('Processing photons '+str(j*size)+' to '+ str((j+1)*size))
        plottascatterfreccinexz(c, (j*size, (j+1)*size))
        # plottascatterfreccineyz(c, (j*size, (j+1)*size))
        # plt.savefig(savedir+'anim/' + str(j+100)+ '.png')
        plt.savefig(savedir+'anim/xz/' + str(j+100)+ '.png')
        plt.close()
#Esempio: plotta l'immagine dell'anellino>2° su 50 mpc di coordinate di arrivo:  kdeplot(X[np.degrees(np.sqrt(Px0**2+Py0**2))>2], Y[np.degrees(np.sqrt(Px0**2+Py0**2))>2], 50)
def kdeplot (x_full,y_full, radius):
    x=x_full[(np.sqrt(x_full**2+y_full**2)<radius*1.5)]
    y=y_full[(np.sqrt(x_full**2+y_full**2)<radius*1.5)]
    from scipy.stats import gaussian_kde
    data = np.vstack([x, y])
    kde = gaussian_kde(data)
    xgrid = np.linspace(-radius, radius, radius*5)
    ygrid = np.linspace(-radius, radius, radius*5)
    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))
    plt.imshow(Z.reshape(Xgrid.shape),
           origin='lower', aspect='auto',
           extent=[-radius, radius, -radius, radius],
           cmap='hot')
    plt.xlabel('X[Mpc]')
    plt.ylabel('Y[Mpc]')
    #cb = plt.colorbar()
    #cb.set_label("density")
    plt.show()

def plottamappa(x_full,y_full,px_full,py_full,radius,centerx,centery):
    x=x_full[(np.sqrt((x_full-centerx)**2+(y_full-centery)**2)<radius*1.5)]
    px=px_full[(np.sqrt((x_full-centerx)**2+(y_full-centery)**2)<radius*1.5)]
    y=y_full[(np.sqrt((x_full-centerx)**2+(y_full-centery)**2)<radius*1.5)]
    py=py_full[(np.sqrt((x_full-centerx)**2+(y_full-centery)**2)<radius*1.5)]
    from scipy.stats import gaussian_kde
    data = np.vstack([np.rad2deg(np.arcsin(px)), np.rad2deg(np.arcsin(py))])
    kde = gaussian_kde(data)
    xgrid = np.linspace(-5, 5, 300)
    ygrid = np.linspace(-5, 5, 300)
    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))
    plt.imshow(Z.reshape(Xgrid.shape),
           origin='lower', aspect='auto',
           extent=[-radius, radius, -radius, radius],
           cmap='hot')
    plt.xlabel('$\Theta_x$[degrees]')
    plt.ylabel('$\Theta_y$[degrees]')
    cb = plt.colorbar()
    #cb.set_label("density")
    plt.show()        

# binspace=np.logspace(np.log10(10),np.log10(2000),num=20)
# plt.hist(a[4][idx*cascadeidx]*1e9,bins=binspace)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()
# 
# Salva il file completo a meno che non esista già
if os.path.exists(savedir+'final'+str(Brms)+'.txt'):
   print ('Global data file already exists, doing nothing')
else:
   np.savetxt(savedir+'final'+str(Brms)+'.txt',a)



