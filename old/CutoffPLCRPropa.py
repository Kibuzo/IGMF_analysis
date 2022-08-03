##Achtung: ricordati prima di generare i file dello spettro con MontecarloCutoff.py##
from crpropa import *
import os
import numpy as np
import time
#Plotta E in funzione di E0, deve venire quadratico
#Partendo da neronov ceca di capire quali possono essere i plot diagnostici, Neronov2009
# --------------   Create Folders -----------------# 
savedir='./BPL/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
savedirtmp='./BPL/TEMP/'
if not os.path.exists(savedirtmp):
    os.mkdir(savedirtmp)
samplefile='./MCsample/mc0.txt'
spettro=np.loadtxt(samplefile)
nphot=len(spettro)

# --------------------   Simulation parameters    ----------------------- #
randomSeed = 1987#int(time.clock())     #Sceglie un seed pseudocasuale
#randomSeed=320
z=0.21
D=redshift2ComovingDistance(z)
#randomSeed = 40
#nphot=100
Brms=1e-9
gridpoints=256 #512 Suggerito nel forum di crpropa, meglio non scendere sotto, abbi 8GB di ram liberi
gridsize=1*Mpc
gridspacing=gridsize/gridpoints
minStep = 1e-4*kpc#2*gridspacing
maxStep = 500*kpc#20*gridspacing
minco=2*gridspacing
maxco=gridsize/2
vgrid = VectorGrid(Vector3d(0,0,D), gridpoints, gridspacing)
#vgrid = VectorGrid(Vector3d(0,0,0), gridpoints, gridspacing)
Bfield = MagneticFieldGrid(vgrid)
initTurbulence(vgrid, Brms*nG, minco, maxco, -11./3, randomSeed)
#initTurbulence(vgrid, *nG, minco, maxco, -11./3, randomSeed)
Bfield2 = PeriodicMagneticField(MagneticFieldGrid(vgrid), Vector3d(1.2 * D))
print('Turbulent correlation length: '+str(np.round(turbulentCorrelationLength(minco,maxco)/Mpc,decimals=2))+' Mpc\n')
tol = 1e-4
sz = D
dz = 0.005
CutoffTeV=2.08

# ------  Monte Carlo section  ------------#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
Cut=CutoffTeV*1000                                      #In GeV
Nzero=3.83e-14
nsample=nphot                                             #Montecarlo samples
extent=[100.,1e4]                                       #estremi in energia della PDF
gridpointsMC=int(extent[1]-extent[0])*10                  #numero di punti della griglia per la discretizzazione della CDF
estep=(extent[1]-extent[0])/gridpointsMC                  #passo della griglia
egrid=np.linspace(extent[0],extent[1],num=gridpointsMC)   #definizione della griglia

# ------  Monte Carlo functions -----------#


def GridToE(gridpoint):
    return (extent[0]+gridpoint*estep)
    
def Cutoff(E,N0,gamma,Ecut):
    return ((N0*(E/300.)**-gamma)*np.exp(-E/Ecut))

#Segmenta con rettangoli la funzione e  restituisce la CDF grigliata grid
def rettangoli(N0,Gamma,Ecut):
    CDF=[0]
    PDF=[]
    for x in egrid:
        areasup=estep*Cutoff(x,Nzero,Gamma,Ecut)
        areainf=estep*Cutoff(x+estep,Nzero,Gamma,Ecut)
        CDF.append((areasup+areainf)/2+CDF[-1])
        PDF.append((areasup+areainf)/2)
    return(CDF,PDF)

#Gli do uno xi generato a caso tra emin ed emax e lui mi ricava l'energia corrispondente nella CDF
def eval(xi,CDF):
    j=np.min(np.where(CDF>=xi)) #Funziona ed è veloce, che tu ci creda o no
    return (GridToE(j))

    
# ------ Modulo da chiamare per eseguire la simulazione n-esima di un fotone a energia e -------#
def SimulaB(n, e):
    tol = 1e-4
    sz = D
    dz = 0.005
    CutoffTeV=2.08
    obs = Observer()
    #obs.add(ObserverLargeSphere(Vector3d(0,0,D),D))
    obs.add(ObserverSmallSphere(Vector3d(0,0,0),30*Mpc))
    filename=savedirtmp+'photon_electron_output'+str(Brms)+str(n)+'.txt'
    if os.path.exists(filename):
        os.remove(filename)
    t = TextOutput(filename,Output.Everything)
    obs.onDetection(t)
    source = Source()
    source.add(SourcePosition(Vector3d(0,0,D))) #Vettore 3d della distanza dell'oggetto
    source.add(SourceRedshift(z))
    source.add(SourceParticleType(22)) #22 significa fotone, e' uno standard definito da pdg credo
    #source.add(SourcePowerLawSpectrum(100 * GeV, CutoffTeV * TeV, -1.5) ) #Emin, Emax, Index della power law
    source.add(SourceEnergy(float(e)*GeV))
    source.add(SourceEmissionCone(Vector3d(0,0,-1),(5/360.)*np.pi))
    sim = ModuleList()
    sim.add(MinimumEnergy(1e9 * eV))
    sim.add(FutureRedshift())
    sim.add(EMPairProduction(IRB_Franceschini08,True))
    bound=SphericalBoundary(Vector3d(0,0,D),D+maxStep)
    sim.add(bound)
    sim.add(PropagationBP(Bfield, tol, minStep, maxStep))
    sim.add(EMInverseComptonScattering(CMB,True))
    sim.add(obs)
    tic=time.time()
    sim.run(source,1,True) #quel numero sono il numero di gamma sparati
    toc=time.time()
    print (toc-tic)

# # ---------------MC starts ------------------ #
# print 'Generating Monte Carlo Sample of'+str(nphot)+' photons'
# tic=time.time()
# CDF,PDF=rettangoli(Nzero,1.5,Cut) #CDF, PDF
# CDF=CDF/max(CDF)                  #Normalized CDF
# spettro=[]
# for j in range (0,int(nsample)):
#     xi=np.random.random()        
#     spettro.append(eval(xi,CDF))
# toc=time.time()
# print 'Monte Carlo Sample generated in '+str(np.round(toc-tic,decimals=1))+' seconds'
# 
# # --------------- MC diagnostic plot ------------------#
# plt.close()
# nbin=nphot/5
# data,bins=np.histogram(spettro,bins=nbin) #I bin adesso sono in energia perché eval returna energie
# bincenters = 0.5*(bins[1:]+bins[:-1])
# e=np.linspace(min(bincenters),max(bincenters),num=nbin)
# func=Cutoff(e,Nzero,1.5,Cut)
# 
# #Confronto plot di PDF valutata al centro del bin e istogramma montecarlo
# plt.loglog(e,func,label='Ideal BPL Flux')
# plt.loglog(bincenters,data/(np.sum(data)/np.sum(func)),label='Monte Carlo flux, rescaled')
# 
# plt.legend(loc='lower left')
# plt.savefig(savedir+'MC_SAMPLE.png')
# plt.close()
# 
# i=0
tic=time.time()
print 'Starting crpropa simulation\n'

for i in range (0,nphot):
    print '\nPropagating photon '+str(i+1)+' of '+str(nphot)+' with energy '+str(spettro[i])+' GeV'
    SimulaB(i,spettro[i])

toc=time.time()
print 'Simulation ended after'+str(toc-tic)+' seconds'

# --- Diagnostic plot: turbulent magnetic field appearence ------- #
sizeplot=int(10*maxco/Mpc)#int(D/(100))
Bf=np.zeros((sizeplot*1000,sizeplot*1000))
for x in range (0,sizeplot*1000):
    for y in range (0,sizeplot*1000):
        xfield=Bfield.getField(Vector3d(x*Mpc/1000,y*Mpc/1000,0))[0]/nG
        yfield=Bfield.getField(Vector3d(x*Mpc/1000,y*Mpc/1000,0))[1]/nG
        zfield=Bfield.getField(Vector3d(x*Mpc/1000,y*Mpc/1000,0))[2]/nG
        Bf[x,y]=np.sqrt(xfield**2+yfield**2+zfield**2)

plt.imshow(Bf)#,norm=mpl.colors.LogNorm())
# for x in range (0,sizeplot/5):
#     for y in range (0,sizeplot/5):
#         xfield=Bfield.getField(Vector3d(x*5*Mpc,y*5*Mpc,0))[0]
#         yfield=Bfield.getField(Vector3d(x*5*Mpc,y*5*Mpc,0))[1]
#         plt.arrow(x*10,y*10,2*xfield/np.sqrt(xfield**2+yfield**2),2*yfield/np.sqrt(xfield**2+yfield**2),color='red',head_width=1)
plt.title('Simulated magnetic field slice (z=0)')
plt.xlabel('x coordinate [kpc]')
plt.ylabel('y coordinate [kpc]')
plt.colorbar(label="Magnetic field strength [nG]")
#plt.savefig(savedir+'Turbulent magnetic field.png')
plt.show()
#plt.close()
        
    
        
        
        
    
