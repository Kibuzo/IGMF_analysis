from crpropa import *
import os
import numpy as np
import time
#Plotta E in funzione di E0, deve venire quadratico
#Partendo da neronov ceca di capire quali possono essere i plot diagnostici, Neronov2009
#savedir='/home/kibuzo/TesiMagistrale/Python/Outfigs/' 
savedir='./Data/'
if not os.path.exists(savedir):
    os.mkdir(savedir)
z=0.21
D=redshift2ComovingDistance(z)
#randomSeed = 40
nphot=10
def SimulaB(Brms, nphot,i):
    randomSeed = 1987#int(time.clock())     #Sceglie un seed pseudocasuale
    Brms=1e-7      #commentare se si vogliono tutti Brms diversi forniti in input
    gridpoints=256 #512 Suggerito nel forum di crpropa, meglio non scendere sotto, abbi 8GB di ram liberi
    gridsize=200*Mpc
    gridspacing=gridsize/gridpoints
    minStep = 1e-4*kpc#1*gridspacing
    maxStep = 500*kpc#10*gridspacing
    minco=2*gridspacing
    maxco=gridsize/2
    vgrid = VectorGrid(Vector3d(0,0,D), gridpoints, gridspacing)
    Bfield = PeriodicMagneticField(MagneticFieldGrid(vgrid), Vector3d(1.2 * D))
    initTurbulence(vgrid, Brms*nG, minco, maxco, -11./3, randomSeed)
    print('Turbulent correlation length: '+str(np.round(turbulentCorrelationLength(minco,maxco)/Mpc,decimals=2))+' Mpc\n')
    
    tol = 1e-2
    sz = D
    dz = 0.005
    
    # Tolti i seguenti moduli:
    # pp su URB (che dovrebbe essere Universal Radio Background), nell'ipotesi che abbia energia troppo bassa
    # Triple e Double pair production, che dovrebbe essere rilevante solo per energie sopra il PeV
    
    CutoffTeV=2.08
    obs = Observer()
    #obs.add(ObserverLargeSphere(Vector3d(0,0,D),D))
    obs.add(ObserverSmallSphere(Vector3d(0,0,0),40*Mpc))
    #filename='photon_electron_output'+str(Brms)'.txt'  #usare per BRMS diversi
    filename=savedir+'photon_electron_output'+str(Brms)+str(i)+'.txt'  #Usare per BRMS uguali
    if os.path.exists(filename):
        os.remove(filename)
    t = TextOutput(filename,Output.Everything)
    obs.onDetection(t)
    source = Source()
    source.add(SourcePosition(Vector3d(0,0,D))) #Vettore 3d della distanza dell'oggetto
    source.add(SourceRedshift(z))
    source.add(SourceParticleType(22)) #22 significa fotone, e' uno standard definito da pdg credo
    #source.add(SourcePowerLawSpectrum(100 * GeV, CutoffTeV * TeV, -1.5) ) #Emin, Emax, Index della power law
    source.add(SourceEnergy(15*TeV))
    source.add(SourceEmissionCone(Vector3d(0,0,-1),(5/360.)*np.pi))
    
    
    sim = ModuleList()
    sim.add(MinimumEnergy(1e9 * eV))
    sim.add(FutureRedshift())
    sim.add(EMPairProduction(IRB_Franceschini08,True))
    bound=SphericalBoundary(Vector3d(0,0,D),D+maxStep)
    sim.add(bound)
    #sim.add(PropagationBP(Bfield, tol, minStep, maxStep))
    sim.add(PropagationBP(Bfield, tol,minStep,maxStep))
    sim.add(EMInverseComptonScattering(CMB,True))
    sim.add(obs)
    tic=time.time()
    sim.run(source,nphot,True) #quel numero sono il numero di gamma sparati
    toc=time.time()
    print (toc-tic)


for i in range (3,9):
    print('\n'+'Running simulation '+str(i-2)+ ' of ' +str(9-3))
    Brms=float('1e-'+str(i)) 
    SimulaB(Brms, nphot,i)
