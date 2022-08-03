from crpropa import *
import os
import numpy as np
import time
#Creo il campo magnetico: crea una grid size e definisci il campo magnetico turbolento con spettro (Kolmogorov -11/3) e rms in questo caso di 8nG
#Plotta E in funzione di E0, deve venire quadratico
#Partendo da neronov ceca di capire quali possono essere i plot diagnostici, Neronov2009
savedir='/home/kibuzo/TesiMagistrale/Python/Outfigs/'
z=0.21
D=redshift2ComovingDistance(z)
#randomSeed = 18
def SimulaB(Brms, nphot,i):
    gridpoints=256 #512 Suggerito nel forum di crpropa, meglio non scendere sotto
    gridsize=400*Mpc
    gridspacing=gridsize/gridpoints
    minStep = gridspacing/100#20 * kpc
    maxStep = gridspacing/10#400 * kpc
    minco=2*gridspacing
    maxco=gridsize/2
    vgrid = VectorGrid(Vector3d(0,0,D), gridpoints, gridspacing)
    Bfield = MagneticFieldGrid(vgrid)
    initTurbulence(vgrid, Brms*nG, minco, maxco, -11.7, randomSeed)
    print('Turbulent correlation length: '+str(np.round(turbulentCorrelationLength(minco,maxco)/Mpc,decimals=2))+' Mpc\n')
    
    tol = 1e-4
    sz = D
    dz = 0.005
    import time
    
    
    CutoffTeV=2.08
    obs = Observer()
    obs.add(ObserverLargeSphere(Vector3d(0,0,D),D))
    #obs.add(ObserverSmallSphere(Vector3d(0,0,0),30*Mpc))
    filename='electron_output'+str(i)+'.txt'
    if os.path.exists(filename):
        os.remove(filename)
    t = TextOutput(filename,Output.Everything)
    obs.onDetection(t)
    source = Source()
    source.add(SourcePosition(Vector3d(0,0,D))) #Vettore 3d della distanza dell'oggetto
    source.add(SourceRedshift(z))
    source.add(SourceParticleType(11)) #22 significa fotone, e' uno standard definito da pdg credo
    #source.add(SourcePowerLawSpectrum(1 * GeV, CutoffTeV * TeV, -1.5) ) #Emin, Emax, Index della power law
    source.add(SourceEnergy(5*TeV))
    source.add(SourceEmissionCone(Vector3d(0,0,-1),(5/360.)*np.pi))
    
    
    sim = ModuleList()
    sim.add(FutureRedshift())
    bound=SphericalBoundary(Vector3d(0,0,D),D+maxStep)
    sim.add(bound)
    sim.add(PropagationBP(Bfield, tol, minStep, maxStep))
    sim.add(obs)
    tic=time.time()
    sim.run(source,nphot,True) #quel numero sono il numero di gamma sparati
    toc=time.time()
    print (toc-tic)
   # filename='photon_electron_output'+str(Brms)+'.txt'
   # a=loadtxt(filename)
   # print np.arctan(np.max(np.sqrt(a[:,17]**2+a[:,18]**2)/np.sqrt(a[:,19]**2+a[:,17]**2+a[:,18]**2)))*180/np.pi

for i in range (3,9):
    randomSeed = int(time.clock())
    print('\n'+'Running simulation '+str(i-2)+ ' of ' +str(9-3))
    #Brms=float('1e-'+str(i)) 
    Brms=1e-8
    #Brms+=i*Brms/100
    nphot=int(1e4)
    SimulaB(Brms, nphot,i)