from crpropa import *
#Creo il campo magnetico: crea una grid size e definisci il campo magnetico turbolento con spettro (Kolmogorov -11/3) e rms in questo caso di 8nG
z=0.21
D=redshift2ComovingDistance(z)
randomSeed = 400
# gridpoints=100
# gridsize=50*Mpc
# gridspacing=gridsize/gridpoints
# vgrid = VectorGrid(Vector3d(D), gridpoints, gridspacing)
# initTurbulence(vgrid, 8*nG, 2*gridspacing, 20*gridspacing, -11./3., randomSeed)
# Bfield = MagneticFieldGrid(vgrid)
# minStep = 80. * kpc
# maxStep = 1000 * kpc
# tol = 1e-4
# sz = D
dz = 0.0005
    
sim = ModuleList()
sim.add(SimplePropagation())
sim.add(Redshift())
sim.add(EMPairProduction(IRB_Franceschini08,True))
sim.add(EMInverseComptonScattering(IRB_Franceschini08,True))
sim.add(EMInverseComptonScattering(CMB,True))
sim.add(EMInverseComptonScattering(URB_Protheroe96,True))


CutoffTeV=2.08
obs = Observer()
#obs.add(ObserverPoint()) # Passato a largeSphere
obs.add(ObserverRedshiftWindow(-dz,dz))
obs.add(ObserverSurface(Sphere(Vector3d(D,D,D),D)))
#obs.add(ObserverInactiveVeto()) #non so cosa sia
t = TextOutput("photon_electron_output.txt",Output.Event3D)
obs.onDetection(t)
source = Source()
source.add(SourcePosition(Vector3d(D,D,D))) #Vettore 3d della distanza dell'oggetto
source.add(SourceEmissionCone(Vector3d(0,0,0),.5))
source.add(SourceRedshift1D())
source.add(SourceParticleType(22)) #22 significa fotone, e' uno standard definito da pdg credo
source.add(SourcePowerLawSpectrum(1 * GeV, CutoffTeV * TeV, -1.5) ) #Emin, Emax, Index della power law
sim.add(obs)
sim.run(source,100000,True) #quel numero sono il numero di gamma sparati

#Sezione per creare il plot del flusso
from pylab import *

t.close()
figure(figsize=(6,6))

a = loadtxt("photon_electron_output.txt")
E = logspace(9,np.log10(CutoffTeV*10**12),71)
idx = a[:,1] == 22
photons = a[idx,2] * 1e18
idx = fabs(a[:,1]) == 11
ep = a[idx,2] * 1e18
data,bins = histogram(photons,E)
bincenter = (E[1:] -E[:-1])/2 + E[:-1]
plot(bincenter, data,label="photons")
data,bins = histogram(ep,E)
plot(bincenter, data, label="electrons / positrons")
grid()
loglog()
#xlim(1e16, 1e21)
#ylim(1e2, 1e4)
legend(loc="upper right")
xlabel("Energy [eV]")
ylabel("Number of Particles")
show()