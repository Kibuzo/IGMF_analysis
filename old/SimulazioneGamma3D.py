from crpropa import *
#Creo il campo magnetico: crea una grid size e definisci il campo magnetico turbolento con spettro (Kolmogorov -11/3) e rms in questo caso di 8nG
#Plotta E in funzione di E0, deve venire quadratico
#Partendo da neronov ceca di capire quali possono essere i plot diagnostici, Neronov2009
z=0.21
D=redshift2ComovingDistance(z)
randomSeed = 40
minStep = 10. * kpc
maxStep = 1000 * kpc
gridpoints=250
gridsize=5000*kpc
gridspacing=gridsize/gridpoints
vgrid = VectorGrid(Vector3d(D/2), gridpoints, gridspacing)
Bfield = MagneticFieldGrid(vgrid)
initTurbulence(vgrid, 0*nG, 100*kpc , int(1e3)*kpc, -11./3., randomSeed)

# minStep=D/1000
# maxStep=D/100
tol = 1e-4
sz = D
dz = 0.005
nphot=100 #conta il numero di fotoni ammessi dai tuoi dati.
import time

# Tolti i seguenti moduli:
# pp su URB (che dovrebbe essere Universal Radio Background), nell'ipotesi che abbia energia troppo bassa
# Triple e Double pair production, che dovrebbe essere rilevante solo per energie sopra il PeV

CutoffTeV=2.08
obs = Observer()
D0=25.5*Mpc
#obs.add(ObserverSurface(Sphere(Vector3d(D,D,D),D)))
#obs.add(ObserverSmallSphere(Vector3d(0,D,D),D0))  #6400*km))
obs.add(ObserverLargeSphere(Vector3d(0,0,D),D))
#obs.add(ObserverInactiveVeto()) #non so cosa sia
t = TextOutput("photon_electron_output.txt",Output.Everything)
#obs.add(ObserverRedshiftWindow(-dz,dz))
obs.onDetection(t)
source = Source()
source.add(SourcePosition(Vector3d(0,0,D))) #Vettore 3d della distanza dell'oggetto
source.add(SourceRedshift(z))
source.add(SourceParticleType(22)) #22 significa fotone, e' uno standard definito da pdg credo
source.add(SourcePowerLawSpectrum(1 * GeV, CutoffTeV * TeV, -1.5) ) #Emin, Emax, Index della power law
source.add(SourceEmissionCone(Vector3d(0,0,-1),1/18.)) #Achtung! la direzione è definita come posizione0-posizione1
#source.add(SourceUniformSphere(Vector3d(D,D,D),5*Mpc))

sim = ModuleList()
#sim.add(MinimumEnergy(1e9*eV))
sim.add(Redshift())
sim.add(EMPairProduction(IRB_Franceschini08,True))
#bound=SphericalBoundary(Vector3d(D,D,D),D+minStep)
bound=SphericalBoundary(Vector3d(0,0,D),D+maxStep)
sim.add(bound)#Spherical Boundary serve a flaggare come inattiva una particella che fa overshooting e supera l'osservatore. L'ho fatta delle dimensioni minime dello step di propagazione altrimenti ne muore troppe.
sim.add(PropagationCK(Bfield, tol, minStep, maxStep))
sim.add(EMInverseComptonScattering(CMB,True))
sim.add(obs)
tic=time.time()
sim.run(source,nphot,True) #quel numero sono il numero di gamma sparati
toc=time.time()
print (toc-tic)
#Sezione per creare il plot del flusso


#------------------------Parte 2: entra il campo magnetico------------------------------#
initTurbulence(vgrid, 1e-7*nG, 2*gridspacing, 200*gridspacing, -11./3., randomSeed)
print('Turbulent correlation length: '+str(np.round(turbulentCorrelationLength(2*gridspacing,20*gridspacing)/Mpc,decimals=2))+' Mpc')
obs2 = Observer()
#obs.add(ObserverSurface(Sphere(Vector3d(D,D,D),D)))
obs2.add(ObserverSmallSphere(Vector3d(0,0,D),D0))#6400*km))
#obs.add(ObserverInactiveVeto()) #non so cosa sia
t = TextOutput("photon_electron_output_BFIELD.txt",Output.Everything)
obs2.add(ObserverRedshiftWindow(-dz,dz))
obs2.onDetection(t)
source2 = Source()
source2.add(SourcePosition(Vector3d(D,D,D))) #Vettore 3d della distanza dell'oggetto
source2.add(SourceRedshift(z))
source2.add(SourceParticleType(22)) #22 significa fotone, e' uno standard definito da pdg credo
source2.add(SourcePowerLawSpectrum(100 * GeV, CutoffTeV * TeV, -1.5) ) #Emin, Emax, Index della power law
#source.add(SourceEnergy(10*TeV))
source2.add(SourceEmissionCone(Vector3d(0,0,1),5.))


simb = ModuleList()
simb.add(Redshift())
simb.add(EMPairProduction(IRB_Franceschini08,True))
#bound=SphericalBoundary(Vector3d(D,D,D),D+minStep)
bound=SphericalBoundary(Vector3d(0,D,D),D0+maxStep)
simb.add(bound)#Spherical Boundary serve a flaggare come inattiva una particella che fa overshooting e supera l'osservatore. L'ho fatta delle dimensioni minime dello step di propagazione altrimenti ne muore troppe.
simb.add(PropagationCK(Bfield, tol, minStep, maxStep))
simb.add(EMInverseComptonScattering(CMB,True))
simb.add(obs2)
tic=time.time()
simb.run(source2,nphot,True) #quel numero sono il numero di gamma sparati
toc=time.time()
print (toc-tic)
#Sezione per creare il plot del flusso


from pylab import *

t.close()
figure(figsize=(8,8))
nbins=30

a = loadtxt("photon_electron_output.txt")
b = loadtxt("photon_electron_output_BFIELD.txt")
E = logspace(9,np.log10(CutoffTeV*10**12),nbins)
# ------------Output3D Plot------------------#
# idx = a[:,1] == 22 #Seleziona gli indici che contengono fotoni
# idxb = b[:,1] == 22
# photons = a[idx,2] * 1e18
# photonsrepridx=np.logical_and(~np.isnan(a[:,6]),idx)
# photonsrepr=a[photonsrepridx,2]*1e18
# photonsb= b[idxb,2] * 1e18
# photonsbrepridx=np.logical_and(~np.isnan(b[:,6]),idxb)

# --------------------Everything plot------------------#
idx = a[:,3] == 22 #Seleziona gli indici che contengono fotoni
idxb = b[:,3] == 22
photons = a[idx,4] * 1e18
photonsrepridx=np.abs(a[:,21])==11#np.logical_and(np.abs(a[:,21])==11,idx)
photonsrepr=a[photonsrepridx,2]*1e18
photonsb= b[idxb,4] * 1e18
photonsbrepridx=np.abs(b[:,21])==11#np.logical_and(np.abs(b[:,21])==11,idxb)
photonsbrepr=b[photonsbrepridx,2]*1e18
#idx = fabs(a[:,1]) == 11
#ep = a[idx,2] * 1e18
data,bins = histogram(photons,E)
dataB,binsB = histogram(photonsb,E)
datarepr,binsrepr = np.histogram(photonsrepr,E)
databrepr,binsbrepr = np.histogram(photonsbrepr,E)
bincenter = (E[1:] -E[:-1])/2 + E[:-1]
plot(bincenter, data,label="no magnetic field")
plot(bincenter, dataB,label="magnetic field")
plot(bincenter, datarepr,label="no magnetic field, reprocessed")
plot(bincenter, databrepr,label="magnetic field, reprocessed")
#data,bins = histogram(ep,E)
#plot(bincenter, data, label="electrons / positrons")
grid()
loglog()
#xlim(1e16, 1e21)
#ylim(1e2, 1e4)
legend(loc="upper right")
xlabel("Energy [eV]")
ylabel("Number of Photons/bin")
title("flux of particles in 35 logarithmically spaced bins")
print len (a[:,1]) / nphot
show()

def countinside(radius, a):
    return np.sum(np.sqrt(a[:,6]**2+a[:,5]**2)<radius)

vecs=[]
vettore=np.linspace(0,100,num=100)

for j in vettore:
    vecs.append(countinside(j,a))
    
    
contained=np.sqrt(a[:,6]**2+a[:,5]**2)<10 #Contenuti nel cerchio di 10 megaparsec

