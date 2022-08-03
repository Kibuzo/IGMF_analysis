from crpropa import *


#Definisci come propagare la particella, che moduli della simulazione attivare e poi crei particella e osservatore



# simulation: a sequence of simulation modules
sim = ModuleList()

# add propagator for rectalinear propagation
sim.add( SimplePropagation() )

# add interaction modules
sim.add( PhotoPionProduction(CMB) )
sim.add( ElectronPairProduction(CMB) )
sim.add( NuclearDecay() )
sim.add( MinimumEnergy( 1 * EeV) )

# define particle to propagate
cosmicray = Candidate(nucleusId(1,1), 200 * EeV, Vector3d(100 * Mpc, 0, 0))

sim.run(cosmicray)
print cosmicray
print 'Propagated distance', cosmicray.getTrajectoryLength() / Mpc, 'Mpc'

# add an observer: la particella parte dalla coordinata positiva definita prima e si muove fino a quella messa da te
obs = Observer()
obs.add( ObserverPoint() )  # observer at x = 0
sim.add(obs)
print obs

# adesso se vuoi puoi salvare tutte le informazioni su un file di output
# trajectory output
output1 = TextOutput('trajectories.txt', Output.Trajectory1D)
#sim.add(output1)  # generates a lot of output

#output1.disable(Output.RedshiftColumn)  # don't save the current redshift
#output1.disableAll()  # disable everything to start from scratch
#output1.enable(Output.CurrentEnergyColumn)  # current energy
#output1.enable(Output.CurrentIdColumn)      # current particle type

#Volendo puoi aggiungere un altro file di output
# event output
output2 = TextOutput('events.txt', Output.Event1D)
obs.onDetection(output2)

#sim.run(cosmicray)
#output2.close()

#Invece del singolo fotone puoi definire una sorgente di cosmici
# cosmic ray source
source = Source()
source.add( SourcePosition(100 * Mpc) )
source.add( SourceParticleType(nucleusId(1, 1)) )
source.add( SourcePowerLawSpectrum(1 * EeV, 200 * EeV, -1) )
print source

#Adessi fai girare tutto 
sim.setShowProgress(True)  # switch on the progress bar
sim.run(source, 10000)

#plotti
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

output2.close()  # close output file before loading
data = np.genfromtxt('events.txt', names=True)
print 'Number of events', len(data)

logE0 = np.log10(data['E0']) + 18
logE  = np.log10(data['E']) + 18

plt.figure(figsize=(10, 7))
h1 = plt.hist(logE0, bins=25, range=(18, 20.5), histtype='stepfilled', alpha=0.5, label='At source')
h2 = plt.hist(logE,  bins=25, range=(18, 20.5), histtype='stepfilled', alpha=0.5, label='Observed')
plt.xlabel('log(E/eV)')
plt.ylabel('N(E)')
plt.legend(loc = 'upper left', fontsize=20)