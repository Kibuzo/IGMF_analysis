from IO import _crpropa
import matplotlib.pyplot as plt
import numpy

toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'

my_sim = _crpropa.simulation(1000, toy_sim_dir, '2data_final_largecone', '1e-08')
print ((my_sim.data['E']))
my_sim.dermer_cut(90)
plt.hist(my_sim.data['P0x'], bins = 100)
plt.xscale('log')
plt.yscale('log')
plt.show()
