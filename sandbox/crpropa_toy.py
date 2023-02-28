from IO import _crpropa
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
from math_.base import radians_to_degree
import logging
toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'
#toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/1e-6/BPL'


my_sim = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-08', overwrite = True)
non_interacting = (my_sim.data['SN'] == my_sim.data['SN0'])
nphot = len(numpy.unique(my_sim.data['SN0']))


my_sim.plot_radial_profile(dermer = True, bins = 100)

plt.show()





#my_sim.coordinate_cut()
#my_sim.plot_map()
#plt.show()
