from IO import _crpropa
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
from math_.base import radians_to_degree
import logging
toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'
#toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/1e-6/BPL'

sim9 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-09')
sim8 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-08')
sim7 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-07')
sim6 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-06')



sim9.plot_radial_profile(dermer = True, bins = 100)
sim8.plot_radial_profile(dermer = True, bins = 100)
sim7.plot_radial_profile(dermer = True, bins = 100)
sim6.plot_radial_profile(dermer = True, bins = 100)
sim9.plot_radial_profile(dermer = False, bins = 100)
sim8.plot_radial_profile(dermer = False, bins = 100)
sim7.plot_radial_profile(dermer = False, bins = 100)
sim6.plot_radial_profile(dermer = False, bins = 100)


plt.show()





#my_sim.coordinate_cut()
#my_sim.plot_map()
#plt.show()
