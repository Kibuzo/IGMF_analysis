from IO import _crpropa
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
from math_.base import radians_to_degree
import logging
import os
toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'
#toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/1e-6/BPL'

fig_dir = '/home/kibuzo/TesiPlot'
sim9 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-09')
sim8 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-08')
sim7 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-07')
sim6 = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-06')



SIM9D = sim9.plot_radial_profile(dermer = True, bins = 100)
plt.close()
SIM8D = sim8.plot_radial_profile(dermer = True, bins = 100)
plt.close()
SIM7D = sim7.plot_radial_profile(dermer = True, bins = 100)
plt.close()
SIM6D = sim6.plot_radial_profile(dermer = True, bins = 100)
plt.close()
SIM9 = sim9.plot_radial_profile(dermer = False, bins = 100)
plt.close()
SIM8 = sim8.plot_radial_profile(dermer = False, bins = 100)
plt.close()
SIM7 = sim7.plot_radial_profile(dermer = False, bins = 100)
plt.close()
SIM6 = sim6.plot_radial_profile(dermer = False, bins = 100)
plt.close()

plt.figure(f'Radial profile (1e-9 G)')
plt.plot(SIM9[1], SIM9[0], label = 'momenta')
plt.plot(SIM9D[1], SIM9D[0], label = 'dermer')
plt.title(f'Radial profile (1e-9 G)')
plt.xlabel ('Radial distange (degrees)')
plt.ylabel ('Area-weighted counts')
plt.legend()
plt.yscale('log')
plt.savefig (os.path.join(fig_dir,'1e-9.png'))

plt.figure(f'Radial profile (1e-8 G)')
plt.plot(SIM8[1], SIM8[0], label = 'momenta')
plt.plot(SIM8D[1], SIM8D[0], label = 'dermer')
plt.title(f'Radial profile (1e-8 G)')
plt.xlabel ('Radial distange (degrees)')
plt.ylabel ('Area-weighted counts')
plt.legend()
plt.yscale('log')
plt.savefig (os.path.join(fig_dir,'1e-8.png'))


plt.figure(f'Radial profile (1e-7 G)')
plt.plot(SIM7[1], SIM7[0], label = 'momenta')
plt.plot(SIM7D[1], SIM7D[0], label = 'dermer')
plt.title(f'Radial profile (1e-7 G)')
plt.xlabel ('Radial distange (degrees)')
plt.ylabel ('Area-weighted counts')
plt.legend()
plt.yscale('log')
plt.savefig (os.path.join(fig_dir,'1e-7.png'))


plt.figure(f'Radial profile (1e-6 G)')
plt.plot(SIM6[1], SIM6[0], label = 'momenta')
plt.plot(SIM6D[1], SIM6D[0], label = 'dermer')
plt.title(f'Radial profile (1e-6 G)')
plt.xlabel ('Radial distange (degrees)')
plt.ylabel ('Area-weighted counts')
plt.legend()
plt.yscale('log')
plt.savefig (os.path.join(fig_dir,'1e-6.png'))


plt.show()

