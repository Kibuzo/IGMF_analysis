from IO import _crpropa
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
from math_.base import radians_to_degree
import logging
toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'
#toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/1e-6/BPL'


my_sim = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-08')
non_interacting = (my_sim.data['SN'] == my_sim.data['SN0'])
nphot = len(numpy.unique(my_sim.data['SN0']))
logging.info(f'The simulations contains a total of {nphot} source photons')
logging.info(f'out of which, a total of {numpy.sum(non_interacting)} '
             'have not interacted')
theta_x = radians_to_degree(my_sim.data['Px'])
theta_y = radians_to_degree(my_sim.data['Py'])
theta_r = numpy.sqrt(theta_x**2+theta_y**2)
x_arr = my_sim.data['X']
y_arr = my_sim.data['Y']
r_arr = numpy.sqrt(x_arr**2+y_arr**2)
filter_p = numpy.logical_and(numpy.abs(theta_x)<2, numpy.abs(theta_y)<2)
filter_x = numpy.logical_and(numpy.abs(x_arr)<2, numpy.abs(y_arr)<2)

plt.figure ('Emission spectrum')
plt.hist(my_sim.data['E0'])
plt.yscale('log')

plt.figure ('Angular distribution')
plt.hist2d(theta_x[filter_p], theta_y[filter_p], bins = 100, 
           norm = mpl.colors.LogNorm())
plt.xlabel('$\\theta_{x}$[degrees]')
plt.ylabel('$\\theta_{y}$[degrees]')
plt.colorbar()

plt.figure('Spatial Distribution')
plt.hist2d(x_arr[filter_x], y_arr[filter_x], bins = 100, 
           norm = mpl.colors.LogNorm())
plt.xlabel('$\\theta_{x}$[degrees]')
plt.ylabel('$\\theta_{y}$[degrees]')
plt.colorbar()

my_sim.plot_radial_profile(dermer = True, bins = 100)
plt.show()

plt.show()





#my_sim.coordinate_cut()
#my_sim.plot_map()
#plt.show()
