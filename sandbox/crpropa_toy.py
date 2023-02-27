from IO import _crpropa
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
from math_.base import radians_to_degree

toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'
#toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/1e-6/BPL'


my_sim = _crpropa.simulation(10000, toy_sim_dir, '2data_final_largecone', '1e-09')

theta_x = radians_to_degree(my_sim.data['Px'])
theta_y = radians_to_degree(my_sim.data['Py'])
theta_r = numpy.sqrt(theta_x**2+theta_y**2)
x_arr = my_sim.data['X']
y_arr = my_sim.data['Y']
r_arr = numpy.sqrt(x_arr**2+y_arr**2)
filter_p = numpy.logical_and(numpy.abs(theta_x)<2, numpy.abs(theta_y)<2)
filter_x = numpy.logical_and(numpy.abs(x_arr)<2, numpy.abs(y_arr)<2)

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

plt.figure ('radial plot')
plt.hist(theta_r[filter_p], weights=1./theta_r[filter_p], bins = 100)
r_max = numpy.max(theta_r[filter_p]/numpy.sqrt(2))
plt.xlim(0,r_max)
plt.yscale('log')
plt.ylabel('Area-normalized counts')
plt.xlabel('Angular distance taken from beam center[degrees]')

plt.figure ('radial plot[coordinates]')
plt.hist(r_arr[filter_x], weights=1./r_arr[filter_x], bins = 100)
r_max = numpy.max(r_arr[filter_x]/numpy.sqrt(2))
plt.xlim(0,r_max)
plt.yscale('log')
plt.ylabel('Area-normalized counts')
plt.xlabel('Arrival coordinate distribution taken from beam center[degrees]')

plt.figure('Radial profile (Dermer plot)')
radial = (radians_to_degree(my_sim.dermer_theta))
plt.hist(radial, weights=1./radial, bins=1000)
plt.yscale('log')
plt.xlabel('Angular distance [degrees]')
plt.show()



plt.show()


#my_sim.coordinate_cut()
#my_sim.plot_map()
#plt.show()
