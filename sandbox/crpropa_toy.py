from IO import _crpropa
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
from math_.base import radians_to_degree

toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'

my_sim = _crpropa.simulation(1000, toy_sim_dir, '2data_final_largecone', '1e-08')
theta_x = radians_to_degree(my_sim.data['Px'])
theta_y = radians_to_degree(my_sim.data['Py'])
filter = numpy.logical_and(numpy.abs(theta_x)<1, numpy.abs(theta_y)<1)
plt.hist2d(theta_x[filter], theta_y[filter], bins = 100, norm = mpl.colors.LogNorm())
plt.colorbar()
plt.show()


#my_sim.coordinate_cut()
#my_sim.plot_map()
#plt.show()
