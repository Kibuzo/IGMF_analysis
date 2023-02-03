from IO import _crpropa
import numpy

toy_sim_dir='/media/kibuzo/80a7f701-fd11-4c0c-993a-76b511ae8b86/BackupHESS/BPL/Cone_final/TEMP'

my_sim = _crpropa.simulation(10, toy_sim_dir, '2data_final_largecone', '1e-08')
print (len(my_sim.data['E']))
