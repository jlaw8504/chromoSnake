import colorGenerator as cG
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# noinspection SpellCheckingInspection
mySim = cG.ChromoSim('/home/test/Documents/chromoShake/varLp/5nm_trim.out')
# dictionary method
cohesin_array = np.empty((0,3), dtype=np.float)
time_dict = mySim.sim_dict['0']
for mass_idx in time_dict:
    if time_dict[mass_idx]['label'] == 'chr1_top':
        mass_coords = [
            time_dict[mass_idx]['x'],
            time_dict[mass_idx]['y'],
            time_dict[mass_idx]['z'],
        ]
        mass_coords = list(map(np.float, mass_coords))
        mass_array = np.array(mass_coords) * 10**9
        cohesin_array = np.vstack((cohesin_array, mass_array))
cohesin_array[np.argsort(cohesin_array[:, 2], axis=0), :]
cohesin_array[np.argsort(cohesin_array[:, 1], axis=0), :]
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(cohesin_array[:, 0], cohesin_array[:, 1], cohesin_array[:, 2])
# ax.scatter(mySim.cohesin_masses_array[:,3,:].flatten().astype(np.double)*1e9,
#     mySim.cohesin_masses_array[:,4,:].flatten().astype(np.double)*1e9,
#     mySim.cohesin_masses_array[:,5,:].flatten().astype(np.double)*1e9
#     )
plt.show()
