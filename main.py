import colorGenerator as cG
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# noinspection SpellCheckingInspection
mySim = cG.ChromoSim('/home/test/Documents/chromoShake/bin/Video/WTspindle_0_5_first.out')
mass_array = np.array(mySim.mass_list)
spring_array = np.array(mySim.spring_list)
hinge_array = np.array(mySim.hinge_list)
super_mass_indexes = mySim.super_bool_array
condensin_spring_indexes = mySim.condensin_spring_indexes
condensin_array = spring_array[condensin_spring_indexes, :]
mba = np.logical_or(np.array(mySim.mass_labels) == 'chr1_top', np.array(mySim.mass_labels) == 'super')
mba = np.logical_or(mba, np.array(mySim.mass_labels) == 'chr2_top')
print(sum(mba))

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(mass_array[mba, 3].astype(np.double)*1e9, mass_array[mba, 4].astype(np.double)*1e9, mass_array[mba, 5].astype(np.double)*1e9)
# ax.scatter(mySim.cohesin_masses_array[:,3,:].flatten().astype(np.double)*1e9,
#     mySim.cohesin_masses_array[:,4,:].flatten().astype(np.double)*1e9,
#     mySim.cohesin_masses_array[:,5,:].flatten().astype(np.double)*1e9
#     )
plt.show()
