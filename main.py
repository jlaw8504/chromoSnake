import colorGenerator as cG
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# noinspection SpellCheckingInspection
mySim = cG.ChromoSim('/home/test/Documents/chromoShake/bin/Video/WTspindle_0_5_first.out')
mass_array = np.array(mySim.mass_list)
spring_array = np.array(mySim.spring_list)
hinge_array = np.array(mySim.hinge_list)
super_mass_indexes = mySim.super_mass_indexes
condensin_spring_indexes = mySim.condensin_spring_indexes
condensin_array = spring_array[condensin_spring_indexes, :]
print(condensin_array.shape)

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(mass_array[:, 3].astype(np.double)*1e9, mass_array[:, 4].astype(np.double)*1e9, mass_array[:, 5].astype(np.double)*1e9)
ax.scatter(mass_array[0, 3].astype(np.double)*1e9, mass_array[0, 4].astype(np.double)*1e9, mass_array[0, 5].astype(np.double)*1e9)
plt.show()
