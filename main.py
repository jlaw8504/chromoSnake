from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import tensionSim as tS
import colorGenerator as cG
from datetime import datetime
import numpy as np

start = datetime.now()
colorSim = cG.ChromoSim('/home/test/Documents/chromoShake/varLp/50nm_trim.out')
colorSim.mass_labels.reverse()
colorSim.loop_labels.reverse()
label_bool_list = ['chr3_top' in label for label in colorSim.mass_labels]
cohesin_bool_list = ['cohesin' in label for label in colorSim.mass_labels]
label_bool_array = np.array(label_bool_list)
cohesin_loop_array = np.array(cohesin_bool_list)
loop_array = np.array(colorSim.loop_labels)
axis_array = np.logical_and(np.logical_not(loop_array), np.logical_not(cohesin_loop_array))
idx = axis_array
print(datetime.now() - start)
start = datetime.now()
mySim = tS.TensionSim('/home/test/Documents/chromoShake/varLp/50nm_trim.out')
print(datetime.now() - start)

# create the 3d axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# create list of time keys without struct
time_keys_list = list(mySim.sim_dict.keys())
time_keys_list = time_keys_list[1:]

# first frame
my_list_list = [list(mySim.sim_dict[time_keys_list[0]][mass_key].values()) for mass_key in
                mySim.sim_dict[time_keys_list[0]]]
my_array = np.array(my_list_list, dtype=float)


def plot_frame(time_key):
    plot_list_list = [list(mySim.sim_dict[time_key][mass_key].values()) for mass_key in mySim.sim_dict[time_key]]
    plot_array = np.array(plot_list_list, dtype=float)
    my_scatter._offsets3d = (plot_array[idx, 0], plot_array[idx, 1], plot_array[idx, 2])
    plot_color_array = scale_array(plot_array[idx, 3])
    my_scatter.set_array(plot_color_array)
    print(time_key)
    return my_scatter


def scale_array(array):
    sub_array = array - array.min()
    return np.round(sub_array / sub_array.max() * 255)


color_array = scale_array(my_array[idx, 3])
my_scatter = ax.scatter(my_array[idx, 0], my_array[idx, 1], my_array[idx, 2], cmap='inferno', c=color_array)
ani = FuncAnimation(fig, plot_frame, frames=time_keys_list, interval=30)
ani.save('50nm_trim_axis.mp4')
