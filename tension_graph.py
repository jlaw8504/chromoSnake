import colorGenerator as cG
import tensionSim as tS
import numpy as np
import matplotlib.pyplot as plt

file_list = ['WTspindle_0_5.out', 'WTspindle_1_5.out']
size_list = ['6', '15']
for filename, loop_size in zip(file_list, size_list):
    print('Parsing {0} ...'.format(filename))
    file_string = '/home/test/Documents/chromoShake/varSize/{0}'.format(filename)
    colorSim = cG.ChromoSim(file_string)
    colorSim.mass_labels.reverse()
    colorSim.loop_labels.reverse()
    label_bool_list = ['chr3_top' in label for label in colorSim.mass_labels]
    cohesin_bool_list = ['cohesin' in label for label in colorSim.mass_labels]
    super_bool_list = ['super' in label for label in colorSim.mass_labels]
    label_bool_array = np.array(label_bool_list)
    cohesin_loop_array = np.array(cohesin_bool_list)
    super_bool_array = np.array(super_bool_list)
    loop_array = np.array(colorSim.loop_labels)
    axis_array = np.logical_and(np.logical_not(loop_array), np.logical_not(cohesin_loop_array))
    idx = np.logical_and(axis_array, np.logical_not(super_bool_array))
    mySim = tS.TensionSim(file_string)

    # built a tension and idx nd-array with time
    tension_array = np.zeros(
        (
            np.sum(idx),  # number of masses to parse
            2,  # idx and tension
            len(mySim.sim_dict) - 1  # time-points, excluding struct
        )
    )
    for time_idx_sub_me, time_key in enumerate(mySim.sim_dict):
        if time_key == 'struct':
            continue
        time_idx = time_idx_sub_me - 1
        for mass_idx, mass_key in enumerate(list(np.where(idx)[0])):
            # must correct time_idx for excluded 'struct' dict
            tension_array[mass_idx, 0, time_idx] = float(mass_key)
            tension_array[mass_idx, 1, time_idx] = mySim.sim_dict[time_key][str(mass_key)]['tension']
    # mySim is huge so delete it and delete colorSim for good measure
    del mySim
    del colorSim
    plt.plot(np.mean(tension_array[:, :, 2001:], axis=2)[:, 1])
    plt.title('Mean Tension, Final 1000 Timepoints, loop size = {0} kb'.format(loop_size))
    plt.xlabel('Bead Index Position on Axis')
    plt.ylabel('Tension (N)')
    plt.savefig('{0}kb_axis_only_final_1000.png'.format(loop_size))
    plt.close('all')
