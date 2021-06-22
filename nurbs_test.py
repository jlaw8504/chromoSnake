# HARD-CODED VARIABLES
chromoSnake_dir_path = ''
simulation_file_path = ''
render_dir_path = ''

import sys
sys.path.append(chromoSnake_dir_path)
import os
import bpy
import numpy as np
from colorsys import hsv_to_rgb
import colorGenerator as cG


def parse_mass_coords(time_dict, mass_key):
    mass_coords = [
        time_dict[mass_key]['x'],
        time_dict[mass_key]['y'],
        time_dict[mass_key]['z']
    ]
    mass_array = np.array(list(map(np.float, mass_coords)))
    mass_array = mass_array * 10 ** 6
    return mass_array


def parse_coord_list(sim, time_str, label_str):
    # collect all mass points for chr1_top
    first_array = np.empty(shape=(0, 3))
    second_array = np.empty(shape=(0, 3))
    # specify first or second strand
    strand_idx = 0
    # going to have keep track of sides of C-loop
    in_label = False
    dist = -1
    time_dict = sim.sim_dict[time_str]
    for mass_key in time_dict:
        if time_dict[mass_key]['label'] == label_str:
            if in_label is False:  # this is first coord
                in_label = True
                dist = 0
                mass_array = parse_mass_coords(time_dict, mass_key)
                prev_array = mass_array

            if in_label:
                mass_coords = [
                    time_dict[mass_key]['x'],
                    time_dict[mass_key]['y'],
                    time_dict[mass_key]['z']
                ]
                mass_array = np.array(list(map(np.float, mass_coords)))
                mass_array = mass_array * 10 ** 6
                dist = np.linalg.norm(mass_array - prev_array)
                if dist > 0.02:
                    strand_idx = 1

            if strand_idx == 0:
                first_array = np.vstack((first_array, mass_array))
            else:
                second_array = np.vstack((second_array, mass_array))

            prev_array = mass_array
    # now need to split the coord_array on the third axis
    if 'top' in my_label:
        second_array = np.flipud(second_array)
    elif 'bottom' in my_label:
        first_array = np.flipud(first_array)

    coord_array = np.vstack((first_array, second_array))
    coord_list = coord_array.tolist()
    return coord_list


# parse simname
sim_path_tuple = os.path.split(simulation_file_path)
sim_name_tuple = os.path.splitext(sim_path_tuple[-1])
sim_basename = sim_name_tuple[0]

# import simulation
print('Starting sim loading...')
mySim = cG.ChromoSim(simulation_file_path)
print('Sim loaded')

# generate labels
dna_label_list = []
# iterate through all position labels
for i in range(1, 17):
    for position in ['top', 'bottom']:
        my_label = 'chr{0}_{1}'.format(str(i), position)
        dna_label_list.append(my_label)

coh_label_list = [mySim.sim_dict['0'][mass_key]['label']
                  for mass_key in mySim.sim_dict['0']
                  if 'cohesin' in mySim.sim_dict['0'][mass_key]['label']]

# append cohesin complex labels to label_list
coh_label_set = set(coh_label_list)
label_list = dna_label_list + list(coh_label_set)

# make cohesin material
coh_mat = bpy.data.materials.new('Material.cohesin')
coh_mat.diffuse_color = (1, 1, 1, 1)

for my_label in label_list:
    my_coord_list = parse_coord_list(mySim, '0', my_label)
    # make a new curve
    crv = bpy.data.curves.new(my_label, 'CURVE')
    crv.dimensions = '3D'

    # make a new spline in that curve
    spline = crv.splines.new(type='NURBS')

    # a spline point for each point
    spline.points.add(len(my_coord_list) - 1)  # theres already one point by default

    # assign the point coordinates to the spline points
    for p, new_co in zip(spline.points, my_coord_list):
        p.co = (new_co + [1.0])  # (add nurbs weight)

    # make a new object with the curve
    obj = bpy.data.objects.new(my_label, crv)
    bpy.context.scene.collection.objects.link(obj)
    bpy.data.objects[my_label].data.bevel_depth = 0.005
    # generate and assign material
    if 'cohesin' in my_label:
        bpy.data.objects[my_label].active_material = coh_mat
        continue

    if 'top' in my_label:
        idx = int(my_label.split('chr')[-1].split('_')[0])
        h = 0.4 + (0.5 * (idx % 2)) + (1 / 15) * idx
    elif 'bottom' in my_label:
        idx = int(my_label.split('chr')[-1].split('_')[0])
        h = (0.4 + (0.5 * (idx % 2)) + (1 / 15) * idx) + 0.5  # for comp color

    rgb = hsv_to_rgb(h, 0.85, 1)
    rgba = (rgb[0], rgb[1], rgb[2], 1)
    mat = bpy.data.materials.new('Material.' + str(i))
    mat.diffuse_color = rgba
    bpy.data.objects[my_label].active_material = mat

for i, key in enumerate(mySim.sim_dict.keys()):
    time_dict = mySim.sim_dict[key]
    for my_label in label_list:
        # select existing spline in blender
        crv = bpy.data.curves[my_label]
        spline = crv.splines[0]
        points = spline.points
        # parse corresponding points in sim_dict
        my_coord_list = parse_coord_list(mySim, str(key), my_label)
        for point, coord in zip(points, my_coord_list):
            for idx in range(3):
                point.co[idx] = coord[idx]
    bpy.context.scene.render.filepath = render_dir_path + sim_basename + '_' + str(i).zfill(4)
    bpy.ops.render.render(write_still=True)
