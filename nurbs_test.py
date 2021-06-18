import sys
import bpy
import numpy as np
from colorsys import hsv_to_rgb

sys.path.append('/home/test/PycharmProjects/chromoSnake/')
import colorGenerator as cG

print('Starting sim loading...')
mySim = cG.ChromoSim('/home/test/Documents/chromoShake/varLp/5nm_trim.out')
print('Sim loaded')

top_mat_list = []
bottom_mat_list = []
for i in range(1, 17):
    h = 0.4 + (0.5 * (i % 2)) + (1 / 15) * i
    rgb = hsv_to_rgb(h, 0.85, 1)
    rgba = (rgb[0], rgb[1], rgb[2], 1)
    mat = bpy.data.materials.new('Material.' + str(i))
    mat.diffuse_color = rgba
    top_mat_list.append(mat)
    # complementary color section
    comp_h = h + 0.5
    comp_rgb = hsv_to_rgb(comp_h, 0.85, 1)
    comp_rgba = (comp_rgb[0], comp_rgb[1], comp_rgb[2], 1)
    comp_mat = bpy.data.materials.new('Comp_Material.' + str(i))
    comp_mat.diffuse_color = comp_rgba
    bottom_mat_list.append(comp_mat)


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


# keep labels for later
label_list = []
# iterate through all position labels
for i in range(1, 17):
    for position in ['top', 'bottom']:
        my_label = 'chr{0}_{1}'.format(str(i), position)
        label_list.append(my_label)
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

        if 'top' in my_label:
            h = 0.4 + (0.5 * (i % 2)) + (1 / 15) * i
        elif 'bottom' in my_label:
            h = (0.4 + (0.5 * (i % 2)) + (1 / 15) * i) + 0.5  # for comp color

        rgb = hsv_to_rgb(h, 0.85, 1)
        rgba = (rgb[0], rgb[1], rgb[2], 1)
        mat = bpy.data.materials.new('Material.' + str(i))
        mat.diffuse_color = rgba
        bpy.data.objects[my_label].active_material = mat

# iterate over time to create animation
scene = bpy.context.scene
scene.frame_start = 1
scene.frame_end = len(mySim.sim_dict)

for i, key in enumerate(mySim.sim_dict.keys()):
    scene.frame_set(i)
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

            point.keyframe_insert(data_path='co')

    print('Frame {0} completed...'.format(i))
