import sys
import numpy as np
import bpy
from colorsys import hsv_to_rgb
sys.path.append('/home/test/PycharmProjects/chromoSnake')
import colorGenerator as cg

print('Loading Sim')
mySim = cg.ChromoSim('/home/test/Documents/chromoShake/varLp/5nm_trim.out')
time_dict = mySim.sim_dict['0']
print('Timepoint chosen')
del(mySim)
print('Sim deleted')
# create materials for all chromatid types
top_mat_list = []
bottom_mat_list= []
for i in range(1,17):
    h = 0.4 + (0.5 * (i %2)) + (1/15)*i
    rgb = hsv_to_rgb(h, 0.85, 1)
    rgba = (rgb[0], rgb[1], rgb[2], 1)
    mat = bpy.data.materials.new('Material.' + str(i))
    mat.diffuse_color = rgba
    top_mat_list.append(mat)
    #complementary color section
    comp_h = h + 0.5
    comp_rgb = hsv_to_rgb(comp_h, 0.85, 1)
    comp_rgba = (comp_rgb[0], comp_rgb[1], comp_rgb[2], 1)
    comp_mat = bpy.data.materials.new('Comp_Material.' + str(i))
    comp_mat.diffuse_color = comp_rgba
    bottom_mat_list.append(comp_mat)

print('Materials made')
# create super and cohesin mats
super_mat = bpy.data.materials.new('super')
super_mat.diffuse_color = (0,0.3,0,1) #dark green
cohesin_mat = bpy.data.materials.new('cohesin')
cohesin_mat.diffuse_color = (1,1,1,1) #white
print('Iterating over masses')
for mass_key in time_dict:
    mass_location = [
        time_dict[mass_key]['x'],
        time_dict[mass_key]['y'],
        time_dict[mass_key]['z']
    ]
    mass_location = np.array(list(map(np.float, mass_location))) * 10**6
    bpy.ops.mesh.primitive_ico_sphere_add(radius=0.0055, location=mass_location)
    sphere = bpy.context.active_object
    if time_dict[mass_key]['label'] == 'super':
        sphere.active_material = super_mat
    elif time_dict[mass_key]['label'] == 'cohesin':
        sphere.active_material = cohesin_mat
    else:
        label_idx = int(time_dict[mass_key]['label'].split('chr')[-1].split('_')[0])
        label_pos = time_dict[mass_key]['label'].split('_')[-1]
        if label_pos == 'top':
            sphere.active_material = top_mat_list[label_idx - 1]
        elif label_pos == 'bottom':
            sphere.active_material = bottom_mat_list[label_idx - 1]
        else:
            raise Exception('Label position: ' + label_pos + ' not found')
    print('Generating mass: ' + str(mass_key))
