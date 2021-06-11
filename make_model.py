import sys
import numpy as np
import bpy
from colorsys import hsv_to_rgb
sys.path.append('/home/test/PycharmProjects/chromoSnake')
import colorGenerator as cg

mySim = cg.ChromoSim('/home/test/Documents/chromoShake/bin/Video/WTspindle_0_5_first.out')

# create materials for all chromatid types
top_mat_list = []
bottom_mat_list= []
for i in range(1,17):
    h = 0.4 + (1/15)*(i-1)
    rgb = hsv_to_rgb(h, 0.8, 1)
    rgba = (rgb[0], rgb[1], rgb[2], 1)
    mat = bpy.data.materials.new('Material.' + str(i))
    mat.diffuse_color = rgba
    top_mat_list.append(mat)
    #complementary color section
    comp_h = h + 0.5
    comp_rgb = hsv_to_rgb(comp_h, 0.8, 1)
    comp_rgba = (comp_rgb[0], comp_rgb[1], comp_rgb[2], 1)
    comp_mat = bpy.data.materials.new('Comp_Material.' + str(i))
    comp_mat.diffuse_color = comp_rgba
    bottom_mat_list.append(comp_mat)

# create super and cohesin mats
super_mat = bpy.data.materials.new('super')
super_mat.diffuse_color = (0,0,0,1) #black
cohesin_mat = bpy.data.materials.new('cohesin')
cohesin_mat.diffuse_color = (1,1,1,1) #white

for label, mass in zip(mySim.mass_labels, mySim.mass_list):
    mass_location = np.array(mass[3:6]).astype(np.float) * 10**6
    bpy.ops.mesh.primitive_ico_sphere_add(radius=0.0055, location=mass_location)
    sphere = bpy.context.active_object
    if label == 'super':
        sphere.active_material = super_mat
    elif label == 'cohesin':
        label == 'cohesin'
        sphere.active_material = cohesin_mat
    else:
        label_idx = int(label.split('chr')[-1].split('_')[0])
        label_pos = label.split('_')[-1]
        if label_pos == 'top':
            sphere.active_material = top_mat_list[label_idx - 1]
        elif label_pos == 'bottom':
            sphere.active_material = bottom_mat_list[label_idx - 1]
        else:
            raise Exception('Label position: ' + label_pos + ' not found')
    print('Generating mass: ' + mass[1])
