import sys
import bpy
import numpy as np
sys.path.append('/home/test/PycharmProjects/chromoSnake/')
import colorGenerator as cg

print('Staring sim loading...')
mySim = cg.ChromoSim('/home/test/Documents/chromoShake/varLp/5nm_trim.out')
print('Sim loaded')
time_dict = mySim.sim_dict['0']
print('Time selected')
del(mySim)
print('Sim cleared')

# collect all mass points for chr1_top
first_array = np.empty(shape=(0,3))
second_array = np.empty(shape=(0,3))
# specify first or second strand
strand_idx = 0
# going to have keep track of sides of C-loop
in_label = False
dist = -1


def parse_mass_coords(time_dict, mass_key):
    mass_coords = [
        time_dict[mass_key]['x'],
        time_dict[mass_key]['y'],
        time_dict[mass_key]['z']
    ]
    mass_array = np.array(list(map(np.float, mass_coords)))
    mass_array = mass_array * 10 ** 6
    return mass_array


for mass_key in time_dict:
    if time_dict[mass_key]['label'] == 'chr1_top':
        if in_label is False: #this is first coord
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
            mass_array = mass_array * 10**6
            dist = np.linalg.norm(mass_array - prev_array)
            if dist > 0.02:
                strand_idx = 1
        
        print(str(dist))
        print(time_dict[mass_key]['label'] + '_strand_' + str(strand_idx))
        
        if strand_idx == 0:
            first_array = np.vstack((first_array, mass_array))
        else:
            second_array = np.vstack((second_array, mass_array))
        
        prev_array = mass_array

# now need to split the coord_array on the third axis
second_array = np.flipud(second_array)
coord_array = np.vstack((first_array, second_array))
coord_list = coord_array.tolist()
# make a new curve
crv = bpy.data.curves.new('crv', 'CURVE')
crv.dimensions = '3D'

# make a new spline in that curve
spline = crv.splines.new(type='NURBS')

# a spline point for each point
spline.points.add(len(coord_list)-1) # theres already one point by default

# assign the point coordinates to the spline points
for p, new_co in zip(spline.points, coord_list):
    p.co = (new_co + [1.0]) # (add nurbs weight)

# make a new object with the curve
obj = bpy.data.objects.new('object_name', crv)
bpy.context.scene.collection.objects.link(obj)
