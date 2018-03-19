
# coding: utf-8

# In[1]:

def number_of_masses(outfile):
    """Returns the number of masses in an outfile"""
    import re
    
    cnt = 0
    massexp = re.compile('mass ') #marks the mass lines for counting
    with open(outfile) as f:
        for line in f:
            if massexp.search(line):
                cnt = cnt + 1
    f.close()
    return cnt

def super_masses(outfile):
    """Return a list of the super massive beads from chromoShake pericentromere simulation outfiles"""
    import re
    
    mass_nums = []
    regexp = re.compile('3.38889e-015') #mass of the super massive beads
    endexp = re.compile('}') #this is the end of the section we will need to parse
    cnt = 0
    massexp = re.compile('mass ') #marks the mass lines for counting
    with open(outfile) as file_object:
        for line in file_object:
            if endexp.search(line):
                break
            elif regexp.search(line):
                mass_list = line.split()
                mass_nums.append(mass_list[1])
            elif massexp.search(line):
                cnt = cnt + 1
    file_object.close()
    return mass_nums

def super_mass_coords(outfile):
    """Returns a numpy array of the super massive bead coordinates"""
    import numpy as np
    import re
    
    regexp = re.compile('3.38889e-015') #mass of the super massive beads
    endexp = re.compile('}') #this is the end of the section we will need to parse
    mass_count = len(super_masses(outfile)) #number of super massive beads
    coord_array = np.empty([mass_count, 3]) #random array that we will index to and then return
    row = 0 #counter that increments as for loop progresses
    with open(outfile) as file_object:
        for line in file_object:
            if endexp.search(line):
                break
            if regexp.search(line):
                mass_list = line.split()
                for col in range(3):
                    coord_array[row][col] = mass_list[col+3]
                row = row + 1               
    return coord_array

def super_mass_indexs(config_file):
    """Returns a list of the super massive bead indexes"""
    import re
    
    regexp = re.compile('3.38889e-015') #mass of the super massive beads
    endexp = re.compile('}') #this is the end of the section we will need to parse
    coord_idxs = list()
    with open(config_file) as file_object:
        for line in file_object:
            if endexp.search(line):
                break
            if regexp.search(line):
                mass_list = line.split()
                coord_idxs.append(mass_list[1])
    return coord_idxs

def labels_from_end(mass_nums, increment):
    """Create a new mass color list by indicating how many beads from the end you want to label"""
    alpha_list = list() #mass indexs that start each strand
    omega_list = list() #mass indexs that end each strand
    for idx, value in enumerate(mass_nums):
        if idx % 2 == 0:
            alpha_list.append(value)
        else:
            omega_list.append(value)
    alpha_array = [ int(x) for x in alpha_list ]
    omega_array = [ int(x) for x in omega_list ]
    alpha_new = [ x + increment for x in alpha_array]
    omega_new = [ x - increment for x in omega_array]
    new_mass_nums = alpha_new + omega_new
    new_mass_nums = [ str(x) for x in sorted(new_mass_nums) ]
    return new_mass_nums

def recolor_outfile(outfile, increment, filename):
    """Outputs a recolored centromere simulation based on how many beads from the end you want labeled"""
    import re
    
    mass_nums = super_masses(outfile)
    cnt = number_of_masses(outfile)
    new_mass_nums = labels_from_end(mass_nums, increment)
    color_list = list();
    for number in range(cnt):
        if str(number) in new_mass_nums:
            color_list.append(5)
        else:
            color_list.append(1)
    sub_toggle = False
    time_toggle = True
    timeexp = re.compile('Time')
    rev_color_list = list(reversed(color_list))
    colorexp = re.compile('MassColors')
    outfile = open(filename,'w')
    with open('default_spindle.out') as infile:
        for line in infile:
            if sub_toggle:
                for number in rev_color_list:
                    outfile.write(str(number)+'\n')
                outfile.write('\n')
                sub_toggle = False
                time_toggle = False
            elif timeexp.search(line):
                time_toggle = True
                outfile.write(line)
            elif colorexp.search(line):
                sub_toggle = True
                outfile.write(line)
            elif time_toggle:
                outfile.write(line)
    outfile.close()

def src_generator(outfile, increment, SRC_filename):
    """Outputs a someRandomColor file in the correct order for the microscope simulator"""
    mass_nums = super_masses(outfile)
    new_mass_nums = labels_from_end(mass_nums, increment)
    cnt = number_of_masses(outfile)
    color_list = list();
    for number in range(cnt):
        if str(number) in new_mass_nums:
            color_list.append(4)
        else:
            color_list.append(1)
    rev_color_list = list(reversed(color_list))
    src = open(SRC_filename,'w')
    for number in rev_color_list:
        src.write(str(number)+'\n')
    src.close()

def outfile_convert(outfile, convert_filename):
    """Converts the output coordinates to microns and removes header portion of outfile"""
    import re
    
    cf = open(convert_filename,'w')
    time_toggle = False
    regexp = re.compile('Time ')
    with open(outfile) as f:
        for line in f:
            if time_toggle and not regexp.search(line) and line.strip():
                line_list = line.split()
                line_list = [ float(x)*1000000 for x in line_list]
                cf.write(str(line_list[0]) + ' ' + str(line_list[1]) + ' ' + str(line_list[2]) + '\n')
            elif time_toggle and not regexp.search(line) and not line.strip():
                cf.write(line)
            if regexp.search(line):
                time_toggle = True
                cf.write(line)
    cf.close()

def most_proximal_bead_list(coords):
    import numpy as np
    
    """Parse the array of super massive bead coordinates and return list of tuples of proximal beads"""
    #Loop through coordinate arrays to find most proximal coordinate
    idx_list = list() #list to keep track of found idxs to prevent repeat entries in tuple_list
    tuple_list = list() #place proximal coords here
    for idx in range(np.shape(coords)[0]):
        if not(idx in idx_list):
            x = np.subtract(coords, coords[idx])
            x_sq = np.square(x)
            sq_sum = np.sum(x_sq, axis=1)
            dist = np.sqrt(sq_sum)
            dist_del = np.delete(dist,idx)
            min_idx_del = np.argmin(dist_del)
            if min_idx_del >= idx:
                min_idx = min_idx_del + 1
            elif min_idx_del < idx:
                min_idx = min_idx_del
            tuple_list.append((idx,min_idx))
            idx_list.append(min_idx)
    return(tuple_list)

def new_centromere_positions(coords, tuple_list):
    import numpy as np
    
    """Parse a xyz coordinate array with a list of tuples that refer to indexes of 
    the super massive beads and returns an array of the new centromere bead positions"""
    #loop through tuple_list to generate average x and y coordinates and extend z-coordinate by 10 nm
    coords_shape = np.shape(coords)
    new_beads = np.empty([int(coords_shape[0]/2), coords_shape[1]])
    row = 0
    for pair in tuple_list:
        new_xy = np.divide(np.add(coords[pair[0],:2], coords[pair[1],:2]),2)
        old_z = coords[pair[0],-1]
        if old_z > 0:
            new_z = np.add(old_z, float(1e-08))
        else:
            new_z = np.subtract(old_z, float(1e-08))
        new_bead = np.append(new_xy, new_z)
        new_beads[row] = new_bead
        row = row + 1
    return(new_beads)

def new_kmt_positions(coords, tuple_list, distance):
    import numpy as np
    
    """Parse a xyz coordinate array with a list of tuples that refer to indexes of 
    the super massive beads and returns an array of kMT bead positions shifted by distance"""
    #loop through tuple_list to generate average x and y coordinates and extend z-coordinate by 
    coords_shape = np.shape(coords)
    new_beads = np.empty([int(coords_shape[0]/2), coords_shape[1]])
    row = 0
    for pair in tuple_list:
        new_xy = np.divide(np.add(coords[pair[0],:2], coords[pair[1],:2]),2)
        old_z = coords[pair[0],-1]
        if old_z > 0:
            new_z = np.add(old_z, distance)
        else:
            new_z = np.subtract(old_z, distance)
        new_bead = np.append(new_xy, new_z)
        new_beads[row] = new_bead
        row = row + 1
    return(new_beads)

def get_mass_separation(config_file):
    """Parse the configuration for the mass separation"""
    import re
    
    spring_re = re.compile("spring *")
    with open(config_file) as cf:
        for line in cf:
            if spring_re.search(line):
                split_line = line.split()
                if len(split_line) == 5:
                    mass_sep = split_line[3]
                    break
    return(mass_sep)

def get_spring_constant(config_file):
    """Parse the configuration file for the spring constant and return it as a string"""
    import re
    
    spring_re = re.compile("spring *")
    with open(config_file) as cf:
        for line in cf:
            if spring_re.search(line):
                split_line = line.split()
                if len(split_line) == 5:
                    ks = split_line[4]
                    break
    return(float(ks))

def centromere_spring_list(config_file):
    """Generate a master list containing two spring lines for every centromere mass"""
    #define number of existing masses as a new index
    mass_num = number_of_masses(config_file)
    #define centromere linking springs
    spring_lines_list = list()
    mass_sep = get_mass_separation(config_file)
    coords = super_mass_coords(config_file)
    ks = get_spring_constant(config_file)
    tuple_list = most_proximal_bead_list(coords)
    bead_idxs = super_mass_indexs(config_file)
    for pair in tuple_list:
        spring_lines = list()
        spring_lines.append('spring ' + str(bead_idxs[pair[0]]) + ' ' + str(mass_num) + ' ' + mass_sep + ' ' + str(ks) + ' \n')
        spring_lines.append('spring ' + str(bead_idxs[pair[1]]) + ' ' + str(mass_num) + ' ' + mass_sep + ' ' + str(ks) + ' \n')
        spring_lines_list.append(spring_lines)
        mass_num = mass_num + 1
    return(spring_lines_list)

def get_hinge_force(config_file):
    """Parse hinge force from a configuration file and return it as a string"""
    import re
    
    hinge_re = re.compile("hinge *")
    with open(config_file) as cf:
        for line in cf:
            if hinge_re.search(line):
                line_split = line.split()
                if len(line_split) == 5:
                    hinge_force = line_split[-1]
                    break
    return(hinge_force)
                
def super_partner_dict(config_file):
    """Returns a dictionary with super massive bead indices as keys and the
    index of the binding partner as a value. Used for making hinges in the
    centromere_hinge_list function"""
    import re
    
    bead_idxs = super_masses(config_file)
    spring_re = re.compile('^spring ')
    partner_idxs = list()
    with open(config_file) as f:
        for line in f:
            if spring_re.search(line.strip()):
                if line.split()[1] in bead_idxs:
                    partner_idxs.append(line.split()[2])
                elif line.split()[2] in bead_idxs:
                    partner_idxs.append(line.split()[1])
    partner_dict = dict(zip(bead_idxs, partner_idxs))
    return(partner_dict)

def centromere_hinge_list(config_file):
    """Generate a master list of two hinge lines for every centromere mass"""
    #define number of existing masses as a new index
    mass_num = number_of_masses(config_file)
    #define centromere linking hinges
    hinge_lines_list = list()
    coords = super_mass_coords(config_file)
    hinge_force = get_hinge_force(config_file) ## RETURN HERE
    tuple_list = most_proximal_bead_list(coords)
    bead_idxs = super_mass_indexs(config_file)
    partner_dict = super_partner_dict(config_file)
    for pair in tuple_list:
        hinge_lines = list()
        hinge_lines.append('hinge ' + str(bead_idxs[pair[0]]) + ' ' + partner_dict[bead_idxs[pair[0]]] + ' ' + str(mass_num) + ' ' + hinge_force + ' \n')
        hinge_lines.append('hinge ' + str(bead_idxs[pair[1]]) + ' ' + partner_dict[bead_idxs[pair[1]]] + ' ' + str(mass_num) + ' ' + hinge_force + ' \n')
        hinge_lines.append('hinge ' + str(bead_idxs[pair[0]]) + ' ' + str(bead_idxs[pair[1]]) + ' ' + str(mass_num) + ' ' + hinge_force + ' \n')
        hinge_lines_list.append(hinge_lines)
        mass_num = mass_num + 1
    return(hinge_lines_list)

def add_centromeres_spindle(config_file, altered_file, distance):
    """Parse a configuration file of a spindle simulation and create an altered configuration file with centromeres added"""
    import re
    
    cf = open(altered_file,'w')
    regexp = re.compile('}')
    mass_re = re.compile('mass *')
    #define number of existing masses as a new index
    mass_num = number_of_masses(config_file)
    spring_lines_list = centromere_spring_list(config_file)
    spring_cnt = 0
    kmt_spring_cnt = 0
    #define centromere bead postions
    coords = super_mass_coords(config_file)
    tuple_list = most_proximal_bead_list(coords)
    new_beads = new_centromere_positions(coords, tuple_list)
    #define hinge lines
    hinge_lines_list = centromere_hinge_list(config_file)
    #define kMT bead positions
    kmt_beads = new_kmt_positions(coords, tuple_list, distance)
    ks = get_spring_constant(config_file)
    #begin looping through config file
    with open(config_file) as f:
        for line in f:
            if mass_re.search(line): #remove all other supermassive beads
                line_list = line.split()
                if len(line_list) == 7:
                    cf.write(line_list[0] + ' ' + line_list[1] + '\t' + '3.38889e-020' + '\t'                             + line_list[3] + ' ' + line_list[4] + ' ' + line_list[5]                             + ' ' + line_list[6] + '\n')
                else:
                    cf.write(line)
            elif regexp.search(line): 
                    for bead in new_beads: # add in new centromeres
                        new_line = 'mass ' + str(mass_num) + '\t' + '3.38889e-020' + '\t'                        + '{:g}' + ' ' + '{:g}' + ' ' + '{:g}' + ' 5\n'
                        mass_num = mass_num + 1
                        cf.write(new_line.format(bead[0],bead[1],bead[2]))
                        #write spring and hinge lines
                        cf.write(spring_lines_list[spring_cnt][0])
                        cf.write(hinge_lines_list[spring_cnt][0])
                        cf.write(spring_lines_list[spring_cnt][1])
                        cf.write(hinge_lines_list[spring_cnt][1])
                        cf.write(hinge_lines_list[spring_cnt][2])
                        spring_cnt = spring_cnt + 1                    
                    for kmt_bead in kmt_beads: # add in the kMT beads
                        kmt_line = 'mass ' + str(mass_num) + '\t' + '3.38889e-015' + '\t'                        + '{:g}' + ' ' + '{:g}' + ' ' + '{:g}' + ' 3\n'
                        corr_spring_line = spring_lines_list[kmt_spring_cnt][0]
                        corr_split = corr_spring_line.split()
                        cen_idx = corr_split[2]
                        kmt_spring_line = 'spring ' + cen_idx + ' ' + str(mass_num) + ' ' +                         str(distance) + ' ' + str(ks) + '\n'
                        cf.write(kmt_line.format(kmt_bead[0],kmt_bead[1],kmt_bead[2]))
                        cf.write(kmt_spring_line)
                        mass_num = mass_num + 1
                        kmt_spring_cnt = kmt_spring_cnt + 1
                    cf.write(line) # writes the end bracket
            else:
                cf.write(line)
    cf.close()

def alter_kmt_lengths(outfile, distance_nm, newfile):
    """Parses the configuration portion of a chromoShake centromere simulation
    with kmt springs and alters the springs' rest lengths by given distance in nanometers"""
    #WARNING:This function assumes all the supermassive beads are kmt ends!!!
    #This function will use subprocess.Popen, wc, and grep, so your system needs grep
    #grep and wc for the timelines and coordinate section of an outfile
    #WARNING:This function uses grep -x and assumes that you are using newlines of \n (unix)
    #not \r\n (dos). If using chromoShake files from a Windows machine, convert
    #all newline characters to /n using dos2unix, or tr -d '\r' < input > output, or
    #perl -pi -e 's/\r\n/\n/g' input, or sed -i 's/^M$//' input
    
    import re
    import subprocess as sp
    
    #use system wc -l command to generate number of lines in file
    wc_line = sp.Popen(["wc", "-l", outfile], stdout=sp.PIPE).communicate()[0]
    wc = wc_line.decode('utf-8').split()[0]
    print('Number of lines: '+ wc + '\n')
    #find first timeline
    time_b = sp.Popen(["grep", "-m1", "Time ", outfile], stdout=sp.PIPE).communicate()[0]
    time = time_b.decode('utf-8')
    print('First timeline: ' + time.strip())
    #use system grep command to parse outfile for existing coordinates
    coords_b = sp.Popen(["grep", "-x", time.strip(), "-A", wc, outfile], stdout=sp.PIPE).communicate()[0]
    coords = coords_b.decode('utf-8').strip()
    #find supremassive bead indexes as list
    mass_idx = super_masses(outfile)
    #parse config portion of outfile
    spring_re = re.compile('^spring')
    end_re = re.compile('Time \d')
    #open temp file
    new_f = open(newfile, "w")
    with open(outfile) as f:
        for line in f:
            if spring_re.search(line):
                #split the line to grab connecting mass indexes
                line_list = line.split()
                if line_list[1] in mass_idx or line_list[2] in mass_idx:
                    new_line = line.replace(line_list[3], str(format(float(line_list[3]) + (distance_nm*10**-9), '.3g')))
                    new_f.write(new_line)
                else:
                    new_f.write(line)
            elif end_re.search(line):
                break
            else:
                new_f.write(line)
    new_f.write(coords.strip())
    new_f.write('\n')
    new_f.close()

def cen_cen_distances(outfile):
    """Returns a numpy array of the cen-cen distances in models with kmts added"""
    import numpy as np
    import chromoSnake as cs
    import re

    #going to need to grab the centromere data
    #the centromeres will be attached to the kmt beads via springs
    #we can grab the kmt ids using super_mass_indexes
    mass_idx = list(cs.super_mass_indexs(outfile))
    #loop through the config header portion of the outfile to grab spring lines
    spring_re = re.compile('^spring')
    cens = []
    with open(outfile) as f:
        for line in f:
            if spring_re.search(line):
                line_list = line.split()
                mass1 = line_list[1]
                mass2 = line_list[2]
                if mass1 in mass_idx:
                    cens.append(mass2)
                elif mass2 in mass_idx:
                    cens.append(mass1)
    #Since chromoshake flips the bead order upon printout we need to convert config indices to timeline indices
    num_masses = cs.number_of_masses(outfile)
    cens_idx = [str((num_masses - 1) - int(cen)) for cen in cens]
    #new re-loop through file and pull out the coordinates FROM THE TIME/COORD section NOT the header
    cnt = 0
    time_re = re.compile("Time ")
    time_bool = False
    cen_coords = []
    with open(outfile) as f:
        for line in f:
            if time_bool:
                if str(cnt) in cens_idx:
                    cen_coords.append(line.split())
                    cnt = cnt + 1
                else:
                    cnt = cnt + 1
            if time_re.search(line):
                time_bool = True
    #calculate the distances assuming pairings are 0 to 15, 1 to 16 etc
    dists = []
    for cen in range(int(len(cens)/2)):
        dists.append(np.linalg.norm(np.array(cen_coords[cen]).astype('float')-np.array(cen_coords[cen+15]).astype('float')))
    return(np.array(dists))

def cohesin_mass_idxs(config_file):
    """Search a spindle config file for cohesin masses assuming cohesin masses are only white/5 beads"""
    import re
    
    start_re = re.compile('^mass')
    end_re = re.compile('5$')
    cohesin_idxs = list()
    with open(config_file) as f:
        for line in f:
            if start_re.search(line) and end_re.search(line):
                split_line = line.split()
                cohesin_idxs.append(split_line[1])
    return(cohesin_idxs)
