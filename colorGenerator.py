import numpy as np
import re
from math import sqrt


class ChromoSim:
    """
    A class containing a reference to a chromoShake simulation
    """

    def __init__(self, filepath, SUPER_MASS=None):
        self.SUPER_MASS = SUPER_MASS
        try:
            open(filepath)
        except IOError as err:
            print(err)
        self.filepath = filepath
        self.mass_list, self.spring_list, self.hinge_list = self.parse_header()
        self.super_bool_array = self.get_super_masses()
        self.condensin_spring_indexes, self.cohesin_masses_array = self.get_smc_springs_masses()
        self.mass_labels = self.gen_mass_labels()
        self.sim_dict = self.parse_timepoints()

    def parse_header(self):
        # set of toggles
        in_meta = True
        in_struct = False
        in_color = False
        in_output = False
        # set of lists
        mass_list = list()
        spring_list = list()
        hinge_list = list()
        # compile re patterns
        meta_pattern = re.compile('^meta ')
        color_pattern = re.compile('^MassColors')
        struct_pattern = re.compile('^structure')
        out_pattern = re.compile('^Time ')
        mass_pattern = re.compile('^mass ')
        spring_pattern = re.compile('^spring ')
        hinge_pattern = re.compile('^hinge ')
        with open(self.filepath) as f:
            for line in f:
                line = line.strip()
                meta_match = meta_pattern.match(line)
                color_match = color_pattern.match(line)
                struct_match = struct_pattern.match(line)
                out_match = out_pattern.match(line)
                if in_color and out_match is not None:
                    in_output = True
                    in_color = False
                if in_struct and color_match is not None:
                    in_color = True
                    in_struct = False
                if in_meta and struct_match is not None:
                    in_struct = True
                    in_meta = False
                if in_meta and meta_match is None:
                    raise Exception('File does not start with metadata section. Is this a chromoShake outfile?')
                if in_struct:
                    if mass_pattern.match(line) is not None:
                        mass_list.append(line.split())
                    elif spring_pattern.match(line) is not None:
                        spring_list.append(line.split())
                    elif hinge_pattern.match(line) is not None:
                        hinge_list.append(line.split())
                if in_output:
                    f.close()
                    break
        return mass_list, spring_list, hinge_list

    def get_super_masses(self):
        """
        Identify super masses
        params:
            super_mass : A numpy double constant
        returns: A numpy array of boolean values to index mass_list
        """

        mass_array = np.array(self.mass_list)
        mass_masses = mass_array[:, 2].astype(np.double)
        if self.SUPER_MASS is not None and self.SUPER_MASS.dtype.name == 'float64':
            super_masses = mass_masses == self.SUPER_MASS
        else:
            super_bool_array = mass_masses == mass_masses.max()

        return super_bool_array

    def get_smc_springs_masses(self):
        """
        Identify condensin springs
        params: None
        returns: A numpy array of boolean values to index spring_list
        """

        spring_array = np.array(self.spring_list)
        start_idxs = spring_array[:, 1].astype(np.int)
        end_idxs = spring_array[:, 2].astype(np.int)
        diff_idxs = end_idxs - start_idxs
        condensin_spring_indexes = np.logical_and(diff_idxs != 1, diff_idxs != -15)
        # cohesin springs will have -15 in diff_idxs
        mass_array = np.array(self.mass_list)
        cohesin_spring_indexes = diff_idxs == -15
        cohesin_springs_start_end_idx = spring_array[cohesin_spring_indexes, 1:3].astype(np.int)
        cohesin_masses_array = np.zeros((16, mass_array.shape[1], cohesin_springs_start_end_idx.shape[0]),
                                        dtype=mass_array.dtype)
        for i, (end, start) in enumerate(cohesin_springs_start_end_idx):
            cohesin_ring = mass_array[start:end + 1, :]
            cohesin_masses_array[:, :, i] = cohesin_ring

        return condensin_spring_indexes, cohesin_masses_array

    def gen_mass_labels(self):
        """
        Create a list of labels for each mass in a simulation. Masses can be "super", "chrX_top", "chrX_bottom", or
        cohesin_X where X is a chromatid or complex index number
        params: self
        returns: A list of labels for each mass in a simulation
        """

        # super massive masses
        super_mass_indexes = np.array(np.where(self.super_bool_array))
        super_mass_indexes = super_mass_indexes.reshape([-1, 2])
        strand_length_array = super_mass_indexes[:, 1] - super_mass_indexes[:, 0] + 1
        if len(np.unique(strand_length_array)) == 1:
            strand_length = np.unique(strand_length_array)
        else:
            raise Exception("Calculated pericentromere strand lengths are not consistent.")

        if strand_length % 2 == 0:
            mid_mass_indexes = super_mass_indexes[:, 0] + strand_length / 2
        else:
            mid_mass_indexes = super_mass_indexes[:, 0] + np.floor(strand_length / 2)

        # pair strands for C-loop indexing
        top_list = []
        bottom_list = []
        for pair in super_mass_indexes:
            mid_bool = np.logical_and(mid_mass_indexes > pair[0], mid_mass_indexes < pair[1])
            top_half_indexes = [top for top in range(pair[0], mid_mass_indexes[mid_bool].astype(np.int)[0])]
            bottom_half_indexes = [bottom for bottom in
                                   range(mid_mass_indexes[mid_bool].astype(np.int)[0], pair[1] + 1)]
            top_list.append(top_half_indexes)
            bottom_list.append(bottom_half_indexes)
        top_array = np.array(top_list)
        bottom_array = np.array(bottom_list)

        # ensure that the first column of super_mass_indexes all have same z-value
        mass_array = np.array(self.mass_list)
        if len(np.unique(mass_array[super_mass_indexes[:, 0], 5])) != 1:
            raise Exception("Super massive masses are not indexed as expected.")

        # pair strands of super massive beads to correctly label C-loops
        super_mass_ends = mass_array[super_mass_indexes[:, 0], 3:5].astype(np.float)
        c_loop_pair_list = []
        for i, index in enumerate(super_mass_indexes[:, 0]):
            ref_xy = mass_array[index, 3:5].astype(np.float)
            dist_array = np.sum((super_mass_ends - ref_xy) ** 2, axis=1)
            dist_array[0 == dist_array] = np.nan
            c_loop_pair_list.append([i, np.nanargmin(dist_array)])

        c_loop_indexes = []
        for index in range(16):
            top_c_loop = np.concatenate((top_array[c_loop_pair_list[index][0], :],
                                         top_array[c_loop_pair_list[index][1]]))
            bottom_c_loop = np.concatenate((bottom_array[c_loop_pair_list[index][0], :],
                                            bottom_array[c_loop_pair_list[index][1]]))
            c_loop_indexes.append(top_c_loop.tolist())
            c_loop_indexes.append(bottom_c_loop.tolist())

        # cohesin masses
        cohesin_mass_indexes = self.cohesin_masses_array[:, 1, :].flatten().astype(np.int)
        cohesin_bool_list = []
        for mass in self.mass_list:
            if int(mass[1]) in cohesin_mass_indexes:
                cohesin_bool_list.append(True)
            else:
                cohesin_bool_list.append(False)

        cohesin_bool_array = np.array(cohesin_bool_list)
        # condensin-bound masses
        spring_array = np.array(self.spring_list)
        condensin_spring_array = spring_array[np.where(self.condensin_spring_indexes)]
        condensin_mass_indexes = condensin_spring_array[:, 1:3].astype(np.int)
        # create strand, cohesin, condensin-bound, and super mass labels

        mass_labels = []
        coh_cnt = 0
        coh_idx = 0
        for i, _ in enumerate(self.mass_list):
            if i in super_mass_indexes:
                mass_labels.append('super')
            elif i in cohesin_mass_indexes:
                mass_labels.append('cohesin_' + str(coh_idx))
                coh_cnt += 1
                if coh_cnt == 16:
                    coh_cnt = 0
                    coh_idx += 1

            else:
                for loop_index, index_list in enumerate(c_loop_indexes):
                    if i in index_list:
                        if loop_index % 2 == 0:
                            position = 'top'
                            chromosome = loop_index / 2 + 1
                        else:
                            position = 'bottom'
                            chromosome = (loop_index - 1) / 2 + 1
                mass_labels.append('chr' + str(int(chromosome)) + '_' + position)

        return mass_labels

    def parse_timepoints(self):
        """
        Generate a simulation dictionary of time-points and mass positions
        params : self
        outputs : sim_dict : A dictionary of dictionaries containing mass positions for each timepoint
        """

        # regexp patterns
        timeline_pattern = re.compile('^Time ')
        # toggles
        in_timeline = False
        # instantiate sim_dict
        sim_dict = dict()
        # get number of masses (must parse header prior to this)
        num_masses = len(self.mass_list)
        with open(self.filepath) as f:
            for line in f:
                line = line.strip()
                timeline_match = timeline_pattern.match(line)
                if timeline_match and not in_timeline:
                    # this is the first time line
                    in_timeline = True

                if timeline_match and in_timeline:
                    # store the time
                    time = line.strip().split()[-1]
                    mass_cnt = 0
                    sim_dict[time] = dict()

                if not timeline_match and in_timeline and line:
                    mass_idx = num_masses - 1 - mass_cnt  # must flip cnt to get index to match header
                    x_string, y_string, z_string = line.strip().split()
                    sim_dict[time][mass_cnt] = {
                        'x': x_string,
                        'y': y_string,
                        'z': z_string,
                        'label': self.mass_labels[mass_idx]
                    }
                    mass_cnt += 1
        f.close()
        return sim_dict

    def calc_tension(self):
        """
        Calculates the tension on each mass per timepoint and appends tension to sim_dict
        params: self
        returns : updated sim_dict
        """

        def calc_paired_distance(time_dict, mass_idx1, mass_idx2):
            """
            Parse the time_dict with mass_idx1 and mass_idxs vars to calculate distance between 2 masses
            """
            x1 = time_dict[mass_idx1]['x']
            x2 = time_dict[mass_idx2]['x']
            y1 = time_dict[mass_idx1]['y']
            y2 = time_dict[mass_idx2]['y']
            z1 = time_dict[mass_idx1]['z']
            z2 = time_dict[mass_idx2]['z']
            return sqrt(
                (float(x1) - float(x2)) ** 2 +
                (float(y1) - float(y2)) ** 2 +
                (float(z1) - float(z2)) ** 2
            )

        for time_key in self.sim_dict:
            my_time_dict = self.sim_dict[time_key]
            idx_tension_list = []
            for line_list in self.spring_list:
                idx1 = int(line_list[1])
                idx2 = int(line_list[2])
                distance = calc_paired_distance(my_time_dict, idx1, idx2)
                rest_length = float(line_list[3])
                spring_const = float(line_list[4])
                tension = (distance - rest_length) * spring_const
                idx_tension_list.append([idx1, idx2, tension])
            idx_tension_array = np.array(idx_tension_list)
            for mass_key in my_time_dict:
                mass_array_1 = idx_tension_array[idx_tension_array[:, 0] == mass_key, :]
                mass_array_2 = idx_tension_array[idx_tension_array[:, 1] == mass_key, :]
                if mass_array_1.size > 0 and mass_array_2.size > 0:
                    mass_tension = (mass_array_1[0, -1] + mass_array_2[0, -1]) / 2.0
                elif mass_array_1.size == 0 and mass_array_2.size > 0:
                    mass_tension = mass_array_2[0, -1]
                elif mass_array_1.size > 0 and mass_array_2.size == 0:
                    mass_tension = mass_array_1[0, -1]
                else:
                    raise Exception("Could not find mass {0} in idx_tension_array".format(mass_key))
                self.sim_dict[time_key][mass_key]['tension'] = mass_tension
        return self
