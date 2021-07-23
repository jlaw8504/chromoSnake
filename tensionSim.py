import re
from math import sqrt
from colorGenerator import ChromoSim


class TensionSim:
    """
    A class for exploring tension in chromoShake simulations
    """

    def __init__(self, filepath, supermass=None):
        self.supermass = supermass
        try:
            open(filepath)
        except IOError as err:
            print(err)
        self.filepath = filepath

    def parse_sim(self, outfile):
        """
        Parse an entire chromoShake outfile in one pass
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

        # set of toggles
        in_meta = True
        in_struct = False
        in_color = False
        in_timeline = False
        # compile re patterns
        meta_pattern = re.compile('^meta ')
        color_pattern = re.compile('^MassColors')
        struct_pattern = re.compile('^structure')
        timeline_pattern = re.compile('^Time ')
        # use ChromoSim class to parse header
        my_sim_header = ChromoSim(self.filepath)
        # begin parsing the file
        with open(self.filepath) as f:
            with open(outfile, 'w') as f_out:
                for line in f:
                    line = line.strip()
                    meta_match = meta_pattern.match(line)
                    color_match = color_pattern.match(line)
                    struct_match = struct_pattern.match(line)
                    timeline_match = timeline_pattern.match(line)
                    if in_color and timeline_match is not None:
                        in_timeline = True
                        num_masses = len(my_sim_header.mass_list)
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
                        # this is now taken care of by the my_sim_header object
                        continue
                    if timeline_match and in_timeline:
                        # store the time
                        time = line.strip().split()[-1]
                        # reset the mass counter
                        mass_cnt = 0
                        # reset the dictionary
                        time_dict = dict()

                    if not timeline_match and in_timeline and line:
                        mass_idx = str(num_masses - 1 - mass_cnt)  # must flip cnt to get index to match header
                        x_string, y_string, z_string = line.strip().split()
                        time_dict[mass_idx] = {
                            'x': x_string,
                            'y': y_string,
                            'z': z_string,
                        }
                        mass_cnt += 1
                        # once all masses are parsed, calculate tension
                        if mass_cnt == num_masses:
                            print(f'time: {time}')
                            for spring_list_info in my_sim_header.spring_list:
                                spring_mass_idx_1 = spring_list_info[1]
                                spring_mass_idx_2 = spring_list_info[2]
                                rest_length = float(spring_list_info[3])
                                spring_k = float(spring_list_info[4])
                                dist = calc_paired_distance(time_dict, spring_mass_idx_1, spring_mass_idx_2)
                                delta_dist = dist - rest_length
                                tension = delta_dist * spring_k
                                f_out.write(
                                    '{0}, {1}, {2}, {3}, {4}, {5}, {6}\n'.format(
                                        time,
                                        spring_mass_idx_1,
                                        spring_mass_idx_2,
                                        tension,
                                        my_sim_header.mass_labels[int(spring_mass_idx_1)],
                                        my_sim_header.mass_labels[int(spring_mass_idx_2)],
                                        my_sim_header.loop_labels[int(spring_mass_idx_1)],
                                        my_sim_header.loop_labels[int(spring_mass_idx_2)]
                                    )
                                )
