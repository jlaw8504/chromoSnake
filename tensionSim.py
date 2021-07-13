import re
from math import sqrt
import numpy as np

class TensionSim():
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
        self.sim_dict = self.parse_sim()

    def parse_sim(self):
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
        mass_pattern = re.compile('^mass ')
        spring_pattern = re.compile('^spring ')
        hinge_pattern = re.compile('^hinge ')
        # instantiate nested dictionary
        sim_dict = dict()
        sim_dict['struct'] = dict()
        # begin parsing the file
        with open(self.filepath) as f:
            for line in f:
                line = line.strip()
                meta_match = meta_pattern.match(line)
                color_match = color_pattern.match(line)
                struct_match = struct_pattern.match(line)
                timeline_match = timeline_pattern.match(line)
                if in_color and timeline_match is not None:
                    in_timeline = True
                    num_masses = len(sim_dict['struct'])
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
                        mass_list = line.split()
                        mass_idx = mass_list[1]
                        sim_dict['struct'][mass_idx] = {
                            'mass_kg': float(mass_list[2]),
                            'spring_idx_list': [],
                            'hinge_idx_list': []
                        }

                    elif spring_pattern.match(line) is not None:
                        spring_list = line.split()
                        sim_dict['struct'][spring_list[1]]['spring_idx_list'].append(
                            [spring_list[1], spring_list[2], spring_list[3], spring_list[4]]
                        )
                        sim_dict['struct'][spring_list[2]]['spring_idx_list'].append(
                            [spring_list[1], spring_list[2], spring_list[3], spring_list[4]]
                        )

                    elif hinge_pattern.match(line) is not None:
                        hinge_list = line.split()
                        sim_dict['struct'][hinge_list[2]]['hinge_idx_list'].append(
                            [hinge_list[1], hinge_list[2], hinge_list[3]]
                        )

                if timeline_match and in_timeline:
                    # store the time
                    time = line.strip().split()[-1]
                    mass_cnt = 0
                    parsed_mass_list = []
                    sim_dict[time] = dict()

                if not timeline_match and in_timeline and line:
                    mass_idx = str(num_masses - 1 - mass_cnt)  # must flip cnt to get index to match header
                    x_string, y_string, z_string = line.strip().split()
                    sim_dict[time][mass_idx] = {
                        'x': x_string,
                        'y': y_string,
                        'z': z_string,
                    }
                    mass_cnt += 1
                    # once all masses are parsed, calculate tension
                    if mass_cnt == num_masses:
                        for key in sim_dict['struct']:
                            tension_list = []
                            for spring_list_info in sim_dict['struct'][key]['spring_idx_list']:
                                dist = calc_paired_distance(sim_dict[time], spring_list_info[0], spring_list_info[1])
                                rest_length = float(spring_list_info[2])
                                delta_dist = dist - rest_length
                                tension_list.append(delta_dist * float(spring_list_info[3]))
                            tension_array = np.array(tension_list)
                            sim_dict[time][key]['tension'] = tension_array.mean()

            return sim_dict
