class ChromoSim:
    """
    A class containing a reference to a chromoShake simulation
    """

    def __init__(self, filepath, SUPER_MASS=None):
        try:
            open(filepath)
        except IOError as err:
            print(err)
        self.filepath = filepath
        self.mass_list, self.spring_list, self.hinge_list = self.parse_header()
        self.super_mass_indexes = self.get_super_masses()
        self.condensin_spring_indexes = self.get_condensin_springs()

    def parse_header(self):
        import re
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

    def get_super_masses(self, SUPER_MASS=None):
        """
        Identify super masses
        params:
            super_mass : A numpy double constant
        returns: A numpy array of boolean values to index mass_list
        """
        import numpy as np

        mass_array = np.array(self.mass_list)
        mass_masses = mass_array[:, 2].astype(np.double)
        if SUPER_MASS is not None and SUPER_MASS.dtype.name == 'float64':
            super_masses = mass_masses == SUPER_MASS
        else:
            super_mass_indexes = mass_masses == mass_masses.max()

        return super_mass_indexes

    def get_condensin_springs(self):
        """
        Identify condensin springs
        params: None
        returns: A numpy array of boolean values to index spring_list
        """
        import numpy as np

        spring_array = np.array(self.spring_list)
        start_idxs = spring_array[:, 1].astype(np.int)
        end_idxs = spring_array[:, 2].astype(np.int)
        diff_idxs = end_idxs - start_idxs
        condensin_spring_indexes = diff_idxs != 1
        return condensin_spring_indexes
