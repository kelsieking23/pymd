class StaticLigand:

    def __init__(self, structure, name=None, covalent_id=None, model=None, ignh=False, alt=None):
        '''
        Ligand class holds data related to ligand structure. Only holds data related to one model. This model can be specified. 
        Arguments:
        ** structure: string. Path to structure file. Currently, must be a .pdbqt file. 
        ** name (optional): string. Name to be associated with this ligand. The name dictates how it will be indexed in output .csv files and internally. If no argument is
           passed, will default to the file name (NOT file path).
        ** covalent_id (optional): string. If the ligand is covalent, or otherwise in the protein structure file, this indicates which lines of data are associated with the 
           ligand in the PDB file. ID syntax is <residue_name><residue_number> (ex: FHR405 (see PDBID: 6LZE for reference))
        ** model (optional): string. If ligand structure file has multiple models, specify which to use as coordinates. Defaults to MODEL 1. 
        ** alt (optional):  string. Gives ligand an alternative name, which can be used in certain functions. 
        '''
        self.atoms = []
        self.structure = structure
        self.covalent_id = covalent_id
        if name is not None:
            self.name = name
        else:
            try:
                self.name = os.path.basename(self.structure)[:-4]
            except:
                self.name = self.structure

        if self.covalent_id is not None:
            self.coordinates = self.getCovalentCoordinates(ignh=ignh)
        else:
            self.coordinates = self.getLigandCoordinates(model=model, ignh=ignh)
        
        # self.log = self.getLog()
        if alt is not None:
            self.alt = alt
        
        self._minimum_distance = None
        self.atomic_masses = {
            'C':12.011,
            'A':12.011,
            'N':14.007,
            'H':1.008,
            'O':15.999,
            'S':32.06
        }
    
        
    def getLog(self):
        directory = os.path.dirname(self.structure)
        files = os.listdir(directory)
        for filename in files:
            if 'log.txt' in filename:
                log_path = os.path.join(directory, filename)
                return log_path
        return 0

    def getCOM(self):
        '''
        Gets residue center of mass for all residues. 
        ** Returns: dict, where key = residue_id, and value = xyz coordinates of COM (tuple)
        '''
        pass
        x = []
        y = []
        z = []
        i = 0
        ligand_mass = 0
        for atom in self.atoms:
            mass = self.atomic_masses[atom.type]
            ligand_mass += mass
            x.append(atom.coordinates[0] * mass)
            y.append(atom.coordinates[1] * mass)
            z.append(atom.coordinates[2]* mass)
            i += 1
        com_x = sum(x) / ligand_mass
        com_y = sum(y) / ligand_mass
        com_z = sum(z) / ligand_mass
        com = (com_x, com_y, com_z)
        return com

    @property
    def log(self):
        if self.getLog() == 0:
            return None
        else:
            return self.getLog()

    @property
    def minimum_distance(self):
        return self._minimum_distance

    def getCovalentCoordinates(self, ignh):
        '''
        Gets coordinates from protein structure file. 
        Returns: list
        '''
        coordinates = []
        f = open(self.structure, 'r')
        for line in f:
            line_parts = line.split()
            if 'HETATM' in line_parts[0]:
                atom_type = line_parts[-1][0]
                if ignh == False:
                    if atom_type == 'H':
                        continue
                residue_name = line_parts[3]
                residue_number = line_parts[5]
                residue_id = residue_name + residue_number
                if residue_id == self.covalent_id:
                    atom_coordinates = list(map(float, line_parts[6:9]))
                    coordinates.append(atom_coordinates)
        f.close()
        return coordinates
    
    def getLigandCoordinates(self, model, ignh):
        '''
        Takes coordinate data from self.split() and returns appropriate coordinate data, depending on if the structure file passed has multiple models.
        Arguments:
        ** model (optional): string. If structure file has multiple models, specifies which model to take coordinates of. 
        Returns: list
        '''
        coordinates = self.split(ignh)
        if model is None:
            model = 'MODEL 1'
        if isinstance(coordinates, list):
            return coordinates
        if isinstance(coordinates, dict):
            return coordinates[model]

    def split(self, ignh):
        '''
        Splits structure file into multiple models.
        Returns: list or dict (depending on if the structure file has multiple models or not)
        '''
        model_coordinates = {}
        ligand_coordinates = []
        current_model = None
        f = open(self.structure, 'r')
        for line in f:
            line_parts = line.split()
            if 'HETATM' in line_parts[0]:
                if ignh == True:
                    atom_type = line_parts[-1][0]
                    if atom_type == 'H':
                        continue
                if line_parts[4].isalpha():
                    coords = list(map(float, line_parts[6:9]))
                else:
                    coords = list(map(float, line_parts[5:8]))                    
                ligand_coordinates.append(coords)
            if 'MODEL' in line:
                current_model = line.strip()
                if current_model not in model_coordinates.keys():
                    model_coordinates[current_model] = []
            if 'ENDMDL' in line:
                model_coordinates[current_model] = ligand_coordinates
                ligand_coordinates = []
        f.close()
        if model_coordinates != {}:
            return model_coordinates
        else:
            return ligand_coordinates

    def atomCount(self, heavy=True):
        '''
        Gets atom count for ligand. 
        Arguments:
        ** heavy (optional): boolean. Specifies if atom count should include heavy atoms only (all atoms except H) or include H. 
           If true, only counts heavy atoms. If false, counts all atoms. 
        '''
        if heavy == False:
            return(len(self.coordinates))
        if heavy == True:
            atoms = 0
            if self.covalent_id is not None:
                f = open(self.structure, 'r')
                for line in f:
                    line_parts = line.split()
                    if 'HETATM' in line_parts[0]:
                        residue_name = line_parts[3]
                        residue_number = line_parts[5]
                        residue_id = residue_name + residue_number
                        if residue_id == self.covalent_id:
                            atom_type = line_parts[-1][0]
                            if atom_type == 'H':
                                continue
                            else:
                                atoms += 1
                    if 'ENDMDL' in line_parts:
                        break
                f.close()
            else:
                f = open(self.structure, 'r')
                for line in f:
                    line_parts = line.split()
                    if 'HETATM' in line_parts[0]:
                        atom_type = line_parts[-1][0]
                        if atom_type == 'H':
                            continue
                        else:
                            atoms += 1
                    if 'ENDMDL' in line_parts:
                        break
                f.close()
        return atoms