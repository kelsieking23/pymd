class Solvent:
    def __init__(self, pdb):
        self.pdb = pdb
        self.molecules = self.getMolecules()
    
    def getMolecules(self):
        i = 0
        k = 0
        molecules = []
        f = open(self.pdb, 'r')
        contents = f.readlines()
        f.close()
        for line in contents:
            if not line.startswith('ATOM'):
                continue
            line_parts = line.split()
            atom_name = line_parts[2]
            if (atom_name == 'OW') or (atom_name == 'OM'):
                data = [line, contents[i+1], contents[i+2]]
                molecule = SolventMolecule(index=k, data=data)
                molecules.append(molecule)
                k += 1
            if (atom_name == 'NA') or (atom_name == 'CL'):
                data = [line]
                molecule = SolventMolecule(index=k, data=data)
                molecules.append(molecule)
                k += 1
            i += 1
        return molecules

class SolventMolecule:
    def __init__(self, index, data):
        self.index = index
        self.data = data
        self.atoms = self.getAtoms()
    
    def getAtoms(self):
        atoms = []
        i = 0
        for line in self.data:
            atom = SolventAtom(index=i, data=line, molecule=self)
            atoms.append(atom)
            i += 1
        return atoms

class SolventAtom:
    def __init__(self, index, data, molecule):
        self.index = index
        self.molecule = molecule
        self.id = str(self.molecule.index) + str(self.index)
        self.data = data
        self.coordinates = self.getCoordinates()
    
    def getCoordinates(self):
        return list(map(float, self.data.split()[5:8]))