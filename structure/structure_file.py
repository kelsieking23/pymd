import os
from dataclasses import dataclass

class StructureFile:

    def __init__(self, inp):
        self.inp = inp

    def structureReader(self):
        '''
        Iterates through structure file and yeilds lines.
        '''
        f = open(self.inp, 'r')
        contents = f.readlines()
        f.close()
        for line in contents:
            yield line
    
    def pdb(self):
        atoms = []
        model = 1
        last_chain = None
        last_residue = None
        atom_index = 0
        chain_index = -1
        residue_index = -1
        box = (0,0,0)
        for line in self.structureReader():
            line_parts = line.split()
            if len(line_parts) == 0:
                continue
            if line.startswith('ENDMDL'):
                model = model + 1
                atom_index = 0
                chain_index = 0
                residue_index = 0
                last_chain = None
                last_residue = None
            elif line.startswith('ATOM'):
                atom_number = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:21].strip()
                chain = line[21]
                residue_number = int(line[22:26].strip())
                residue_id = residue_name + str(residue_number)
                if (residue_name + str(residue_number) != last_residue):
                    residue_index += 1
                if (last_chain != chain):
                    chain_index += 1
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                segid = line[72:76].strip()
                elem = line[76:78].strip()
                atom = AtomData(atom_number, atom_index, atom_name, residue_name, residue_id, 
                                chain, chain_index, residue_number, residue_index,
                                x, y, z, segid, elem, model, box)
                atom_index += 1
                last_chain = chain
                last_residue = residue_name + str(residue_number)
                yield atom
            elif line.startswith('HETATM'):
                #TODO: add hetatm
                pass
            else:
                continue
        return atoms

    def gro(self):
        return []

    def sdf(self):
        return []
    
    def mol(self):
        return []

                

@dataclass
class AtomData:

    atom_number: int
    atom_index: int
    atom_name: str
    residue_name: str
    residue_id: str
    chain: str
    chain_index: int
    residue_number: int
    residue_index: int
    x: float
    y: float
    z: float
    segid: str
    elem: str
    model: float
    box: tuple


