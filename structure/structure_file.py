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

    def get_atom_data_pdb(self, line):
        pass
    
    def crd(self, writer=False):
        atom_index = 0
        residue_index = -1
        last_residue = None
        last_chain = None
        chain_index = -1
        chain = 'A'
        charge = ''
        box = (0,0,0)
        model = 1
        for line in self.structureReader():
            if (line.startswith('*')) or ('EXT' in line) or (line.strip() == ''):
                continue
            else:
                line_parts = line.strip().split()
                atom_number = line_parts[0]
                atom_name = line_parts[3]
                x = float(line_parts[4])*10
                y = float(line_parts[5])*10
                z = float(line_parts[6])*10
                segid = line_parts[7]
                residue_name = line_parts[2]
                residue_number = int(line_parts[1])
                residue_id = residue_name + residue_name + str(residue_number)
                if residue_id != last_residue:
                        residue_index += 1
                elem = atom_name[0]
                charge = ''
                if segid != last_chain:
                    chain_index += 1
                    if chain != 'A':
                        chain = chr(ord(chain) + 1)
                temp = 0.00
                occ = 0.00
                if not writer:
                    atom = AtomData(atom_number, atom_index, atom_name, residue_name, residue_id, chain, chain_index, residue_number, residue_index, x, y, z, occ, temp, segid, elem, charge, model, box)
                    yield atom
                else:
                    atom = ['ATOM', str(atom_number), atom_name, residue_name, chain, str(residue_number), x, y, z, '1.00', '0.00', segid, elem]
                    yield atom

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
                if (chain == '') or (chain == ' '):
                    _chain_index = -1
                else:
                    _chain_index = chain_index
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occ = float(line[54:60].strip())
                temp = float(line[60:66].strip())
                segid = line[72:76].strip()
                elem = line[76:78].strip()
                _charge = line.strip()[-1]
                if (_charge == '+') or (_charge == '-'):
                    charge = _charge
                else:
                    charge = ''
                atom = AtomData(atom_number, atom_index, atom_name, residue_name, residue_id, 
                                chain, _chain_index, residue_number, residue_index,
                                x, y, z, occ, temp, segid, elem, charge, model, box)
                atom_index += 1
                last_chain = chain
                last_residue = residue_name + str(residue_number)
                yield atom
            elif line.startswith('HETATM'):
                #TODO: add hetatm
                pass
            else:
                yield line
        return atoms

    def gro(self):
        atoms = []
        atom_index = 0
        residue_index = 0
        i = 0
        for line in self.structureReader():
            if i < 2:
                i += 1
                continue
                # yield line
            try:
                residue_number = int(line[0:5].strip())
                residue_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                atom_number = int(line[15:20].strip())
                residue_id = '{}{}'.format(residue_name, residue_number)
                x = float(line[20:29].strip())
                y = float(line[29:37].strip())
                z = float(line[37:46].strip())
                if (len(line) > 45):
                    vx = float(line[46:55].strip())
                    vy = float(line[55:64].strip())
                    vz = float(line[64:73].strip())
                else:
                    vx = 0.0
                    vy = 0.0
                    vz = 0.0
                atom = GroAtomData(atom_number, atom_index, atom_name, residue_name, residue_id,
                                   residue_number, residue_index, x, y, z, vx, vy, vz)
                atom_index += 1
                residue_index += 1
                yield atom
            except:
                pass
            i += 1
        return atoms

    def sdf(self):
        return []
    
    def mol(self):
        return []

                
@dataclass
class GroAtomData:
    atom_number: int
    atom_index: int
    atom_name: str
    residue_name: str
    residue_id: str
    residue_number: int
    residue_index: int
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float


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
    occ:float
    temp:float
    segid: str
    elem: str
    charge:str
    model: float
    box: tuple
