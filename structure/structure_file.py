import os
from dataclasses import dataclass
from typing import Union
# from pymd.structure.library import _canonical

class StructureFile:

    def __init__(self, inp, atom_data=[]):
        self.inp = inp
        self.atoms = atom_data
        if self.inp is not None:
            _, ext = os.path.splitext(self.inp)
            self.ext = ext[1:]
            if self.ext != 'gro':
                self.atoms = [atom for atom in self.read()]
        else:
            self.ext = None
        # self.atoms = [atom for atom in self.read()]

    def _structureReader(self):
        with open(self.inp, 'r') as f:
            for line in f:
                yield line

    def _atomDataIterator(self):
        for atom in self.atoms:
            yield atom
        return self.atoms

    def read(self):
        '''
        Iterates through structure file and yeilds lines.
        '''
        if self.ext is None:
            return self._atomDataIterator()
        if self.ext == 'pdb':
            return self.pdb()


    def write(self, out):
        if self.ext == 'pdb':
            self.write_pdb(self.atoms, out)
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
        for line in self._structureReader():
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
                    atom = AtomData(atom_number, atom_index, atom_name, residue_name, residue_id, chain, chain_index, 
                                    residue_number, residue_index, x, y, z, occ, temp, segid, elem, charge, model, box, 
                                    line, False, '', 'crd')
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
        for line in self._structureReader():
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
            elif (line.startswith('ATOM')) or (line.startswith('HETATM')):
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
                                x, y, z, occ, temp, segid, elem, charge, model, box, line, True, line[0:6].strip(), 'pdb', self)
                atom_index += 1
                last_chain = chain
                last_residue = residue_name + str(residue_number)
                yield atom
            else:
                continue
        return atoms

    def gro(self):
        atoms = []
        atom_index = 0
        residue_index = 0
        i = 0
        for line in self._structureReader():
            if i < 2:
                i += 1
                continue
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
                    vx = float(line[46:53].strip())
                    vy = float(line[53:61].strip())
                    vz = float(line[61:].strip())
                else:
                    vx = 0.0
                    vy = 0.0
                    vz = 0.0
                atom = GroAtomData(atom_number, atom_index, atom_name, residue_name, residue_id,
                                   residue_number, residue_index, x, y, z, vx, vy, vz, line)
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

    @staticmethod
    def write_gro(atom_data, out, title='Title', box=[0.0, 0.0, 0.0]):
        with open(out, 'w') as f:
            f.write(f'{title}\n')
            f.write(f' {len(atom_data)}\n')
            for atom in atom_data:
                ld = [atom.residue_number, atom.residue_name, atom.atom_name, atom.atom_number, atom.x, atom.y, atom.z,
                    atom.vx, atom.vy, atom.vz]
                line = '{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.4f}{:>8.4f}{:>8.4f}\n'.format(*ld)
                f.write(line)
            f.write('   {:.7f}   {:.7f}   {:.7f}\n'.format(*box))
    
    @staticmethod
    def write_pdb(atom_data, out):
        with open(out, 'w') as f:
            for atom in atom_data:
                f.write(atom.line)
            
    @classmethod
    def fromAtomData(cls, atom_data):
        return cls(inp=None, atom_data=atom_data)
    
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
    line: str

    def update_line(self):
        ld = [self.residue_number, self.residue_name, self.atom_name, self.atom_number, self.x, self.y, self.z,
              self.vx, self.vy, self.vz]
        line = '{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.4f}{:>8.4f}{:>8.4f}\n'.format(*ld)
        self.line = line
        return self


class AtomData:

    def __init__(self, atom_number: int, atom_index: int, atom_name: str, residue_name: str, residue_id: str,
                 chain: str, chain_index: int, residue_number: int, residue_index: int, x: float, y: float, z: float,
                 occ:float, temp: float, segid: str, elem: str, charge: str, model: float, box: tuple, line: str, 
                 is_pdb: bool, pdb_label: str, ext: str, parent: Union[StructureFile, None]):
        self._atom_number = atom_number
        self.atom_index = atom_index
        self.atom_name = atom_name
        self.residue_name = residue_name
        self.residue_id = residue_id
        self.chain = chain
        self.chain_index = chain_index
        self._residue_number = residue_number
        self.residue_index = residue_index
        self.x = x
        self.y = y
        self.z = z
        self.occ = occ
        self.temp = temp
        self.segid = segid
        self.elem = elem
        self.charge = charge
        self.model = model
        self.box = box
        self._line = line
        self.is_pdb = is_pdb
        self.pdb_label = pdb_label
        self.ext = ext
        self.parent = parent
        self.update_dict: dict={
                'pdb':self._update_pdb,
                'crd':self._update_crd}
        self._cannonical = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'GLY', 'PRO', 
                            'SER', 'THR', 'ASN', 'GLN', 'CYS', 'HIS', 'HSD', 'ASP', 'GLU', 'ARG', 'LYS']
        
    def __str__(self):
        return f'<pymd.structure.structure_file.AtomData Object>: {self.atom_number} {self.atom_name} {self.residue_number} {self.residue_name} {self.chain}'

    def _update_parent(self):
        self._line = self.update_dict[self.ext]
        self.parent.atoms[self.atom_index] = self

    @property
    def atom_number(self):
        return self._atom_number
    
    @atom_number.setter
    def atom_number(self, number):
        try:
            self._atom_number = int(number)
        except Exception:
            raise ValueError(f'Error setting attribute "atom_number"')
        self._update_parent()
    
    @property
    def residue_number(self):
        return self._residue_number
    
    @residue_number.setter
    def residue_number(self, number):
        try:
            self._residue_number = int(number)
        except Exception:
            raise ValueError(f'Error setting attribute "residue_number"')
        self._update_parent()

    @property
    def line(self):
        if self.ext == 'pdb':
            self._line = self._update_pdb()
            if isinstance(self.parent, StructureFile):
                self.parent.atoms[self.atom_index] = self
        return self._line
    
    def _update_pdb(self):
        line = [self.pdb_label, self.atom_number, self.atom_name, self.residue_name, self.chain, self.residue_number, self.x, self.y, self.z,
                self.occ, self.temp, self.segid, self.elem]
        string = "{:6s}{:5d} {:^4s} {:^4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>10s}  {:<3s}\n".format(*line)
        return string

    def _update_crd(self):
        return ''
    
    def _in_rtp(self, rtp):
        with open(rtp, 'r') as f:
            for line in f:
                if line.startswith('[ '):
                    entry = line.strip()[2:-2]
                    if entry.upper() == self.residue_name.upper():
                        return True
        return False
    
    def mol_type(self, ff='guess'):
        '''
        at this moment super rudimentary, but just want to be able to ID if it is protein or not protein. 
        ff (str) : can be 'guess', which means just go off of _cannonical() which is hard coded
                   ff can also be a path to a force field folder that contains an aminoacids.rtp file
        '''
        if ff == 'guess':
            if self.residue_name in self._cannonical:
                return 'protein'
            else:
                return 'other'
        else:
            if self._in_rtp(os.path.join(ff, 'aminoacids.rtp')):
                return 'protein'
            else:
                for file in os.listdir(ff):
                    if file.endswith('rtp'):
                        if self._in_rtp(os.path.join(ff, file)):
                            return os.path.splitext(file)[0]
                return 'other'
