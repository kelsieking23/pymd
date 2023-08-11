import os
import numpy as np
import pandas as pd
import sys

from pymd.structure.structure_file import StructureFile
from pymd.utilities import library
from pymd.utilities.rewritepdb import addChainID, writePDB, editChainIDResidue, fixBadCoordinates


class Protein:

    def __init__(self, structure, ignore=[]):
        self.structure = os.path.abspath(structure)
        # atomic mass map for center of mass calculations
        self.atomic_masses = {
            'C':12.011,
            'A':12.011,
            'N':14.007,
            'H':1.008,
            'O':15.999,
            'S':32.06,
            'P':30.974
        }
        self._residues = self.loadTopology()
        self.loadModel(1)
        self.ignore = ignore


    def structureReader(self):
        '''
        Iterates through structure file and yeilds lines.
        '''
        if isinstance(self.structure, str):
            f = open(self.structure, 'r')
            contents = f.readlines()
            f.close()
            for line in contents:
                yield line
        else:
            for line in self.structure:
                yield line

    def loadTopology(self):
        struct_data = StructureFile(self.structure)
        if self.structure.endswith('pdb'):
            func = struct_data.pdb
        elif self.structure.endswith('gro'):
            func = struct_data.gro
        elif self.structure.endswith('sdf'):
            func = struct_data.sdf
        elif (self.structure.endswith('mol')) or (self.structure.endswith('mol2')):
            func = struct_data.mol
        else:
            raise ValueError('Only pdb, gro, sdf, mol, or mol2 are accepted')
        res_data = []
        residues = []
        atoms = []
        last_atom = None
        for atom in func():
            if last_atom is None:
                res_data.append(atom)
            else:
                if last_atom.residue_id != atom.residue_id:
                    res = Residue(res_data)
                    residues.append(res)
                    res_data = []
                    res_data.append(atom)
                else:
                    res_data.append(atom)
            last_atom = atom
        return residues
        
    def loadModel(self, model=1):
        self.model = model

    def getResidues(self, model):
        '''
        Gets residue coordinates.
        ** Returns: dict, where key = residue_id, and value = list of residue XYZ coordinates (list) 
        '''
        residues = []
        last_residue_id_chain = None
        data = []
        i = 0
        for line in self.structureReader():
            line_parts = line.split()
            if (len(line_parts) == 0):
                continue
            if ('SOL' in line) or ('SWM4' in line):
                continue
            residue_name = line_parts[3]
            if residue_name not in self.ignore:
                residue_number = line_parts[5]
                chain = line_parts[4]
                residue_id = residue_name + residue_number
                residue_id_chain = residue_id + chain
                if (residue_id_chain == last_residue_id_chain) or (last_residue_id_chain is None):
                    data.append(line_parts)
                else:
                    residue = Residue(data=data, index=i)
                    residues.append(residue)
                    i += 1
                    data = []
                    data.append(line_parts)
                last_residue_id_chain = residue_id_chain
        residue = Residue(data=data, index=i)
        residues.append(residue)
        return residues
    
    @property
    def residues(self):
        # return [res for res in self._residues if res.model == self.model]
        for residue in self._residues:
            if residue.model == self.model:
                yield residue
    
    def residue(self, index):
        return list(self.residues)[index]
    
    @property
    def atoms(self):
        for residue in self.residues:
            for atom in residue.atoms:
                yield atom
    
    def atom(self, index):
        return list(self.atoms)[index]
    
    @property
    def n_chains(self):
        return len(np.unique([res.chain_index for res in self.residues]))

    @property
    def chains(self):
        for i in range(self.n_chains):
            yield [res for res in self.residues if (res.chain_index == i)]

    @property
    def lipids(self):
        return [res for res in self.residues if res.restype == 'lipid']
    
    def select(self, subset=None, res_name=None, res_index=None, res_id=None, res_type=None):
        if subset is None:
            subset = self.residues
        elif (isinstance(subset, str)):
            if subset.startswith('residue'):
                subset = self.residues
            elif (subset.startswith('atom')) or (subset.startswith('protein')):
                subset = self.atoms
            elif subset.startswith('lipid'):
                subset = self.lipids
            else:
                raise ValueError('you are shit out of luck lol')
        else:
            subset = subset
        if res_name is not None:
            return [item for item in subset if item.name == res_name]
        if res_index is not None:
            return [item for item in subset if item.index == res_index]
        if res_id is not None:
            return [item for item in subset if item.id == res_id]
        if res_type is not None:
            return [item for item in subset if item.restype == res_type]
        


class Residue:

    def __init__(self, data):
        self.atomic_masses = {
            'C':12.011,
            'A':12.011,
            'N':14.007,
            'H':1.008,
            'O':15.999,
            'S':32.06,
            'P':30.974
        }
        self.getResidueData(data)
        # self.restype = None
        # if self.isResidue(data[0]):
        #     self.getResidueData(data)
        # elif self.isNucleic(data[0]):
        #     self.getResidueData(data)
        # elif self.isLipid(data[0]):
        #     self.getResidueData(data)
        # else:
        #     # print('idk what this residue is :(')
        #     self.getResidueData(data)

    @property
    def restype(self):
        if self.name in library.residues():
            return 'residue'
        elif self.name in library.dna():
            return 'nucleic'
        elif self.name.startswith(tuple(library.lipid())):
            return 'lipid'
        else:
            return None
        
    
    @staticmethod
    def linecheck(line_parts):
        if not line_parts[4].isalpha():
            chain_id = ''.join([char for char in line_parts[4] if char.isalpha()])
            resnum = ''.join([char for char in line_parts[4] if char.isnumeric()])
            new_line_parts = []
            for k, item in enumerate(line_parts):
                if k < 4:
                    new_line_parts.append(item)
                elif k == 4:
                    new_line_parts.append(chain_id)
                    new_line_parts.append(resnum)
                else:
                    new_line_parts.append(item)
            return new_line_parts
        return line_parts

    def getResidueData(self, data):
        _data = data[0]
        self.name = _data.residue_name
        self.number = _data.residue_number
        self.num = _data.residue_number
        self.index = _data.residue_index
        self.id = _data.residue_id
        self.chain = _data.chain
        self.chain_index = _data.chain_index
        self.atoms = []
        self.model = _data.model
        for _data in data:
            self.atoms.append(Atom(_data, self))

    def getHydrogenConnections(self):
        for hydrogen in self.hydrogens:
            for atom in self.atoms:
                if (not hydrogen is atom):
                    if 'N' in atom.name:
                        if  (self.distance(atom,hydrogen) > 0.98) and (self.distance(atom, hydrogen) < 1.1):
                            hydrogen.hbond = True
                            hydrogen.electronegative = atom
                    if 'O' in atom.name:
                        if (self.distance(atom,hydrogen) > 0.94) and (self.distance(atom,hydrogen) < 1.0):
                            hydrogen.hbond = True
                            hydrogen.electronegative = atom

    @property
    def com(self):
        x = []
        y = []
        z = []
        for atom in self.atoms:
            mass = self.atomic_masses[atom.type]
            x.append(atom.coordinates[0] * mass)
            y.append(atom.coordinates[1] * mass)
            z.append(atom.coordinates[2]* mass)
        com_x = sum(x) / self.mass
        com_y = sum(y) / self.mass
        com_z = sum(z) / self.mass
        com = (com_x, com_y, com_z)
        return com
    
    @property
    def mass(self):
        residue_mass = 0
        for atom in self.atoms:
            mass = self.atomic_masses[atom.type]
            residue_mass = residue_mass + mass
        return residue_mass
    @staticmethod
    def distance(x,y):
        return np.sqrt((x.coordinates[0]-y.coordinates[0])**2 + (x.coordinates[1] - y.coordinates[1])**2 + (x.coordinates[2] - y.coordinates[2])**2)

    def _min(self, dim):
        coords = []
        for atom in self.atoms:
            coords.append(atom.coordinates)
        return np.min(np.array(coords)[dim])
    
    def _max(self, dim):
        coords = []
        for atom in self.atoms:
            coords.append(atom.coordinates)
        return np.max(np.array(coords)[dim])
class Atom:
    def __init__(self, data, residue):
        self.residue = residue
        self.name = data.atom_name
        self.number = data.atom_number
        self.index = data.atom_index
        self.coordinates = (data.x, data.y, data.z)
        self.electronegative = None
        self.type = data.elem
        self.elem = data.elem
        if self.elem == 'H':
            self.hbond = True
        else:
            self.hbond = False