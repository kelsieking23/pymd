####################################################################################################################
# ligand .py                                                                                                       #
# Written by Kelsie M. King                                                                                        #
# Last Edited 06/03/2020                                                                                           #
# Ligand contains Ligand class, which holds data related to docked structures.                                     #
####################################################################################################################
import os
import numpy as np
from pymd.structure.protein import Atom

class Ligand:

    def __init__(self, structure, name=None, model=None, ignh=False):
        self.structure = structure
        self.name = name
        self.ignh = ignh
        model_coords = self.split(ignh)
        if model is not None:
            if isinstance(model_coords, dict):
                self.coordinates = model_coords[model]
            if isinstance(model_coords, list):
                # print('WARNING: model {} was specified, but only one model found in structure {}.'.format(model, structure))
                self.coordinates = model_coords
        else:
            model = 'MODEL 1'
            if isinstance(model_coords, dict):
                self.coordinates = model_coords[model]
            if isinstance(model_coords, list):
                # print('WARNING: model {} was specified, but only one model found in structure {}.'.format(model, structure))
                self.coordinates = model_coords
        self.coordinates = self.getAtomCoordinates()
        self.covalent_id = None
        self.atomic_masses = {
            'C':12.011,
            'A':12.011,
            'N':14.007,
            'H':1.008,
            'O':15.999,
            'S':32.06
        }
        # self.log = self.getLog()
    @property
    def log(self):
        if self.getLog() == 0:
            return None
        else:
            return self.getLog()
            
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

    def getLog(self):
        log_path = None
        directory = os.path.dirname(self.structure)
        files = os.listdir(directory)
        for filename in files:
            if 'log' in filename:
                log_path = os.path.join(directory, filename)
                return log_path
        return log_path

    @staticmethod
    def convertGroLine(line):
        pdb = ['ATOM'] # ATOM
        line_parts = line.split()
        atom_num = line_parts[2] # atom num
        pdb.append(atom_num)
        atom_name = line_parts[1] # atom name
        pdb.append(atom_name)
        res_name = line_parts[0][-3:] # res name
        pdb.append(res_name)
        pdb.append('X') # chain ID
        res_num = line_parts[0][:-3] #res num
        pdb.append(res_num) 
        for coord in line_parts[3:6]:
            c = float(coord) * 10
            c = round(c,3)
            pdb.append(str(c))
        return pdb

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
            print(atom.__dict__)
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

    def getAtomCoordinates(self, data=None):
        if data is None:
            data = open(self.structure, 'r')
        atom_coords = {}
        self.atoms = []
        atom_types = ('ATOM', 'HETATM')
        for line in data:
            line_parts = line.split()
            if not line.startswith(atom_types):
                continue
                # line_parts = self.convertGroLine(line)
            if len(line_parts) > 0:
                atom_num = int(line_parts[1])
                if line_parts[4].isalpha():
                    atom_coords[atom_num] = list(map(float, line_parts[6:9]))
                else:
                    atom_coords[atom_num] = list(map(float, line_parts[5:8]))
            else:
                continue
            if line.startswith(atom_types):
                self.atoms.append(Atom(line_parts, self))
        return atom_coords
    
    def getAtomTypes(self):
        atom_types = {}
        for line in self.lig.data:
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' == line_parts[0]:
                    if 'SOL' in line:
                        break
                    atom_num = int(line_parts[1])
                    if line_parts[2].startswith('CA'):
                        atom_types[atom_num] = 'CA'
                    else:
                        atom_types[atom_num] = line_parts[2][0]
                else:
                    continue
            else:
                continue
        return atom_types
    
    @property
    def atom_types(self):
        return self.getAtomTypes()



    def getAtomNames(self):
        atom_names = {}
        for line in self.lig.data:
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' == line_parts[0]:
                    if 'SOL' in line:
                        break
                    atom_num = int(line_parts[1])
                    atom_names[atom_num] = line_parts[2]
                else:
                    continue
            else:
                continue
        return atom_names

    @property
    def atom_names(self):
        return self.getAtomNames()


    def bonding(self):
        connect, dist = self.connect()
        b = Bonding(connect)
        x = b.findCarbonRings()
        print(x)
        y = b.findCarbonRings(start=x[-1], ref_atom_num=None, visited=x)
        return y

    def connect(self):
        for ligand in self.ligands:
            dist = {}
            for ref_atom_num, ref_atom_coord in ligand.coordinates.items():
                for atom_num, atom_coord in ligand.coordinates.items():
                    d = round(self.distance(atom_coord, ref_atom_coord),3)
                    if ref_atom_num not in dist.keys():
                        dist[ref_atom_num] = {}
                    if atom_num not in dist[ref_atom_num].keys():
                        dist[ref_atom_num][atom_num] = d
            connect = {}
            for ref_atom_num in dist.keys():
                connect[ref_atom_num] = {}
                for atom_num in dist[ref_atom_num].keys():
                    d = dist[ref_atom_num][atom_num]
                    if (d > 1.38) and (d < 1.4):
                        connect[ref_atom_num][atom_num] = 'CA-CA'
                    if (d > 1.22) and (d < 1.23):
                        connect[ref_atom_num][atom_num] = 'C=O'
                    if (d > 1.36) and (d <1.37):
                        connect[ref_atom_num][atom_num] = 'C-O'
                    if (d > 0.99) and (d < 1.01):
                        if (self.atom_types[ref_atom_num] == 'O') or (self.atom_types[atom_num] == 'O'):
                            if (self.atom_types[ref_atom_num] == 'H') or (self.atom_types[atom_num] == 'H'):
                                connect[ref_atom_num][atom_num] = 'O-H'
                    if (d > 1) and (d < 1.1):
                        if (self.atom_types[ref_atom_num] == 'CA') or (self.atom_types[atom_num] == 'CA'):
                            if (self.atom_types[ref_atom_num] == 'H') or (self.atom_types[atom_num] == 'H'):
                                connect[ref_atom_num][atom_num] = 'CA-H'
            break
        return connect, dist

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
                        if line_parts[4].isalpha():
                            residue_number = line_parts[5]
                        else:
                            residue_number = line_parts[4]
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

    @staticmethod
    def distance(x,y):
        return np.sqrt((x[0]-y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2)
    
    @staticmethod
    def theta(v, w): 
        v = np.array(v)
        w = np.array(w)
        from numpy.linalg import norm
        return np.degrees(np.arccos(v.dot(w)/(norm(v)*norm(w))))

class Bonding:

    def __init__(self, connections):
        self.connections = connections

    def findCarbonRings(self, start=None, ref_atom_num=None, visited=None):
        if (start is None) and (ref_atom_num is None):
            visited = []
            for ref_atom_num in self.connections.keys():
                for atom_num in self.connections[ref_atom_num].keys():
                    if (self.connections[ref_atom_num][atom_num] == 'CA-CA') or (self.connections[ref_atom_num][atom_num] == 'CA-O'):
                        visited.append(ref_atom_num)
                        visited.append(atom_num)
                        return self.findCarbonRings(start=ref_atom_num, ref_atom_num=atom_num, visited=visited)
        elif (start is not None) and (ref_atom_num is None) and (visited is not None):
            for atom_num in self.connections[start].keys():
                if (self.connections[start][atom_num] == 'CA-CA') or (self.connections[start][atom_num] == 'CA-O'):
                    if atom_num in visited:
                        continue
                    else:
                        start = atom_num
                        print('!!!!!!!! {}'.format(start))
                        for anum in self.connections[atom_num].keys():
                            if (self.connections[start][anum] == 'CA-CA') or (self.connections[start][anum] == 'CA-O'):
                                if anum not in visited:
                                    visited = [start, anum]
                                    return self.findCarbonRings(start=start, ref_atom_num=anum, visited=visited)
        elif (start is not None) and (ref_atom_num is None) and (visited is None):
            visited = [start]
            for atom_num in self.connections[start].keys():
                if (self.connections[start][atom_num] == 'CA-CA') or (self.connections[start][atom_num] == 'CA-O'):
                        visited.append(atom_num)
                        return self.findCarbonRings(start=start, ref_atom_num=atom_num, visited=visited)
        else:
            print('******{}'.format(start))
            for atom_num in self.connections[ref_atom_num].keys():
                if (self.connections[ref_atom_num][atom_num] == 'CA-CA') or (self.connections[ref_atom_num][atom_num] == 'CA-O'):
                    if (atom_num not in visited):
                        visited.append(atom_num)
                        return self.findCarbonRings(start=start, ref_atom_num=atom_num, visited=visited)
                    if (atom_num == start) and (len(visited) == 6):
                        print(start)
                        return visited
                    if (atom_num in visited):
                        continue
                        # continue
                        # visited.append(atom_num)
                        # return self.findCarbonRings(start=ref_atom_num, ref_atom_num=atom_num, visited=visited)



class SingleLigand(Ligand):

    def __init__(self, data):
        self.data = data
        self.coordinates = self.getAtomCoordinates(data)
        # self.getAtomTypes(data)

    @property
    def atom_types(self):
        return self.getAtomTypes(self.data)




class PredockLigand:

    def __init__(self, structure):
        
        self.structure = structure

