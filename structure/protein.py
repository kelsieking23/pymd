####################################################################################################################
# protein.py                                                                                                       #
# Written by Kelsie M. King                                                                                        #
# Last Edited 10/21/2020                                                                                           #
# Protein contains the protein class, which holds data related to docking receptor.                                #
####################################################################################################################
import os
import numpy as np
import sys

# imports from pymd
sys.path.append(os.getcwd())
from pymd.structure.ligand import Ligand
from pymd.utilities.rewritepdb import addChainID
from pymd.utilities.rewritepdb import writePDB
from pymd.utilities.rewritepdb import editChainIDResidue
from pymd.utilities.pdbhelper import fixBadCoordinates
from pymd.mdanalysis.multixpm import MultiXPM
from pymd.utilities import library


class Protein:

    def __init__(self, structure, id=None, residues=None, chain=None, ignore=None, covalent_ids=None, ligands=None):
        '''
        Protein class contains structure data for protein, such as residue ID's, residue coordinates, residue center of mass.
        Arguments:
        ** structure: string. Path to .pdb or .pdbqt file of protein structure.
        ** id (optional): string. an ID to label the protein with. 
        ** residues (optional): list. Residues to get coordinates of. If none are specified, will get all protein coordinates. 
        ** ignore (optional): list. List of residue ID's to ignore (i.e. a covalent ligand) (Note: automatically ignores solvent and other artifacts)
        ** covalent_ids (optional): list. List of residue ID's of covalent ligands.
        ** ligands(optional): str or list. If given, covalent ligands/
        '''
        self.structure = os.path.abspath(structure)
        self.ligands = []
        # add covalent ligands to ignore list; will be dealt with separately
        if ignore is not None:
            if covalent_ids is not None:
                for item in covalent_ids:
                    if item not in ignore:
                        ignore.append(item)
                    else:
                        continue
        if ignore is not None:
            self.ignore = ignore
        else:
            self.ignore = []
        if ignore is None:
            ignore = []

        if chain is not None:
            ignore.append('notChain{}'.format(chain))
        self.ligand_obj = []
        if ligands is not None:
            if isinstance(ligands, list):
                for item in ligands:
                    ignore.append(item)
                    self.ignore.append(item)
                    self.ligands.append(item)
            else:
                self.ligands.append(ligands)
                ignore.append(ligands)
                self.ignore.append(ligands)
                lig = Ligand(structure=self.structure, name=ligands)
                self.ligand_obj.append(lig)
        
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
        # check structure file
        extension = self.structure[-3:]
        if extension.lower() == 'gro':
            self.structure = self.groHandler()
        # self.structure = self.checkChainID()
        # get residue data
        if residues is None:
            self.ids, ligand_structure = self.getResidueIds(ignore)
        else:
            self.ids = residues
        self.coordinates = self.getResidueCoordinates()
        self.atom_coordinates = self.getAtomCoordinates()
        self.residues = self.getResidues()
        self.atoms = self.getAtoms()
        self.chains = self.getChains()
        self.peptides = len(self.chains.keys())
        # self.structure_data = self.getStructureData()
        self.atom_types = self.getAtomTypes()
        self.masses = self.getResidueMasses()
        
        self.coms = self.getResidueCOM()
        # self.gro = None

        # if ligand is covalently bound, send structure to Ligand class for processing
        if covalent_ids is not None:
            if self.ligands is not None:
                self.ligands = []
            if isinstance(covalent_ids, list):
                for covalent_id in covalent_ids:
                    self.ligands.append(Ligand(structure=self.structure, covalent_id=covalent_id, name=covalent_id))
            if isinstance(covalent_ids, str):
                self.ligands.append(Ligand(structure=self.structure, covalent_id=covalent_ids, name=covalent_ids))

        self.id = id
    @property
    def molecularWeight(self):
        mw = 0
        residues = self.getResidues()
        for residue in residues:
            for atom in residue.atoms:
                mw += self.atomic_masses[atom.type]
        return mw

    @property
    def molecularWeightPerChain(self):
        return self.molecularWeight / len(self.chains.keys())

    @property
    def glycines(self):
        glycines = []
        for chain in self.chains.keys():
            i = 0
            indeces = self.chains[chain]['indeces']
            ids = self.chains[chain]['ids']
            for res_id in ids:
                if 'GLY' in res_id:
                    glycines.append(indeces[i])
                i += 1
        return glycines

    @property
    def nterm(self):
        return self.getTermini()[0]

    @property
    def cterm(self):
        return self.getTermini()[1]

    @property
    def termini(self):
        return self.getTermini()

    @property
    def backbone(self):
        ca = {}
        for chain in self.chains.keys():
            ca[chain] = []
            for residue in self.chains[chain]:
                for atom in residue.atoms:
                    if (atom.name == 'CA') or (atom.name == 'C') or (atom.name == 'N'):
                        ca[chain].append(atom)
        return ca

    def structureReader(self):
        f = open(self.structure, 'r')
        contents = f.readlines()
        f.close()
        for line in contents:
            yield line

    def groHandler(self):
        f = open(self.structure, 'r')
        contents = f.readlines()
        f.close()
        data = []
        i = 0
        k = 0
        residues = []
        sm = None
        for line in contents:
            if i < 2:
                i += 1
                continue
            elif line == contents[-1]:
                break
            else:
                line_parts = line.split()
                res_name = line_parts[0][-3:]
                res_num = line_parts[0][:-3]
                res_id = res_name + res_num
                if (res_name == 'SOL') or (res_name == 'SWM4'):
                    data.append('TER')
                    break
                if self.ligands is not None:
                    if res_name in self.ligands:
                        if k == 0:
                            k = 1
                            sm = res_name
                        else:
                            pass
                    else: residues.append(res_name)
                else:
                    residues.append(res_name)
                atom_type = line_parts[1]
                atom_num = line_parts[2]
                x, y, z = line_parts[3:6]
                x = float(x) * 10
                y = float(y) * 10
                z = float(z) * 10
                newline = ['ATOM', atom_num, atom_type, res_name, 'X', res_num, x, y, z, '1.00', '0.00']
                data.append(newline)
        data.append(['TER'])
        newfilename = self.structure[:-3] + 'pdb'
        writePDB(data, newfilename)
        nterm = residues[0]
        cterm = residues[-1]
        editChainIDResidue(newfilename, newfilename, nterm, cterm, sm)
        return newfilename

    def checkChainID(self):
        linecheck = None
        for line in self.structureReader():
            line_parts = line.split()
            if line_parts == []:
                continue
            if 'ATOM' in line_parts[0]:
                if ('SOL' not in line) and ('SWM4' not in line):
                    linecheck = line_parts
                    break
        if linecheck is not None:
            if not linecheck[4].isalpha():
                nterm, cterm = self.getTermini()
                return self.addChainID(nterm=nterm, cterm=cterm)
            elif linecheck[4] == 'X':
                nterm, cterm = self.getTermini()
                return editChainIDResidue(self.structure, self.structure, nterm=nterm, cterm=cterm)
            else:
                return self.structure
        else:
            return self.structure

    def getTermini(self):
        residue_ids = []
        for line in self.structureReader():
            line_parts = line.split()
            if line.startswith('ATOM'):
                if ('SOL' in line) or ('SWM4' in line):
                    break
                residue_name = line_parts[3]
                if line_parts[4].isalpha():
                    residue_num = line_parts[5]
                else:
                    residue_num = line_parts[4]
                residue_id = residue_name + residue_num
                if residue_name not in self.ligands:
                    if (self.ignore is not None) and (residue_name not in self.ignore):
                        if residue_id not in residue_ids:
                            residue_ids.append(residue_id)
                    else:
                        if residue_id not in residue_ids:
                            residue_ids.append(residue_id)
        nterm = residue_ids[0][:3]
        cterm = residue_ids[-1][:3]
        return (nterm, cterm)
    

    def addChainID(self, nterm, cterm):
        new_filename = self.structure.split('.')[0] + '_chainid.pdb'
        return addChainID(self.structure, new_filename, nterm, cterm)

    def getResidueIds(self, ignore):
        '''
        Get residue ID's for each residue in the PDB. Ignores residues passed in ignore argument, if any. 
        ** Returns: list
        '''
        valid_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD']
        residue_ids = []
        ligands = []
        residue_index = 1
        last_residue_id = None
        residue_id = None
        for line in self.structureReader():
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' in line_parts[0]:
                    if ('SOL' in line) or ('SWM4' in line):
                        break
                    residue_name = line_parts[3]
                    if residue_name not in valid_residues:
                        continue
                    residue_number = line_parts[5]
                    residue_id = residue_name + residue_number
                    chain_id = line_parts[4]
                    if last_residue_id is not None:
                        if last_residue_id != residue_id:
                            residue_index += 1
                    if ignore is not None:
                        if residue_id in ignore:
                            if self.ligands is not None:
                                if residue_id in self.ligands:
                                    ligands.append(line_parts)
                                    continue
                        elif residue_name in ignore:
                            if self.ligands is not None:
                                if residue_name in self.ligands:
                                    ligands.append(line_parts)
                                    continue
                            else:
                                continue
                        elif len([x for x in ignore if x.startswith('notChain')]) > 0:
                            ignore_chain = [x for x in ignore if x.startswith('notChain')][0][-1]
                            if chain_id != ignore_chain:
                                continue
                            else:
                                if residue_id not in residue_ids:
                                    residue_ids.append(residue_id)
                        else:
                            if residue_id not in residue_ids:
                                residue_ids.append(residue_id)
                    elif residue_id not in residue_ids:
                        residue_ids.append(residue_id)
                    else:
                        continue
                else:
                    continue
            if residue_id is not None:
                last_residue_id = residue_id
        if len(ligands) != 0:
            ligands.append(['TER'])
            ligand_filename = os.path.join(os.path.dirname(self.structure), 'ligand.pdb')
            ligand_pdb = writePDB(ligands, ligand_filename)
        else:
            ligand_pdb = None
        ligand_pdb = None
        return residue_ids, ligand_pdb

    def getResidueCoordinates(self):
        '''
        Gets residue coordinates.
        ** Returns: dict, where key = residue_id, and value = list of residue XYZ coordinates (list) 
        '''
        residue_coords = {}
        for line in self.structureReader():
            if 'ENDMDL' in line:
                break
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' == line_parts[0]:
                    if ('SOL' in line) or ('SWM4' in line):
                        break
                    current_residue = line_parts[3]
                    current_residue_number = line_parts[5]
                    residue_id = current_residue + current_residue_number
                    if residue_id in self.ids:
                        if residue_id not in residue_coords.keys():
                            residue_coords[residue_id] = []
                        try:
                            current_res_coords = list(map(float, line_parts[6:9]))
                        except:
                            current_res_coords = fixBadCoordinates(line_parts[6:9])
                        residue_coords[residue_id].append(current_res_coords)
                    else:
                        continue
                else:
                    continue
        filtered = {}
        for res_id in self.ids:
            filtered[res_id] = residue_coords[res_id]
        return residue_coords
    

    def getChains(self):
        chains = {}
        for residue in self.residues:
            if residue.chain not in chains.keys():
                chains[residue.chain] = []
            chains[residue.chain].append(residue)
        return chains
        # chains = {}
        # last_residue_id = None
        # residue_index = 1
        # for line in self.structureReader():
        #     line_parts = line.split()
        #     if 'TER' not in line_parts:
        #         try:
        #             residue_name = line_parts[3]
        #         except:
        #             continue
        #         if 'ATOM' in line_parts[0]:
        #             residue_number = line_parts[5]
        #             residue_id = residue_name + residue_number
        #             if (residue_id in self.ignore) or (residue_name in self.ignore):
        #                 if (residue_id in self.ligands) or (residue_name in self.ligands):
        #                     break
        #                 else: 
        #                     continue
        #             else:
        #                 if last_residue_id is not None:
        #                     if last_residue_id != residue_id:
        #                         residue_index += 1
        #                 chain_id = line_parts[4]
        #                 if chain_id not in chains.keys():
        #                     chains[chain_id] = {
        #                         'indeces':[],
        #                         'ids':[]
        #                     }
        #                 if residue_index not in chains[chain_id]['indeces']:
        #                     chains[chain_id]['indeces'].append(residue_index)
        #                 if residue_id not in chains[chain_id]['ids']:
        #                     chains[chain_id]['ids'].append(residue_id)
        #                 last_residue_id = residue_id
        # return chains

    def renumber(self, chain, start):
        new_chain = []
        for residue in self.chains[chain]:
            name = residue.name
            residue.number = start
            residue.id = '{}{}'.format(name, start)
            new_chain.append(residue)
            start += 1
        self.chains[chain] = new_chain
        return self.chains

    def getAtomCoordinates(self):
        atom_coords = {}
        for line in self.structureReader():
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' == line_parts[0]:
                    if ('SOL' in line) or ('SWM4' in line):
                        break
                    atom_num = int(line_parts[1])
                    residue_name = line_parts[3]
                    residue_num = line_parts[5]
                    res_id = residue_name + residue_num
                    if res_id in self.ids:
                        if line_parts[4].isalpha():
                            try:
                                atom_coords[atom_num] = list(map(float, line_parts[6:9]))
                            except:
                                atom_coords[atom_num] = fixBadCoordinates(line_parts[6:9])
                        else:
                            try:
                                atom_coords[atom_num] = list(map(float, line_parts[5:8]))
                            except:
                                atom_coords[atom_num] = fixBadCoordinates(line_parts[5:8])
                else:
                    continue
            else:
                continue
        return atom_coords

    def getStructureData(self):
        data = []
        for line in self.structureReader():
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' in line_parts[0]:
                    data.append(line_parts)
        return data 

    def getResidueMasses(self):
        '''
        Gets residue atomic masses. 
        ** Returns: dict, where key = residue_id, and value = atomic mass of residue
        '''
        residue_masses = {}
        for line in self.structureReader():
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' in line_parts[0]:
                    if ('SOL' in line) or ('SWM4' in line):
                        break
                    residue_name = line_parts[3]
                    residue_number = line_parts[5]
                    residue_id = residue_name + residue_number
                    if residue_id in self.ids:
                        atom_type = line_parts[-1][0]
                        if not atom_type.isalpha():
                            atom_type = line_parts[2]
                            for char in atom_type:
                                if char.isalpha():
                                    atom_type = char
                                    break
                        mass = self.atomic_masses[atom_type]
                        if residue_id not in residue_masses.keys():
                            residue_masses[residue_id] = mass
                        else:
                            residue_masses[residue_id] += mass
                    else:
                        continue
                else:
                    continue
        return residue_masses
    
    def getAtomTypes(self):
        atom_types = {}
        i = 0
        for line in self.structureReader():
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' in line_parts[0]:
                    if ('SOL' in line) or ('SWM4' in line):
                        break
                    residue_name = line_parts[3]
                    residue_number = line_parts[5]
                    residue_id = residue_name + residue_number
                    if residue_id in self.ids:
                        i += 1
                        atom_type = line_parts[-1][0]
                        if not atom_type.isalpha():
                            atom_type = line_parts[2]
                            for char in atom_type:
                                if char.isalpha():
                                    atom_type = char
                                    break
                        if residue_id not in atom_types.keys():
                            atom_types[residue_id] = [atom_type]
                        else:
                            atom_types[residue_id].append(atom_type)
                    else:
                        continue
                else:
                    continue
        for res_id in atom_types.keys():
            atom_types[res_id] = list(enumerate(atom_types[res_id]))
        return atom_types
    
    def getCOM(self, residues=[]):
        residue_coms = {}
        x = []
        y = []
        z = []
        i = 0
        residue_mass = 0
        for residue in residues:
            residue_id = residue.id
            residue_mass += self.masses[residue_id]
            for atom in residue.atoms:
                atom_type = atom.type
                mass = self.atomic_masses[atom_type]
                x.append(atom.coordinates[0] * mass)
                y.append(atom.coordinates[1] * mass)
                z.append(atom.coordinates[2]* mass)
                i += 1
        com_x = sum(x) / residue_mass
        com_y = sum(y) / residue_mass
        com_z = sum(z) / residue_mass
        com = (com_x, com_y, com_z)
        return com
    
    def getResidueCOM(self):
        '''
        Gets residue center of mass for all residues. 
        ** Returns: dict, where key = residue_id, and value = xyz coordinates of COM (tuple)
        '''
        residue_coms = {}
        for residue_id in self.ids:
            residue_coordinates = self.coordinates[residue_id]
            residue_mass = self.masses[residue_id]
            x = []
            y = []
            z = []
            i = 0
            for coordinate in residue_coordinates:
                atom_type = self.atom_types[residue_id][i][1]
                mass = self.atomic_masses[atom_type]
                x.append(coordinate[0] * mass)
                y.append(coordinate[1] * mass)
                z.append(coordinate[2]* mass)
                i += 1
            com_x = sum(x) / residue_mass
            com_y = sum(y) / residue_mass
            com_z = sum(z) / residue_mass
            com = (com_x, com_y, com_z)
            residue_coms[residue_id] = com
        return residue_coms

    def getProteinCOM(self):
        x = []
        y = []
        z = []
        total_mass = 0
        for residue in self.coordinates:
            residue_coordinates = self.coordinates[residue]
            atom_types = self.atom_types[residue]
            total_mass += self.masses[residue]
            for i in range(len(residue_coordinates)):
                coordinate = residue_coordinates[i]
                atom_type = atom_types[i][1]
                mass = self.atomic_masses[atom_type]
                x.append(coordinate[0] * mass)
                y.append(coordinate[1] * mass)
                z.append(coordinate[2]* mass)
        com_x = sum(x) / total_mass
        com_y = sum(y) / total_mass
        com_z = sum(z) / total_mass
        com = (com_x, com_y, com_z)
        return com
    
    @property
    def proteinCOM(self):
        return self.getProteinCOM()

    @property
    def residueIndeces(self):
        resi = 0
        last_residue_id = None
        indeces = {}
        for line in self.structureReader():
            if not line.startswith('ATOM'):
                continue
            if ('SOL' in line) or ('SWM4' in line):
                break
            line_parts = line.split()
            residue_name = line_parts[3]
            residue_number = line_parts[5]
            residue_id = residue_name + residue_number
            if self.ignore is not None:
                if (residue_name in self.ignore) or (residue_id in self.ignore) or (residue_number in self.ignore):
                    continue
            if residue_id != last_residue_id:
                # indeces[resi] = residue_name
                indeces[resi] = residue_id
                resi += 1
            last_residue_id = residue_id
        return indeces
            


    def getResidues(self):
        residues = []
        residue_index = 0
        last_residue_id_chain = None
        data = []
        i = 0
        for line in self.structureReader():
            line_parts = line.split()
            if not line.startswith('ATOM'):
                continue
            if ('SOL' in line) or ('SWM4' in line):
                break
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

    def getAtoms(self):
        atoms = []
        for residue in self.residues:
            for atom in residue.atoms:
                atoms.append(atom)
        return atoms
    def ligandInteraction(self):
        coms = self.getResidueCOM()
    
    @staticmethod
    def mindistAtoms(x, y):
        _d = None
        for xatom in x:
            for yatom in y:
                d = np.linalg.norm(np.array(xatom.coordinates)-np.array(yatom.coordinates))
                if (_d is None) or d < _d:
                    _d = d
        return _d
    
    @staticmethod
    def distance(x,y):
        return np.sqrt((x[0]-y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2)


class Residue:

    def __init__(self, data, index):
        self.index = index
        first_line = data[0]
        self.atomic_masses = {
            'C':12.011,
            'A':12.011,
            'N':14.007,
            'H':1.008,
            'O':15.999,
            'S':32.06,
            'P':30.974
        }
        self.restype = None
        if self.isResidue(data):
            self.getResidueData(data)
        elif self.isNucleic(data):
            self.getResidueData(data)
        elif self.isLipid(data):
            self.getResidueData(data)
        else:
            # print('idk what this residue is :(')
            self.getResidueData(data)
    def isResidue(self, data):
        residue_name = data[0][3]
        if residue_name not in library.residues():
            return False
        self.restype = 'residue'
        return True
        
    def isNucleic(self, data):
        residue_name = data[0][3]
        if residue_name not in library.dna():
            return False
        self.restype = 'nucleic'
        return True

    def isLipid(self, data):
        residue_name = data[0][3]
        if residue_name.startswith(tuple(library.lipid())):
            return False
        self.restype = 'lipid'
        return True
    
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
        first_line = data[0]
        first_line_clean = self.linecheck(first_line)
        self.name = first_line_clean[3]
        try:
            self.number = int(first_line_clean[5])
        except:
            self.number = int(first_line_clean[4])
        self.num = self.number
        self.id = self.name + str(self.number)
        self.chain = first_line_clean[4]
        self.atoms = []
        self.hydrogens = []
        for line_parts in data:
            atom = Atom(data=self.linecheck(line_parts), residue=self)
            self.atoms.append(atom)
            if 'H' in atom.name:
                self.hydrogens.append(atom)
        self.getHydrogenConnections()

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

    


class Atom:
    def __init__(self, data, residue):
        self.residue = residue
        self.name = data[2]
        self.number = int(data[1])
        self.index = self.number-1
        if not data[4].isnumeric():
            try:
                self.coordinates = list(map(float, data[6:9]))
            except:
                self.coordinates = fixBadCoordinates(data[6:9])
        else:
            try:
                self.coordinates = list(map(float, data[5:8]))
            except:
                self.coordinates = fixBadCoordinates(data[5:8])

        self.electronegative = None
        self.type = data[-1][0]

        if self.name == 'H':
            self.hbond = True
        else:
            self.hbond = False

# path = os.path.abspath('../mdanalysis/test/md_500_600.gro')
# print(os.path.dirname(path))
# pro = Protein(os.path.abspath('../mdanalysis/test/md_500_600.gro'), ligands='MYR', ignore=['ACE', 'NH2'])
# print(pro.getChains())