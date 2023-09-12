import numpy as np
import pandas as pd
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import sys
from pathlib import Path

from pymd.structure.protein import Protein
from pymd.structure.solvent import Solvent, SolventMolecule, SolventAtom
from pymd.utilities.rewritepdb import writePDB
from pymd.utilities.gro import convertGro
from pymd.utilities.rewritepdb import editChainIDResidue


class SolventAnalyzer:

    def __init__(self, structure, ligands=None, trim=True):
        if structure.endswith('gro'):
            self.structure = convertGro(structure, ligands=None)
        else:
            self.structure = structure

        structure, self.solvent = self.sortCoordinates()
        self.protein = Protein(structure, ignore=['ACE', 'NH2'])
        # try:
        self.structure = self.protein.structure
        self.hull = self.getConvexHull()
        # except:
        #     self.hull = None
        #     print('No protein given; will proceed with solvent only.')
        self.solvent_coordinates, self.trimmed = self.getSolventCoordinates(trim=trim)
        if trim == True:
            self.solvent = self.trimmed
        # self.mesh = pymesh.form_mesh(self.hull.vertices, self.hull.simplices)
    
    def sortCoordinates(self):
        protein = []
        solvent = []
        start = 0
        f = open(self.structure, 'r')
        for line in f:
            line_parts = line.split()
            if len(line_parts) > 0:
                if 'ATOM' in line_parts[0]:
                    if ('SOL' in line) or ('SWM4' in line):
                        solvent.append(line_parts)
                    elif ('NA' in line_parts[2]) or ('CL' in line_parts[2]) or ('POT' in line_parts[2]):
                        solvent.append(line_parts)
                    elif ('DPOT' in line_parts[2]):
                        continue
                    else:
                        protein.append(line_parts)
        pro_newfilename = self.structure.split('.')[0] + '_protein.pdb'
        writePDB(protein, pro_newfilename)
        sol_newfilename = self.structure.split('.')[0] + '_sol.pdb'
        writePDB(solvent, sol_newfilename)
        f.close()
        return pro_newfilename, sol_newfilename
    
    def getBoundries(self):
        x = None
        y = None
        z = None
        minx = None
        maxx = None
        miny = None
        maxy = None
        minz = None
        maxz = None
        pairs_x = None
        pairs_y = None
        pairs_z = None
        for atom_num_1, coord_1 in self.protein.atom_coordinates.items():
            for atom_num_2, coord_2 in self.protein.atom_coordinates.items():
                dist_x = coord_2[0] - coord_1[0]
                dist_y = coord_2[1] - coord_1[1]
                dist_z = coord_2[2] - coord_1[2]
                if atom_num_1 != atom_num_2:
                    if (x is None) or (dist_x > x):
                        minx = dist_x
                        pairs_x = (atom_num_1, atom_num_2)
                        if coord_2[0] > coord_1[0]:
                            maxx = coord_2[0]
                            minx = coord_1[0]
                        else:
                            maxx = coord_1[0]
                            minx = coord_2[0]
                    if (y is None) or (dist_y > y):
                        y = dist_y
                        pairs_y = (atom_num_1, atom_num_2)
                        if coord_2[1] > coord_1[1]:
                            maxy = coord_2[1]
                            miny = coord_1[1]
                        else:
                            maxy = coord_1[1]
                            miny = coord_2[1]
                    if (z is None) or (dist_z > z):
                        z = dist_z
                        pairs_z = (atom_num_1, atom_num_2)
                        if coord_2[2] > coord_1[2]:
                            maxz = coord_2[2]
                            minz = coord_1[2]
                        else:
                            maxz = coord_1[2]
                            minz = coord_2[2]
        return (minx, maxx), (miny, maxy), (minz, maxz)



    def getSolventCoordinates(self, trim=True):
        solvent_coordinates = {}
        if trim == True:
            cx, cy, cz = self.getHullCOM()
        lines = []
        self.solvent_atom_index = {}

        if trim is True:
            for atom, coordinate in self.protein.atom_coordinates.items():
                i = 0
                _max = 0
                max_atom1 = None
                max_atom2 = None
                coords = [list(value) for value in self.protein.atom_coordinates.values()]
                atoms = [key for key in self.protein.atom_coordinates.keys()]
                for coord in coords:
                    _d =  np.sqrt((coordinate[0] - coord[0])**2 + (coordinate[1] - coord[1])**2 + (coordinate[2] - coord[2])**2)
                    if _d > _max:
                        _max = _d
                        max_atom1 = atom
                        max_atom2 = atoms[i]
                    i += 1
        self.solvent_obj = Solvent(self.solvent)
        solvent_coordinates = {}
        for mol in self.solvent_obj.molecules:
            for atom in mol.atoms:
                # solvent_coordinates[atom.id] = atom.coordinates
                coordinates = atom.coordinates
                if trim == True:
                    d = np.sqrt((coordinates[0] - cx)**2 + (coordinates[1] - cy)**2 + (coordinates[2] - cz)**2)
                    if d < _max:
                        solvent_coordinates[atom.id] = coordinates
                        lines.append(atom.data.split())
                else:
                    solvent_coordinates[atom.id] = coordinates
        # for line in f:
        #     line_parts = line.split()
        #     atom_num = int(line_parts[1])
        #     if atom_num in solvent_coordinates.keys():
        #         atom_num = int(last_atom_num) + 1
        #     coordinates = list(map(float, line_parts[5:8]))
        #     d = np.sqrt((coordinates[0] - cx)**2 + (coordinates[1] - cy)**2 + (coordinates[2] - cz)**2)
        #     if trim == True:
        #         if d < _max:
        #             solvent_coordinates[atom_num] = coordinates
        #             lines.append(line_parts)
        #     else:
        #         solvent_coordinates[atom_num] = coordinates
        #     self.solvent_atom_index[i] = atom_num
        #     last_atom_num = atom_num
        #     i += 1
        if trim == True:
            newfilename = self.structure.split('.')[0] + '_trimmed.pdb'
            writePDB(lines, newfilename)
        else:
            newfilename = None
        return solvent_coordinates, newfilename
        
    def getSolventMolecules(self):
        solvent_molecules = {}
        self.solvent_residue_index = {}
        f = open(self.solvent, 'r')
        contents = f.readlines()
        f.close()
        i = 0
        for line in contents:
            line_parts = line.split()
            atom_name = line_parts[2]
            res_num = line_parts[4]
            if res_num not in self.solvent_residue_index.keys():
                key = int(res_num)
            else:
                key = last_res_num + 1
            if (atom_name == 'OW') or (atom_name == 'OM'):
                self.solvent_residue_index[i] = int(res_num)
                solvent_molecules[key] = []
                solvent_molecules[key].append(line)
                solvent_molecules[key].append(contents[i+1])
                solvent_molecules[key].append(contents[i+2])
            elif (atom_name == 'HW1') or (atom_name == 'HW2') or (atom_name == 'H1') or (atom_name == 'H2'):
                i += 1
                continue
            else:
                self.solvent_residue_index[i] = res_num
                solvent_molecules[res_num] = [line]
            last_res_num = key
            i += 1
        return solvent_molecules
        
    def getAtomResidueIndecesMap(self):
        f = open(self.solvent, 'r')
        i = 0
        atom_residue_map = {}
        for line in f:
            line_parts = line.split()
            atom_num = int(line_parts[1])
            res_num = int(line_parts[4])
            atom_name = line_parts[2]
            if atom_num not in atom_residue_map.keys():
                key = atom_num
            else:
                key = last_atom_num + 1
            if atom_name == 'OW':
                atom_residue_map[key] = i
            if (atom_name == 'NA') or (atom_name == 'CL'):
                atom_residue_map[key] = i
            i +=1
            last_atom_num = key
        f.close()
        return atom_residue_map


    def getConvexHull(self):
        arr = []
        keys = []
        for key,value in self.protein.atom_coordinates.items():
            arr.append(value)
            keys.append(key)
        X = np.array(arr)
        hull = ConvexHull(X)
        to_write = []
        structure_data = self.protein.getStructureData()
        for point in hull.vertices:
            for line_parts in structure_data:
                atom_num = int(line_parts[1])
                if atom_num == point:
                    to_write.append(line_parts)
        newfilename = self.structure.split('.')[0] + '_vertices.pdb'
        writePDB(to_write, newfilename)
        return hull
    
    def getHullRadius(self):
        verts = self.hull.vertices
        matrix = {}
        largest_d = None
        largest_verts = None
        for vert in verts:
            if vert not in matrix.keys():
                matrix[vert] = {}
            point = self.hull.points[vert]
            for _vert in verts:
                _point = self.hull.points[_vert]
                d = np.sqrt((point[0] - _point[0])**2 + (point[1] - _point[1])**2 + (point[2] - _point[2])**2)
                if (largest_d is None) or (d > largest_d):
                    largest_d = d
                    largest_verts = (vert,_vert)
        r = largest_d / 2
        print(r)
        return r
    
    def getHullCOM(self):
        x = []
        y = []
        z = []
        verts = self.hull.vertices
        for vert in verts:
            point = self.hull.points[vert]
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
        avg_x = sum(x) / len(x)
        avg_y = sum(y) / len(y)
        avg_z = sum(z) / len(z)
        return (avg_x, avg_y, avg_z)
    
    def getInnerSolvent(self, r=None, threshold=0.5, write=True):
        if r is None:
            r = self.getHullRadius()
        com = self.getHullCOM()
        cx, cy, cz = com
        inner_solvent_atoms = []
        i = 0
        len_solv = len(self.solvent_coordinates.keys())
        hull_distances = []
        for point in self.protein.atom_coordinates.values():
        # for point in self.hull.points:
            hull_d = np.sqrt((point[0]-cx)**2 + (point[1]-cy)**2 + (point[2]-cz)**2)
            hull_distances.append(hull_d)
        for key,value in self.solvent_coordinates.items():
            x, y, z = value
            solvent_d = np.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
            in_boundry = 0
            for hull_d in hull_distances:
                if solvent_d <= hull_d:
                    in_boundry += 1
            if (in_boundry / len(self.protein.atom_coordinates.keys())) > threshold:
                inner_solvent_atoms.append(key)
            i += 1
            if i % 1000 == 0:
                percent = round((i/len_solv)*100, 1)
                # print('Search Completion %: {}'.format(percent))
        prev = []
        isl = []
        if len(inner_solvent_atoms) !=0:
            for mol in self.solvent_obj.molecules:
                for atom in mol.atoms:
                    if atom.id in inner_solvent_atoms:
                        if mol.index not in prev:
                            for line in mol.data:
                                isl.append(line.split())
                            prev.append(mol.index)
                        else:
                            continue
                    else:
                        continue
            if write is True:
                newfilename = self.solvent.split('.')[0] + '_inner.pdb'
                writePDB(isl, newfilename)
                print('Wrote inner solvent PDB')
        else:
            print('No atoms in threshold')
        return inner_solvent_atoms
    
    def getOuterSolvent(self, cutoff=3, remove_inner=True, inner_threshold=0.6):
        com = self.getHullCOM()
        cx, cy, cz = com
        outer_solvent_atoms = []
        for atom_coordinate in self.protein.atom_coordinates.values():
            for key, solvent_coordinate in self.solvent_coordinates.items():
                d = np.sqrt((atom_coordinate[0] - solvent_coordinate[0])**2 + (atom_coordinate[1] - solvent_coordinate[1])**2 + (atom_coordinate[2] - solvent_coordinate[2])**2)
                if d <= cutoff:
                    d_com = np.sqrt((solvent_coordinate[0] - cx)**2 + (solvent_coordinate[1] - cy)**2 + (solvent_coordinate[2] - cz)**2)
                    if d < d_com:
                        outer_solvent_atoms.append(key)
        if remove_inner is True:
            inner_solvent_atoms = self.getInnerSolvent(threshold=inner_threshold, write=False)
        prev = []
        osl = []
        if len(outer_solvent_atoms) !=0:
            for mol in self.solvent_obj.molecules:
                for atom in mol.atoms:
                    if atom.id in outer_solvent_atoms:
                        if (mol.index not in prev):
                            if (remove_inner is True) and (atom.id not in inner_solvent_atoms):
                                for line in mol.data:
                                    osl.append(line.split())
                                prev.append(mol.index)
                            if (remove_inner is True) and (atom.id in inner_solvent_atoms):
                                continue
                            if (remove_inner is False):
                                for line in mol.data:
                                    osl.append(line.split())
                        else:
                            continue
                    else:
                        continue

                        
            # print('Writing File...')
            # f = open(self.solvent, 'r')
            # contents = f.readlines()
            # f.close()
            # osl = []
            # i = 0
            # res_nums = []
            # prev = []
            # for line in contents:
            #     line_parts = line.split()
            #     atom_num = line_parts[1]
            #     if atom_num in prev:
            #         atom_num = self.solvent_atom_index[i]
            #     else:
            #         prev.append(atom_num)
            #     res_num = line_parts[4]
            #     if atom_num in outer_solvent_atoms:
            #         if remove_inner is True:
            #             if atom_num in inner_solvent_atoms:
            #                 continue
            #         else:
            #             osl.append(line_parts)
            #     if i % 1000 == 0:
            #         percent = round((i/len(contents))*100, 1)
                    # print('Writing Completion %: {}'.format(percent))
            print('Wrote outer solvent PDB')
            print(len(osl), len(outer_solvent_atoms))
            newfilename = self.solvent.split('.')[0] + '_outer.pdb'
            writePDB(osl, newfilename)
        else:
            print('No atoms in threshold')
        
    def getNumberMolecules(self):
        f = open(self.solvent, 'r')
        res_nums = []
        num_molecules = 0
        for line in f:
            line_parts = line.split()
            res_num = line_parts[4]
            if res_num not in res_nums:
                res_nums.append(res_num)
                num_molecules += 1
            else:
                continue
        f.close()
        return num_molecules


# s = SolventAnalyzer(structure='D:/Work/drude/iapp_drude/1/iapp.drude.0.pdb', trim=False)
# s.getInnerSolvent(threshold=0.3)

# s = SolventAnalyzer(structure='D:/Work/iapp/solventanalysis/myr/md_550_600.pdb', trim=True)
# s.getOuterSolvent(cutoff=3, remove_inner=False)
# filename = 'D:/Work/iapp/solventanalysis/myr/md_550_600_protein.pdb'
# editChainIDResidue(filename=filename, newfilename=filename, nterm='ACE', cterm='NH2')