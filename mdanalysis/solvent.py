import os
import sys
import mdtraj
from mdtraj.core.topology import Residue
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from pymd.mdanalysis.analysis import Analysis
from pymd.utilities.library import residues as canonical
from pymd.utilities.library import ions as _ion_names
from pymd.utilities.library import solvent as _solvent_names


class Solvent(Analysis):

    def __init__(self, inp, top, parent=None, **kwargs):
        '''
        For now, takes a dict of parameters from run.py call
        '''
        self.parent = parent
        self._inp = inp
        self._topfile = top
        # self.top = None
        self.traj = None
        self._traj = None
        self.stride = 1
        self.b = 0
        self.e = -1
        self.res = None
        self._exclude_chain = True
        self.exclude_neighbors = 2
        self._output = 'solvent.csv'
        self.job_name = 'solvent'
        self.nprocs = 'auto'
        self.compress = False
        self.selection = 'all'
        self.df = pd.DataFrame()
        self._iterload = False
        self.method = None
        self.matrix = pd.DataFrame()
        self.verbose = False
        self.__dict__.update(kwargs)
        self.atomic_masses = {
            'C':12.011,
            'A':12.011,
            'N':14.007,
            'H':1.008,
            'O':15.999,
            'S':32.06,
            'P':30.974
        }
    
    def get_residue_coms(self, frame):
        residues = {}
        for i, residue in enumerate(self.top.residues):
            if residue.name not in canonical():
                continue
            mass = 0
            x = []
            y = []
            z = []
            for atom in residue.atoms:
                atom_mass = self.atomic_masses[atom.element.symbol]
                x.append(frame._xyz[0][atom.index,0] * atom_mass)
                y.append(frame._xyz[0][atom.index,1] * atom_mass)
                z.append(frame._xyz[0][atom.index,2] * atom_mass)  
                mass  += atom_mass
            x = sum(x) / mass
            y = sum(y) / mass
            z = sum(z) / mass
            residues[residue.index] = np.array((x, y, z))
        return residues
    
    def get_residue_com(self, frame, residue):
        mass = 0
        x = []
        y = []
        z = []
        for i, atom in enumerate(residue.atoms):
            atom_mass = self.atomic_masses[atom.element.symbol]
            x.append(frame._xyz[0][atom.index,0] * atom_mass)
            y.append(frame._xyz[0][atom.index,1] * atom_mass)
            z.append(frame._xyz[0][atom.index,2] * atom_mass)  
            mass  += atom_mass
        x = sum(x) / mass
        y = sum(y) / mass
        z = sum(z) / mass
        return np.array((x,y,z))
    
    def get_frame_com(self, frame):
        x = []
        y = []
        z = []
        mass = 0
        for i, atom in enumerate(self.top.atoms):
            if atom.residue.name not in canonical():
                continue
            atom_mass = self.atomic_masses[atom.element.symbol]
            x.append(frame._xyz[0][i,0] * atom_mass)
            y.append(frame._xyz[0][i,1] * atom_mass)
            z.append(frame._xyz[0][i,2] * atom_mass)
            mass += atom_mass
        x = sum(x) / mass
        y = sum(y) / mass
        z = sum(z) / mass
        return np.array((x, y, z))

    def get_max_dist_frame(self, frame):
        '''
        get the maximum distance between any two protein atoms in a frame
        '''
        xyz = frame._xyz[0][self.protein_indeces()]
        longest = None
        for i in range(0,3):
            _min = xyz[:,i].min()
            _max = xyz[:,i].max()
            dist = _max - _min
            if longest is None:
                longest = dist
            if dist > longest:
                longest = dist
    #         print(i, dist, longest)
        return longest
    @staticmethod
    def points_within_radius(array, center_point, radius):
        # Calculate Euclidean distances between all points and the center point
        distances = np.linalg.norm(array - center_point, axis=1)

        # Find indexes of points within the specified radius
        within_radius_indices = np.where(distances <= radius)[0]

        return within_radius_indices
    
    def get_HOH_indeces(self, indeces):
        '''
        Get hydrogen indexes for oxygen atoms to make them whole
        '''
        _indeces = []
        for index in indeces:
            res = self.top.atom(index).residue
            for atom in res.atoms:
                _indeces.append(atom.index)
        return np.array(_indeces)
    
    def solvent_names(self):
        pos = ['K', 'NA', 'SOD', 'POT']
        neg = ['CL', 'CLA']
        sol_types = {
            'pos':[],
            'neg':[],
            'sol':[]
        }
        for residue in self.top.residues:
            if (residue.name in _solvent_names()) and (residue.name not in sol_types['sol']):
                sol_types['sol'].append(residue.name)
            elif (residue.name in _ion_names()):
                if (residue.name in pos) and (residue.name not in sol_types['pos']):
                    sol_types['pos'].append(residue.name)
                if (residue.name in neg) and (residue.name not in sol_types['neg']):
                    sol_types['neg'].append(residue.name)
            else:
                continue
        return sol_types
    
    def ion_names(self):
        ions = []
        for residue in self.top.residues:
            if (residue.name in _ion_names()) and (residue.name not in ions):
                ions.append(residue.name)
        return ions

    def solvent_types(self, solvent_indeces):
        return [[index, self.top.atom(index).name] for index in solvent_indeces]


    @property
    def protein_indeces(self):
        return self.top.select("not (water or name {})".format(' or name '.join(self.ion_names())))
    
    @property
    def solvent_indeces(self):
        return self.top.select('(water and name O) or name {}'.format(' or name '.join(self.ion_names())))
    
    def get_solvent_shell(self, frame, radius):
        '''
        Get solvent shell within a given radius of the protein COM at a given frame. 
        frame (traj): frame of trajectory
        radius ('auto', float): radius at which to check
        '''
        all_solv_index = []
        solvent = frame._xyz[0][self.solvent_indeces]
        for residue in frame.top.residues:
            if residue.name not in canonical():
                continue
            com = self.get_residue_com(frame, residue)
            idx = self.points_within_radius(solvent, com, radius)
            solvent_within_radius = self.solvent_indeces[idx]
            all_solv_index.append(solvent_within_radius[:])
        all_solv_concat = np.concatenate(all_solv_index)
        all_solv_unique = np.unique(all_solv_concat)
        all_solvent = self.get_HOH_indeces(all_solv_unique)
        return all_solvent
    

    def run(self, method='shell', radius=1, stride=1, selection='all', chunk=100, output_prefix='solvent', output_path=''):
        self._output = f'{output_prefix}.csv'
        if self.traj_iter is None:
            self.iterloadTrajectory(stride=stride, selection=selection, chunk=chunk)
        frame_idx = 0
        solvent_data = []
        solvent_ndx = []

        for chunk in self.traj_iter: # type: ignore
            first_time = None
            for frame in chunk:
                shell = self.get_solvent_shell(frame, radius)
                solvent_ndx.append(shell)
                solvent_data.append([frame_idx, frame._time[0], len(shell)])
                frame_idx += 1
                if first_time is None:
                    first_time = frame._time[0]
            np.save(os.path.join(output_path, 'solvidx.{}.{}.npy'.format(str(first_time), str(frame._time[0]))), solvent_ndx)
        df = pd.DataFrame(solvent_data, columns=['frame_index', 'time', 'n_solvent'])
        df.to_csv(os.path.join(output_path, '{}.csv'.format(output_prefix)))
                






    



