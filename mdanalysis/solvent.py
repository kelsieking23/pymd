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
        self.traj_iter = None
    
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
    
    def get_protein_com(self, frame):
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

    def get_chain_com(self, frame, chain_index):
        mass = 0
        x = []
        y = []
        z = []
        for i, atom in enumerate(self.top.atoms):
            if atom.residue.chain.index != chain_index:
                continue
            atom_mass = self.atomic_masses[atom.element.symbol]
            x.append(frame._xyz[0][atom.index,0] * atom_mass)
            y.append(frame._xyz[0][atom.index,1] * atom_mass)
            z.append(frame._xyz[0][atom.index,2] * atom_mass)  
            mass  += atom_mass
        x = sum(x) / mass
        y = sum(y) / mass
        z = sum(z) / mass
        return np.array((x,y,z))

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

    def get_max_dist(self, frame, atoms):
        '''
        get max dist between any two atoms
        '''
        xyz = frame._xyz[0][np.array([atom.index for atom in atoms])]
        longest = None
        for i in range(0,3):
            _min = xyz[:,i].min()
            _max = xyz[:,i].max()
            dist = _max - _min
            if longest is None:
                longest = dist
            if dist > longest:
                longest = dist
        return longest
    
    def get_com(self, frame, atoms):
        '''
        Get COM from a collection of atoms
        '''
        mass = 0
        x = []
        y = []
        z = []
        for i, atom in enumerate(atoms):
            atom_mass = self.atomic_masses[atom.element.symbol]
            x.append(frame._xyz[0][atom.index,0] * atom_mass)
            y.append(frame._xyz[0][atom.index,1] * atom_mass)
            z.append(frame._xyz[0][atom.index,2] * atom_mass)  
            mass  += atom_mass
        x = sum(x) / mass
        y = sum(y) / mass
        z = sum(z) / mass
        return np.array((x,y,z))
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
    
    def trim_solvent(self, frame, by='protein'):
        '''
        Get solvent indeces within a certain radius of the protein, or by individual peptides
        returns (trimmed_index, trimmed_xyz)
        '''
        if self.verbose:
            print('>>> Trimming solvent indeces by {}'.format(by))
            print('***')
        element_coms = []
        radii = []
        if by == 'protein':
            element_coms = [self.get_protein_com(frame)]
            radii = [self.get_max_dist_frame(frame)]
        elif (by == 'peptide') or (by == 'chain'):
            for chain_index in self.chain_idx:
                if self.verbose:
                    print('>>>> Chain {}'.format(chain_index))
                atoms = [atom for atom in self.protein_atoms if atom.residue.chain.index == chain_index]
                if self.verbose:
                    print('>>>> Atoms in chain: {}'.format(len(atoms)))
                com = self.get_com(frame, atoms)
                if self.verbose:
                    print('>>>> {}'.format(com))
                element_coms.append(com)
                radius = self.get_max_dist(frame, atoms)
                if radius is None:
                    raise ValueError('Radius for chain index {} was None: is selection empty? Num atoms in selection: {}'.format(chain_index, len(atoms)))
                radii.append(radius + 0.5)
                if self.verbose:
                    print('>>>> Radius {:.3f}, using 0.5 nm buffer ({:.3f})'.format(radius, radius+0.5))
        else:
            raise ValueError('Can only trim solvent by protein, or peptide/chain. Keywords "protein", "peptide", "chain" only accepted')
        if self.verbose:
            print('>>>> Distance calculations for solvent trimming ...')
        solvent = frame._xyz[0][self.solvent_indeces]
        all_solv_index = []
        for com, radius in zip(element_coms, radii):
            idx = self.points_within_radius(solvent, com, radius)
            solvent_within_radius = self.solvent_indeces[idx]
            all_solv_index.append(solvent_within_radius[:])
        if self.verbose:
            print('>>>> Preparing trimmed solvent ...')
        all_solv_concat = np.concatenate(all_solv_index)
        trimmed_index = np.unique(all_solv_concat)
        trimmed_xyz = solvent[trimmed_index]
        if self.verbose:
            print('>>>> Solvent preparation complete.')
            print('>>>> Trimmed solvent has {} atoms'.format(len(trimmed_index)))
            print('***')
        return trimmed_index, trimmed_xyz

    def get_solvent_shell(self, frame, radius, trim_by):
        '''
        Get solvent shell within a given radius of the protein COM at a given frame. 
        frame (traj): frame of trajectory
        radius ('auto', float): radius at which to check
        '''
        _canonical = canonical()
        all_solv_index = []
        if trim_by is None:
            solvent_idx = self.solvent_indeces
            solvent_xyz = frame._xyz[0][self.solvent_indeces]
            if self.verbose:
                print('>>> No solvent trimming. Solvent has {} atoms'.format(len(solvent_idx)))
        else:
            solvent_idx, solvent_xyz = self.trim_solvent(frame, by=trim_by)
        if self.verbose:
            print('>>> Calculating distances...')
        for residue in frame.top.residues:
            if residue.name not in _canonical:
                continue
            if self.verbose:
                print(f'>>> (Frame {frame._time[0]}) {residue.name}')
                print('>>> Getting residue COM ...')
            com = self.get_residue_com(frame, residue)
            if self.verbose:
                print(f'>>> {com}')
                print('>>> Getting solvent indeces within radius of {} ...'.format(radius))
            idx = self.points_within_radius(solvent_xyz, com, radius)
            solvent_within_radius = solvent_idx[idx]
            all_solv_index.append(solvent_within_radius[:])
        if self.verbose:
            print('>>> Concatenating solvent indeces ...')
        all_solv_concat = np.concatenate(all_solv_index)
        if self.verbose:
            print('>>> Finding unique solvent indeces ...')
        all_solv_unique = np.unique(all_solv_concat)
        if self.verbose:
            print('>>> Making solvent whole (finding hydrogens) ...')
        all_solvent = self.get_HOH_indeces(all_solv_unique)
        if self.verbose:
            print(f'    Frame {frame._time[0]} complete.')
            print('>>> Found {} unique solvent atoms'.format(len(all_solv_unique)))
            print('>>> Found {} total solvent atoms (including hydrogens)'.format(len(all_solvent)))
        return all_solvent
    

    def run(self, method='shell', radius=1.0, stride=1, selection='all', chunk=100, output_prefix='solvent', trim_by='protein', output_path=''):
        self._output = f'{output_prefix}.csv'
        if self.traj_iter is None:
            self.iterloadTrajectory(stride=stride, selection=selection, chunk=chunk)
        frame_idx = 0
        solvent_data = []
        solvent_ndx = []
        chunk_idx = 1
        for chunk in self.traj_iter: # type: ignore
            if self.verbose:
                print('> Chunk {}'.format(chunk_idx))
                print('> ', chunk)
                print('***')
            first_time = None
            for frame in chunk:
                if self.verbose:
                    print('>> Frame index {}, time {} ps'.format(frame_idx, frame._time[0]))
                    print('***')
                shell = self.get_solvent_shell(frame, radius, trim_by)
                solvent_ndx.append(shell)
                solvent_data.append([frame_idx, frame._time[0], len(shell)])
                frame_idx += 1
                if first_time is None:
                    first_time = frame._time[0]
            np.save(os.path.join(output_path, 'solvidx.{}.{}.npy'.format(str(first_time), str(frame._time[0]))), solvent_ndx)
            if self.verbose:
                print('Wrote {} '.format('solvidx.{}.{}.npy'.format(str(first_time), str(frame._time[0]))))
        df = pd.DataFrame(solvent_data, columns=['frame_index', 'time', 'n_solvent'])
        df.to_csv(os.path.join(output_path, '{}.csv'.format(output_prefix)))
                






    



