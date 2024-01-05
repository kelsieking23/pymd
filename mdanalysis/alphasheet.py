import json
import multiprocessing as mp
import os
import psutil
import sys

import mdtraj
import numpy as np
import pandas as pd
from statistics import mean
from sklearn.metrics.pairwise import euclidean_distances
from collections.abc import Iterable
from sklearn.preprocessing import StandardScaler
from pymd.mdanalysis.analysis import Analysis
from pymd.structure.peptide import Peptide
from pymd.plot.plot_v2 import Plotter

class AlphaSheet(Analysis):

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
        self._output = 'alphasheet.csv'
        self.job_name = 'alphasheet'
        self.nprocs = 'auto'
        self.compress = False
        self.selection = 'all'
        self.df = pd.DataFrame()
        self._iterload = False
        self.method = None
        self.matrix = pd.DataFrame()
        self.verbose = False
        self.__dict__.update(kwargs)


    def get_angle_data(self, atype, traj=None):
        if traj is None:
            if self.traj is not None:
                traj = self.traj
            else:
                raise ValueError('Trajectory not defined. Load one or pass into function. ')
        top = traj.topology
        if atype == 'phi':
            index, angles = mdtraj.compute_phi(traj)
        elif atype == 'psi':
            index, angles = mdtraj.compute_psi(traj)
        else:
            raise ValueError('atype can only be "phi" or "psi" - {} not valid'.format(atype))
        all_frames = []
        for time_idx, frame_angles in enumerate(angles):
            data = pd.DataFrame()
            time = []
            res_id = []
            res_seq = []
            res_index = []
            chain_index = []
            for group in index:
                if atype == 'phi':
                    res = top.atom(group[1]).residue
                if atype == 'psi':
                    res = top.atom(group[0]).residue
                res_id.append('{}{}'.format(res.name, str(res.resSeq)))
                res_seq.append(res.resSeq)
                res_index.append(res.index)
                chain_index.append(res.chain.index)
                time.append(traj._time[time_idx])
            data['time'] = time
            data['res_id'] = res_id
            data['res_seq'] = res_seq
            data['res_index'] = res_index
            data['chain_index'] = chain_index
            data['angles'] = angles[time_idx]
            data['angles_rad'] = data['angles'] * 180 / np.pi
            if atype == 'phi':
                data = data.drop(data[data['res_seq'] == np.max(np.array(res_seq))].index).reset_index(drop=True)
            if atype == 'psi':
                data = data.drop(data[data['res_seq'] == np.min(np.array(res_seq))].index).reset_index(drop=True)
            all_frames.append(data)
        data = pd.concat(all_frames).reset_index(drop=True)
        return data
    
    def get_phi(self, traj=None):
        return self.get_angle_data('phi', traj)
    
    def get_psi(self, traj=None):
        return self.get_angle_data('psi', traj)
    
    def get_phi_psi(self, traj=None):
        phi = self.get_phi(traj)
        psi = self.get_psi(traj)
        return phi, psi
    
    def classify_ramachandran(self, phi, psi):
        data = pd.DataFrame()
        data['time'] = phi['time']
        data['res_id'] = phi['res_id']
        data['res_seq'] = phi['res_seq']
        data['chain_index'] = phi['chain_index']
        data['phi'] = phi['angles_rad']
        data['psi'] = psi['angles_rad']
        classification = []
        for phi, psi in zip(data['phi'], data['psi']):
            if (0 > phi) and (psi < 45):
                classification.append('left-helix')
            elif (phi > 0) and (psi > -90):
                classification.append('right-helix')
            else:
                classification.append('other')
        data['class'] = classification
        return data

    def find_alphasheets(self, df):
        groups = []
        # Create a mask for rows with 'R'
        r_mask = df['class'] == 'right-helix'
        # Shift the 'class' column to check the rows above and below
        prev_class = df['class'].shift(1)
        next_class = df['class'].shift(-1)
        # Filter rows where the current row has 'R' and both above and below rows have 'L'
        desired_rows_r = df[r_mask & (prev_class == 'left-helix') & (next_class == 'left-helix')]
        for index in desired_rows_r.index:
            groups.append(df.loc[[index-1, index, index+1],:])
        # Create a mask for rows with 'L'
        l_mask = df['class'] == 'left-helix'
        # Shift the 'class' column to check the rows above and below
        prev_class = df['class'].shift(1)
        next_class = df['class'].shift(-1)
        # Filter rows where the current row has 'L' and both above and below rows have 'R'
        desired_rows_l = df[l_mask & (prev_class == 'right-helix') & (next_class == 'right-helix')]
        for index in desired_rows_l.index:
            groups.append(df.loc[[index-1, index, index+1],:])
        allgroups = pd.concat(groups).drop_duplicates()
        self.df = allgroups
        if self.verbose:
            print('Writing {}'.format(self.output))
        self.df.to_csv(self.output)
        return self.df
    
    def run(self, output='alphasheet.csv', traj=None):
        self._output = output
        if self.verbose:
            print('Writing job data...')
        self.save()
        if self.verbose:
            print('Computing phi/psi angles...')
        phi, psi = self.get_phi_psi(traj)
        if self.verbose:
            print('Classifying left- and right-handed helices...')
        df = self.classify_ramachandran(phi, psi)
        if self.verbose:
            print('Finding contiguous alpha-sheets...')
        self.find_alphasheets(df)
        if self.verbose:
            print('Job complete.')
        return self.df

    def get_plot_data(self, aggregate_chains=False):
        counts = {}
        for residue in self.top.residues:
            if aggregate_chains:
                if residue.chain.index > 0:
                    break
            res_id = residue.name + str(residue.resSeq)
            key = res_id
            rows = self.df[(self.df['res_id'] == res_id)]
            if not aggregate_chains:
                rows = rows[rows['chain_index' == residue.chain.index]]
                key = res_id + '_' + str(residue.chain.index)
            counts[key] = len(rows)
        return pd.Series(data=counts.values(), index=counts.keys())
            
    def plot(self, output='alphasheet.png', show=False, aggregate_chains=False, w=10, h=6, 
             x_label='Residue', y_label=r'$\alpha$-Strand Count (Frames)', **kwargs):
        plotter = Plotter(w=w, h=h)
        df = self.get_plot_data(aggregate_chains)
        if 'tick_label' not in kwargs.keys():
            kwargs['tick_label'] = df.index
        plotter.bar(df, show=show, panel=False, output=output, x='index', x_label=x_label, y_label=y_label, **kwargs)

