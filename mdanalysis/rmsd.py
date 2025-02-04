import os
import sys
import mdtraj
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from pymd.mdanalysis.analysis import Analysis
from pymd.utilities.library import residues as canonical
from pymd.utilities.library import ions as _ion_names
from pymd.utilities.library import solvent as _solvent_names
from pymd.plot.plot_v2 import Plotter
import time


class RMSD(Analysis):

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
        self._output = 'rmsd.csv'
        self.job_name = 'rmsd'
        self.nprocs = 'auto'
        self.compress = False
        self.selection = 'all'
        self.df = pd.DataFrame()
        self._iterload = False
        self.method = None
        self.matrix = pd.DataFrame()
        self.verbose = False
        self._top = mdtraj.load(self.topfile).topology
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

    @staticmethod
    def rmsd(frames, ref):
        rms = []
        for frame in frames:
            rms.append(np.sqrt((np.linalg.norm(frame-ref) / frame.shape[0])**2))
        return np.array(rms)

    
    def run(self, selection='backbone', by='protein', ref_idx=0, output='rmsd.csv'):
        self._output = output
        self.save(selection=selection, by=by)
        if self.verbose:
            if isinstance(by, str):
                print('Starting RMSD calculation, by ({}) using selection ({})'.format(by, selection))
            else:
                print('Starting RMSD calculation, by ({}) using selection(s) ({})'.format(by, ', '.join(selection)))
        selections = self.by_selection(by, selection)
        df = pd.DataFrame()
        for selstr, trj, ref in selections:
            frames = trj._xyz
            if self.verbose:
                print('Running for: {}'.format(selstr))
            ref = ref._xyz
            df_column = '_'.join(selstr.split())
            df[df_column] = self.rmsd(frames, ref)
        self.df = df
        if not self.df.empty:
            self.df.index = self.traj.time
        self.df.to_csv(self.output)
        if self.verbose:
            print('RMSD calculation complete.')

    def by_selection(self, by, selection):
        if not isinstance(self.traj, mdtraj.Trajectory):
            raise TypeError('No trajectory loaded')
        selections = []
        if by.lower() == 'protein':
            sele = self.top.select('protein and {}'.format(selection))
            return [('protein_and_{}'.format(selection), self.traj.atom_slice(sele)._xyz)]
        if isinstance(by, str):
            if by.lower() in ('chain', 'peptide', 'chainid'):
                for chain in self.top.chains:
                    selstr = '(chainid {}) and {}'.format(chain.index, selection)
                    sele = self.top.select(selstr)
                    selections.append((selstr, self.traj.atom_slice(sele)._xyz))
            else:
                selstr = '({}) and {}'.format(by, selection)
                sele = self.top.select(selstr)
                selections.append((selstr, self.traj.atom_slice(sele)._xyz))
        else:
            for item in by:
                selstr = '({}) and {}'.format(item, selection)
                sele = self.top.select(selstr)
                selections.append((selstr, self.traj.atom_slice(sele)._xyz))
        return selections
    
    def plot(self, out='rmsd.png', panel=False, ax=None, show=False, titles=[], units='ns', **kwargs):
        plotter = Plotter()
        if units == 'ns':
            self.df.index = [i/1000 for i in self.df.index]
        if 'x_label' not in kwargs.keys():
            kwargs['x_label'] = 'Time ({})'.format(units)
        if 'y_label' not in kwargs.keys():
            kwargs['y_label'] = 'RMSD (nm)'
        out_path = os.path.join(self.root, out)
        plotter.timeseries(self.df, out=out_path, panel=panel, ax=ax, show=show, titles=titles, **kwargs)