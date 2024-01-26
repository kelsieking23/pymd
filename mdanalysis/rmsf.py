import os
import sys
import mdtraj
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from pymd.mdanalysis.analysis import Analysis
from pymd.mdanalysis.rmsd import RMSD
from pymd.utilities.library import residues as canonical
from pymd.utilities.library import ions as _ion_names
from pymd.utilities.library import solvent as _solvent_names
from pymd.plot.plot_v2 import Plotter
import time


class RMSF(Analysis):

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
        self._output = 'rmsf.csv'
        self.job_name = 'rmsf'
        self.method = None
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


    def run(self, output='rmsf.csv', selection='sidechain', average_residues=True, ref_idx=0, **kwargs):
        self._output = output
        if average_residues:
            self.method = 'residue'
        else:
            self.method='atom'
        self.save(ref_idx=ref_idx, **kwargs)
        if self.verbose:
            print('Starting RMSF calculation...')
        ref = self.traj[ref_idx]
        rms = self.rmsf(ref, **kwargs)
        if self.verbose:
            print('Calculation complete.')
            print('Parsing selection:')
            print('   * Selection: {}'.format(selection))
            print('   * Averaging residues: {}'.format(str(average_residues)))
        self.df = self.makeDataFrame(rms, selection, average_residues)
        if self.verbose:
            print('Parsing complete.')
            print('Data Preview (first 10 rows):')
            print(self.df.head(10))
        self.df.to_csv(self.output)
        if self.verbose:
            print('Wrote {}'.format(self.output))
        return self.df

    def rmsf(self, ref, **kwargs):
        self.traj.center_coordinates()
        ref.center_coordinates()
        rms = mdtraj.rmsf(self.traj, ref, precentered=True, **kwargs)
        return rms
    
    def makeDataFrame(self, rms, selection, average_residues=True):
        df = self.top.to_dataframe()[0]
        df['rms'] = rms
        df = df.loc[self.top.select(selection),:].reset_index(drop=True)
        if average_residues:
            resi = [] 
            for residue in self.top.residues:
                data = df[(df['resName'] == residue.name) & (df['resSeq'] == residue.resSeq) & (df['chainID'] == residue.chain.index)]
                resi.append([residue.name, residue.resSeq, residue.name + str(residue.resSeq), residue.chain.index, selection, data['rms'].mean()])
            rmsf = pd.DataFrame(resi, columns=['res_name', 'res_num', 'res_id', 'chain', 'selection', 'rms'])
            return rmsf
        else:
            df['res_id'] = [name + str(seq) for name, seq in zip(df['resName'], df['resSeq'])]
            df['selection'] = [selection] * len(df)
            rmsf = pd.DataFrame([df['serial'], df['name'], df['element'], df['resName'],
                                 df['resSeq'], df['res_id'], df['chainID'], df['selection'],
                                 df['rms']]).T
            rmsf.columns = ['serial', 'atom_name', 'elem', 'res_name', 'res_num', 'res_id', 'chain',
                            'selection', 'rms']
            return rmsf

    def plot(self, output='rmsf.png', average_by=None, plot_by=None, panel=False, nrows=1, ncols=1, w=8, h=6, suptitle='RMSF', titles=[], ax=None, show=False, **kwargs):
        if average_by is not None:
            self.plot_average(output, average_by, show, w, h, **kwargs)
        elif plot_by is not None:
            self.plot_by(output, plot_by, show, panel, nrows, ncols, w, h, titles, **kwargs)
        else:
            self.plot_single(output, show, w, h, **kwargs)

    def plot_average(self, output, average_by, show, w, h, **kwargs):
        plot_df = pd.DataFrame()
        df = self.df.select_dtypes(exclude=['object'])
        df[average_by] = self.df[average_by]
        plot_df['mean'] = df.groupby([average_by]).mean().sort_values(by=average_by, ascending=True)['rms'].values
        plot_df['std'] = df.groupby([average_by]).std().sort_values(by=average_by, ascending=True)['rms'].values
        if 'std' not in kwargs.keys():
            kwargs['std'] = True
        idx = []
        for i, item in enumerate(df.groupby([average_by]).mean().index):
            try:
                idx.append(int(item))
            except:
                idx.append(int(i+1))

        if 'x_label' not in kwargs.keys():
            if self.method == 'atom':
                kwargs['x_label'] = 'Atom Number'
            else:
                kwargs['x_label'] = 'Residue Number'
        if 'y_label' not in kwargs.keys():
            kwargs['y_label'] = 'RMSF (nm)'
        out_path = os.path.join(self.root, output)
        plotter = Plotter(w=w, h=h)
        plotter.timeseries(plot_df, out=out_path, show=show, **kwargs)


    def plot_by(self, output, plot_by, show, panel, nrows, ncols, w, h, titles, **kwargs):
        #TODO: make a panel plot thingy
        pass

    def plot_single(self, output, show, w, h, **kwargs):
        plot_df = pd.DataFrame()
        plot_df['rms'] = self.df['rms']
        if self.method == 'atom':
            plot_df.index = self.df['serial']
        else:
            plot_df.index = self.df['res_num']
        plotter = Plotter(w=w, h=h)
        out_path = os.path.join(self.root, '{}.png'.format(output))
        plotter.timeseries(plot_df, out=out_path, show=show, **kwargs)