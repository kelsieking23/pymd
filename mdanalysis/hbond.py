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
import multiprocessing as mp


class Hbond(Analysis):

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
        self._output = 'hbond.csv'
        self.job_name = 'hbond'
        self.nprocs = 'auto'
        self.compress = False
        self.selection = 'all'
        self.df = pd.DataFrame()
        self._iterload = False
        self.method = None
        self.matrix = pd.DataFrame()
        self.verbose = False
        self.__dict__.update(kwargs)

    def starmap_wrapper(self, arg):
        args, kwargs = arg
        return self.runBakerHubbard(*args, **kwargs)
    
    def launch_mp_bakerHubbard(self, processes, freq_step, freq_max, selection, **kwargs):
        jobs = []
        pool = mp.Pool(processes)
        for i in np.arange(freq_step, freq_max, freq_step):
            jobs.append((i, selection, i+i))
        arg = [(j, kwargs) for j in jobs]
        results = pool.map(self.starmap_wrapper, arg)
        return results
    
    def run(self, selection, freq_step=0.01, freq_max=1, method='baker-hubbard', output='hbond.csv', processes=1, **kwargs):
        self._output = output
        try:
            assert self.traj is not None
        except:
            raise ValueError('No trajectory loaded! Load with hbond.loadTrajectory()')
        if self.verbose:
            print('Writing job data...')
        if isinstance(selection, str):
            selection_tosave = selection
        elif isinstance(selection, (list, tuple, np.ndarray,)):
            selection_tosave = '[list of length {}]'.format(len(selection))
        else:
            raise ValueError('Selection must be a string (mdtraj selection string) or array-like. Selection was {}'.format(selection))
        self.save(selection=selection_tosave, method=method, freq_step=freq_step)
        if self.verbose:
            print('Job Data Written')
            print('Method: {}'.format(method))
        if self.verbose:
            print('Beggining hbond calculation...')
        if (method == 'baker-hubbard'):
            self.df = self.runBakerHubbard(freq_step, selection, freq_max=freq_max, **kwargs)
        else:
            raise ValueError('No method called {}. Only supported method currently is baker-hubbard'.format(method))
        if self.verbose:
            print('Hbond calculation complete.')
        self.df.to_csv(self.output)
        
    def getBlankMatrix(self, top_atoms):
        m = {}
        res_ids = []
        for residue in self.top.residues:
            res_id = residue.name + str(residue.resSeq)
            if res_id not in res_ids:
                res_ids.append(res_id)
        for atom in top_atoms:
            res_id = self.top.atom(atom).residue.name + str(self.top.atom(atom).residue.resSeq)
            if res_id not in res_ids:
                res_ids.append(res_id)
        for res_id in res_ids:
            m[res_id] = {}
            for res in res_ids:
                m[res_id][res] = 0 
        return m
    
    def runBakerHubbard(self, freq_step, selection, freq_max, **kwargs):
        if isinstance(selection, str):
            top_atoms = self.select(selection)
        elif isinstance(selection, (list, tuple, np.ndarray,)):
            top_atoms = selection
        else:
            raise ValueError('Selection must be a string (mdtraj selection string) or array-like. Selection was {}'.format(selection))
        hbond_data = []
        base, ext = tuple(os.path.splitext(self.output))
        freqs = np.arange(freq_step, freq_max, freq_step)
        percents = np.arange(0.1, 1.1, 0.1)
        for i, freq in enumerate(freqs):
            if freq > 1:
                break
            hbond = self.bakerHubbard(freq, **kwargs)
            hbond_pairs_output = base + '.pairs_freq{:1.2f}{}'.format(freq, ext)
            pd.DataFrame(hbond).to_csv(hbond_pairs_output)
            data = self.getBlankMatrix(top_atoms)
            for bond in hbond:
                d = bond[0]
                a = bond[2]
                if (d in top_atoms) and (a in top_atoms):
                    dres = self.top.atom(d).residue.name + str(self.top.atom(d).residue.resSeq)
                    ares = self.top.atom(a).residue.name + str(self.top.atom(a).residue.resSeq)
                    data[dres][ares] = freq
                    data[ares][dres] = freq
            hbond_data.append(pd.DataFrame(data).astype(float))
            hb_data_output = base + '.matrix_freq{:1.2f}{}'.format(freq, ext)
            hbond_data[-1].to_csv(hb_data_output)
            if self.verbose:
                print('({}/{}) complete'.format(i+1, len(freqs)))
                percent = i/len(freqs)
                if percent in percents:
                    print('{:.0f}% complete...'.format(percent*100))
        if self.verbose:
            print('100% complete.')
            print('Finding max fraction for possible pairs...')
        maxed = pd.concat(hbond_data).groupby(level=0, as_index=False).max()
        maxed.index = maxed.columns
        return maxed

    def run_wernet_nilsson(self, selection, **kwargs):
        pass

    
    def bakerHubbard(self, freq, **kwargs):
        hbonds = mdtraj.baker_hubbard(self._traj, freq=freq, **kwargs)
        return hbonds
    


