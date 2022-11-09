import json
import multiprocessing as mp
import os
import psutil
import sys

import mdtraj
import numpy as np
import pandas as pd
from statistics import mean
from collections.abc import Iterable
from sklearn.preprocessing import StandardScaler

class Distance:

    def __init__(self, inp, top, parent=None, **kwargs):
        '''
        For now, takes a dict of parameters from run.py call
        '''
        self.parent = parent
        self._inp = inp
        self._topfile = top
        self.top = None
        self.traj = None
        self._traj = None
        self.stride = 1
        self.b = 0
        self.e = -1
        self.res = None
        self._exclude_chain = True
        self.exclude_neighbors = 2
        self._output = None
        self.job_name = 'dist'
        self.nprocs = 'auto'
        self.compress = False
        self.selection = 'all'
        self.df = pd.DataFrame()
        self._iterload = False
        self.matrix = pd.DataFrame()
        self.__dict__.update(kwargs)
        if self._iterload:
            self.iterload(self.stride)

    def save(self):
        params = {}
        manual_keys = ['parent', 'df', 'matrix']
        for key, value in self.__dict__.items():
            if key not in manual_keys:
                params[key] = value
        filename = os.path.join(self.root, 'job_params.json')
        with open(filename, 'w') as f:
            params_dict = json.dumps(params)
            f.write(params_dict)

    @classmethod
    def from_json(cls, path, inp=None, top=None, parent=None):
        with open(os.path.join(path, 'job_params.json'), 'r') as f:
            params = json.load(f)
        dic = {}
        manual_keys = ['inp', 'top', 'parent', 'df', 'matrix']
        for key, value in params.items():
            if key not in manual_keys:
                dic[key] = value
        for filename in os.listdir(path):
            if filename.endswith('csv'):
                dic['df'] = pd.read_csv(os.path.join(path, filename), index_col=0)
        if parent is not None:
            inp = parent.inp
            top = parent.topfile
        return cls(inp, top, parent, **dic)

    @property
    def inp(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self._inp)  # type: ignore
        return self._inp

    @property
    def topfile(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self._topfile)  # type: ignore
        return self._topfile
    
    @property
    def output(self):
        if self.parent is not None:
            if not os.path.isdir(os.path.join(self.parent.root, self.job_name)): # type: ignore
                os.mkdir(os.path.join(self.parent.root, self.job_name)) # type: ignore
            if self.compress:
                self._output = os.path.splitext(self._output)[0] + '.npy' # type: ignore
            return os.path.join(self.parent.root, self.job_name, self._output) # type: ignore
        return self._output
    
    @property
    def root(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self.job_name)
        else:
            return os.path.join(os.getcwd(), self.job_name)
    
    def load(self, stride=None):
        if stride is not None:
            if stride != 0:
                traj = mdtraj.load(self.inp, top=self.topfile, stride=stride)
            else:
                traj = mdtraj.load(self.inp, top=self.topfile)
        else:
            if self.stride != 0:
                traj = mdtraj.load(self.inp, top=self.topfile, stride=self.stride)
            else:
                traj = mdtraj.load(self.inp, top=self.topfile)
        traj = traj.superpose(traj)
        self._traj = traj
        if (self.e == -1):
            self.traj = traj.center_coordinates()[self.b:]
        else:
            self.traj = traj.center_coordinates()[self.b:self.e]
        if self.selection != 'all':
            sele = self.traj.top.select(self.selection)
            self.traj = self.traj.atom_slice(sele)
        self.frames = self.traj._xyz
        self.top = self.traj.top
        return self

    def iterload(self, stride=None):
        if stride is not None:
            if stride != 0:
                traj = mdtraj.iterload(self.inp, top=self.topfile, stride=stride)
            else:
                traj = mdtraj.iterload(self.inp, top=self.topfile)
        else:
            if self.stride != 0:
                traj = mdtraj.iterload(self.inp, top=self.topfile, stride=self.stride)
            else:
                traj = mdtraj.iterload(self.inp, top=self.topfile)
        self.traj = traj
        self._traj = traj
        if self.selection != 'all':
            sele = self.traj.top.select(self.selection)
            self.traj = self.traj.atom_slice(sele)
        self.top = self.traj.top
        return self

    def getPartitions(self):
        if self.nprocs == 'auto':
            nprocs = int(mp.cpu_count() // 2)
        else:
            nprocs = self.nprocs
        if self.traj is None:
            self.load()
        nframes, _, _ = self.traj._xyz.shape # type: ignore
        interval = int(nframes // nprocs)
        partitions = []
        procid=1
        for i in range(0, nframes, interval):
            data = {
                'b':i,
                'e':i+interval,
                'procid':procid,
            }
            partitions.append(data)
            procid+=1
            if ((i + interval + interval) > nframes) and (i+interval != nframes):
                data = {
                    'b':i+interval,
                    'e':nframes,
                    'procid':procid
                }
                partitions.append(data)
                break
        return partitions, nprocs

    @staticmethod
    def atom_pairs(sele_indices, excluded_neighbors = 2):
        p = []
        for i in range(len(sele_indices)):
            for k in range(i + 1, len(sele_indices)):
                I = sele_indices[i]
                K = sele_indices[k]
                if (I > K):
                    I = sele_indices[k]
                    K = sele_indices[i]
                if (K > I + excluded_neighbors):
                    p.append([I, K])
        return np.array(p)

    @staticmethod
    def residue_pairs(sele_indices, excluded_neighbors = 2):
        p = []
        for i in range(len(sele_indices)):
            for k in range(i + 1, len(sele_indices)):
                I = sele_indices[i]
                K = sele_indices[k]
                if (I > K):
                    I = sele_indices[k]
                    K = sele_indices[i]
                if (K > I + excluded_neighbors):
                    p.append([I, K])
        return np.array(p)

    def exclude_chains(self, pairs=None):
        p = []
        for pair in pairs:
            i = pair[0]
            k = pair[1]
            if self.top.residue(i).chain == self.top.residue(k).chain:
                pass
            else:
                p.append([i,k])
        return np.array(p)

    def run(self, method='residue', output='dist.csv', **kwargs):
        self.save()
        self.load(self.stride)
        if method == 'atom':
            df = self.by_atom()
            return df
        if method == 'residue':
            # pairs = self.pairs()
            pairs = self.residue_pairs([i for i in range(self.traj.n_residues)])
            if self._exclude_chain:
                pairs = self.exclude_chains(pairs)
            df = self.by_residue(contact_pairs=pairs, **kwargs)
            return df
    
    def by_atom(self):
        sele = self.top.select(self.selection)
        pairs = self.atom_pairs(sele, excluded_neighbors=2)
        distances = mdtraj.compute_distances(self.traj, atom_pairs=pairs)
        labels = ['_'.join(list(map(str, x))) for x in pairs]
        self.df = pd.DataFrame(distances, columns=labels)
        if self.output is not None:
            print('Writing {}'.format(self.output))
            print('Shape: {}'.format(self.df.shape))
            np.save(self.output, distances)
        return self.df
    
    def by_residue(self, contact_pairs, **kwargs):
        if self.parent is not None:
            print('Computing residue contacts for {}...'.format(self.parent.id))
        distances, pairs = mdtraj.compute_contacts(self.traj, contacts=contact_pairs, **kwargs) # type: ignore
        distances = pd.DataFrame(distances)
        pairs = pd.DataFrame(pairs)
        distances.columns = ['{}_{}'.format(pairs.loc[index,0], pairs.loc[index, 1]) for index in pairs.index]  # type: ignore
        self.df = distances
        if self.output is not None:
            print('Writing {}'.format(self.output))
            print('Shape: {}'.format(self.df.shape))
            self.df.to_csv(self.output)
        return self.df

    def to_matrix(self, interval=(0,-1)):
        df = self.df
        b, e = interval
        matrix = {}
        if (self.top is None) and (self.parent is not None):
            self.top = self.parent.top
        # init matrix
        n_residues = int(self.top.n_residues / self.top.n_chains) # type:ignore
        matrix = {}
        for i in range(1, n_residues+1):
            matrix[i] = {}
            for j in range(1, n_residues + 1):
                matrix[i][j] = []
        for column in df.columns:
            i,j = list(map(int, column.split('_')))
            iresnum = self.top.residue(i).resSeq # type:ignore
            jresnum = self.top.residue(j).resSeq # type:ignore
            if (e == -1):
                matrix[iresnum][jresnum].append(df.loc[:, column].mean())
                matrix[jresnum][iresnum].append(df.loc[:, column].mean())
            else:
                matrix[iresnum][jresnum].append(df.loc[b:e, column].mean())
                matrix[jresnum][iresnum].append(df.loc[b:e, column].mean())
        for i in range(1, n_residues+1):
            for j in range(1, n_residues + 1):
                if isinstance(matrix[i][j], Iterable):
                    matrix[i][j] = mean(matrix[i][j])
                if isinstance(matrix[i][j], Iterable):
                    matrix[j][i] = mean(matrix[j][i])
        df = pd.DataFrame(matrix)
        df.attrs['ptype'] = 'heatmap'
        df.attrs['interval'] = (b, e)
        if self.parent is not None:
            self.parent.df = df
        self.matrix = df
        return df
    
    def to_matrices(self, interval=200):
        dfs = []
        for i in range(0, self.df.index[-1], interval): # type:ignore
            df = self.to_matrix((i, i+interval))
            dfs.append(df)
        return dfs

    def normalize(self, method='by_residue_atoms'):
        if method == 'sklearn':
            scaler = StandardScaler()
            scaler.fit(self.matrix)
            data = scaler.transform(self.matrix)
            attrs = self.matrix.attrs
            self.matrix = pd.DataFrame(data, index=self.matrix.index, columns=self.matrix.columns)
            self.matrix.attrs = attrs
            if self.parent is not None:
                self.parent.df = self.matrix
            return self.matrix
        if (self.top is None) and (self.parent is not None):
            self.top = self.parent.top 
        for index in self.matrix.index:
            for column in self.matrix.columns:
                i = len(self.top.residue(index)._atoms)
                j = len(self.top.residue(column)._atoms)
                self.matrix.loc[index,column] = self.matrix.loc[index,column]/(i+j)
        if self.parent is not None:
            self.parent.df = self.matrix
        return self.matrix



    
