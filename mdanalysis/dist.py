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
from pymd.mdanalysis.analysis import Analysis

class Distance(Analysis):

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
        self._output = None
        self.job_name = 'dist'
        self.nprocs = 'auto'
        self.compress = False
        self.selection = 'all'
        self.df = pd.DataFrame()
        self._iterload = False
        self.matrix = pd.DataFrame()
        self.__dict__.update(kwargs)

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
        self._output = output
        self.save()
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
    
    def by_residue(self, contact_pairs, squareform=False, **kwargs):
        if self.parent is not None:
            print('Computing residue contacts for {}...'.format(self.parent.id))
        distances, pairs = mdtraj.compute_contacts(self.traj, contacts=contact_pairs, **kwargs) # type: ignore
        # sq = mdtraj.geometry.squareform(distances, pairs)
        distances = pd.DataFrame(distances)
        pairs = pd.DataFrame(pairs)
        distances.columns = ['{}_{}'.format(pairs.loc[index,0], pairs.loc[index, 1]) for index in pairs.index]  # type: ignore
        self.df = distances
        # if squareform:
        #     self.df = sq
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
    
    def frequency(self, cutoff=0.6):
        df = pd.DataFrame()
        frequencies = []
        weighted = pd.DataFrame()
        for col in self.df.columns:
            freq = len(self.df[self.df[col] <= cutoff])/len(self.df.index)
            weighted[col] = self.df[col] * freq
            # i += 1
            # if i == 1:
            #     print(freq)
            frequencies.append(freq)
        df['pairs'] = self.df.columns
        df['frequency'] = frequencies
        weighted.index = self.df.index
        self.df = weighted
        return self.df, df
    
    def squareform(self, frames=(0,-1)):
        pairs = [col.split('_') for col in self.df.columns]
        pairs = np.array([list(map(int, pair)) for pair in pairs])
        distances = self.df.loc[frames[0]:frames[-1],:].to_numpy()
        squares = mdtraj.geometry.squareform(distances, pairs)
        dfs = []
        for sq in squares:
            df = pd.DataFrame(sq, columns=[residue.index for residue in self.top.residues],
                                index=[residue.index for residue in self.top.residues])
            dfs.append(df)
        return dfs




    
