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
        self._output = 'dist.csv'
        self.job_name = 'dist'
        self.nprocs = 'auto'
        self.compress = False
        self.selection = 'all'
        self.df = pd.DataFrame()
        self._iterload = False
        self.method = None
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
    
    @staticmethod
    def _exclude_chains(pairs, top, method='residue'):
        p = []
        for pair in pairs:
            i = pair[0]
            k = pair[1]
            if method == 'residue':
                if top.residue(i).chain == top.residue(k).chain:
                    pass
                else:
                    p.append([i,k])
            else:
                if top.atom(i).residue.chain == top.atom(k).residue.chain:
                    pass
                else:
                    p.append([i,k])
        return np.array(p)

    @staticmethod 
    def euclidean(x,y):
        assert x.shape == y.shape
        distances = []
        for index in range(len(x)):
            j = x[index]
            k = y[index]
            distances.append(np.sqrt((j[0] - k[0])**2 + (j[1] - k[1])**2 + (j[2] - k[2])**2))
        return np.array(distances)

    def exclude_chains(self, pairs, method='residue'):
        p = []
        for pair in pairs:
            i = pair[0]
            k = pair[1]
            if (method == 'residue'):
                if self.top.residue(i).chain.index == self.top.residue(k).chain.index:
                    pass
                else:
                    p.append([i,k])
            else:
                if self.top.atom(i).residue.chain.index == self.top.atom(k).residue.chain.index:
                    pass
                else:
                    p.append([i,k])
        return np.array(p)

    def saltbridge_selection(self):
        sele = self.top.select('(name NZ and resn LYS)')
        sele = np.append(sele, self.top.select('(name OE1 and resn GLU)'))
        pairs = self.atom_pairs(sele)
        if self._exclude_chain:
            pairs = self.exclude_chains(pairs, method='atom')
        true_pairs = []
        for pair in pairs:
            x = pair[0]
            y = pair[1]
            rx = self.top.atom(x).residue
            ry = self.top.atom(y).residue
            if rx.name == ry.name:
                continue
            else:
                true_pairs.append([x, y])
        return np.array(true_pairs)

    def ligand_residue_pairs(self, lig_name):
        valid_res = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD']
        residue_indeces = [res.index for res in self.top.residues if res.name in valid_res]
        ligand_index = [res.index for res in self.top.residues if res.name == '{}'.format(lig_name)][0]
        p = []
        for ri in residue_indeces:
            p.append([ri, ligand_index])
        return np.array(p)

    def run(self, method='residue', output='dist.csv',residues=[], cutoff=0.35, **kwargs):
        self._output = output
        self.method = method
        self.save()
        if method == 'atom':
            df = self.by_atom()
            return df
        elif method == 'residue':
            # pairs = self.pairs()
            if residues == []:
                pairs = self.residue_pairs([i for i in range(self.traj.n_residues)])
            else:
                # for res in self.top.residues:
                    # if (res.resSeq in residues):
                    #     print(res.name, res.resSeq, res.index)
                pairs = self.residue_pairs([res.index for res in self.top.residues if res.resSeq in residues])
            if self._exclude_chain:
                pairs = self.exclude_chains(pairs, method='residue')
            df = self.by_residue(contact_pairs=pairs, **kwargs)
            return df
        elif method == 'saltbridge':
            self.saltbridge()
            self.occupancy(cutoff=cutoff, w=True)
            return self.df
        else:
            raise ValueError('Method must be: residue, atom, saltbridge (strings)')

    def saltbridge(self):
        pairs = self.saltbridge_selection()
        distances = mdtraj.compute_distances(self.traj, atom_pairs=pairs)
        labels = []
        for pair in pairs:
            x, y = pair
            rx = self.top.atom(x).residue
            ry = self.top.atom(y).residue
            cx = self.chain_conversions()[rx.chain.index]
            cy = self.chain_conversions()[ry.chain.index]
            parts = [rx.name, rx.resSeq, cx, ry.name, ry.resSeq, cy]
            label = '{}{}{}_{}{}{}'.format(*parts)
            labels.append(label)
        self.df = pd.DataFrame(distances, columns=labels)
        if self.output is not None:
            print('Writing {}'.format(self.output))
            print('Shape: {}'.format(self.df.shape))
            self.df.to_csv(self.output)
        return self.df

    def by_atom(self):
        if isinstance(self.selection, str):
            sele = self.top.select(self.selection)
        else:
            sele = self.selection
        print(sele)
        pairs = self.atom_pairs(sele, excluded_neighbors=2)
        if self._exclude_chain:
            pairs = self.exclude_chains(pairs)
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

    def ligand_residue(self, contact_pairs, squareform=False, **kwargs):
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

    def occupancy(self, df = None, cutoff = 0.35, w=False):
        if df is None:
            if self.df is not None:
                df = self.df
            else:
                raise ValueError('No dataframe loaded or passed')
        occ = pd.DataFrame()
        n_total = 0
        for column in df.columns:
            x = df[column].to_numpy()
            n = (x <= cutoff).sum()
            occ[column] = [(n / len(x))*100]
            n_total = n_total + n
        n_percent = (n_total / len(df.index))*100
        occf = os.path.join(self.root, 'occupancy.dat')
        f = open(occf, 'w')
        f.write('# Salt Bridge Occupancy\n')
        f.write('# Cutoff = {}\n'.format(cutoff))
        f.write('# Percent (total) = {}\n'.format(n_percent))
        f.close()
        print('Total occupancy: {}'.format(n_percent))
        print('By column:')
        print(occ)
        if w:
            occ.to_csv(occf, mode='a')
            print('Wrote {}'.format(occf))

        return occ, n_percent

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


    def peptideCOM(self, stride, output='dist.csv'):
        self._output = output
        self.save()
        peptides = []
        # p = Peptide(self.topfile, xtc=self.inp, selection='all', name='idk')
        # print(p.topology.to_dataframe())
        for chain in self.top.chains:
            selection = 'chainid {}'.format(chain.index)
            p = Peptide(self.topfile, xtc=self.inp, stride=stride, selection=selection, name='chain {}'.format(chain.index))
            p.com()
            peptides.append(p)
            # print(p.topology.to_dataframe())
        data = []
        labels = []
        for i in range(len(peptides)):
            for j in range(i, len(peptides)):
                if i == j:
                    continue
                label = '{}_{}'.format(i,j)
                labels.append(label)
                p1 = peptides[i]
                p2 = peptides[j]
                # print(p1.data[0])
                # print(p2.data[0])
                d = self.euclidean(p1.data, p2.data)
                data.append(d)
        df = pd.DataFrame(data)
        df = df.T 
        df.columns = labels
        #output
        df.to_csv(self.output)
        self.df = df
        return df




            



    
