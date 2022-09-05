import os
import argparse
from matplotlib import colors
from matplotlib import pyplot as plt
import pandas as pd

import sys
# imports from pymd
sys.path.append(os.getcwd())
from pymd.mdanalysis.xpm import XpmParser
import inspect

class MixedParser(XpmParser):
    def __init__(self, _input, num_residues_1, num_residues_2, num_peptides_1, num_peptides_2, igncaps=True, ):
        if _input is None:
            self.input = None
        elif len(_input) == 1:
            self.input = _input[0]
        elif isinstance(_input, str):
            self.input = os.path.abspath(_input)
        elif isinstance(_input, list):
            self.input = []
            for item in _input:
                self.input.append(os.path.abspath(item))
        else:
            self.input = None
        self.igncaps = igncaps
        self.num_residues_1 = num_residues_1
        self.num_residues_2 = num_residues_2
        self.num_peptides_1 = num_peptides_1
        self.num_peptides_2 = num_peptides_2
        self.order = []
    
    def getMatrixMixed(self, values=None, _input=None):
        temp = self.igncaps
        self.igncaps = False
        if values is None:
            values = self.getValueAssignment(_input=_input)
        matrix = self.getMatrix(values, _input=_input)
        self.igncaps = temp
        x = 1
        drop = []
        ndx_list = [i for i in range(1, len(matrix) + 1)]
        ndx_list.sort(reverse=True)
        matrix.index = ndx_list
        if self.igncaps == True:
            num_residues = self.num_residues_1 + 2
            for i in range(0,self.num_peptides_1):
                drop1 = x
                drop2 = (i + 1) * num_residues
                x = drop2 + 1
                drop.append(drop1)
                drop.append(drop2)
            num_residues = self.num_residues_2 + 2
            for i in range(0, self.num_peptides_2):
                drop1 = x
                drop2 = x + num_residues - 1
                x = drop2 + 1
                drop.append(drop1)
                drop.append(drop2)
        if self.igncaps == True:
            len_df = (self.num_peptides_1 * self.num_residues_1) + (self.num_peptides_2 * self.num_residues_2) + (2 * (self.num_peptides_1 + self.num_peptides_2))
            matrix.columns = [i for i in range(1, len_df+1)]
        matrix.drop(drop, inplace=True, axis=0)
        matrix.drop(drop, inplace=True, axis=1)
        ndx_list = [i for i in range(1, len(matrix.index) + 1)]
        ndx_list.sort(reverse=True)
        matrix.index = ndx_list
        matrix.columns = [i for i in range(1, len(matrix.columns) + 1)]
        return matrix

    def getPeptidesMixed(self, matrix=None):
        if matrix is None:
            values = self.getValueAssignment()
            matrix = self.getMatrixMixed(values)
        _all = []

        # reminder to self: df.loc[[indeces],[coumns]]

        # peptide group 1 and itself
        indeces = self.num_peptides_1 * self.num_residues_1
        lis = [i for i in range(1, indeces + 1)]
        df = matrix.loc[lis, lis]
        new = self.getPeptides(matrix=df, num=self.num_residues_1)
        for df in new:
            _all.append(df)

        # peptide group 2 and itself
        stop  = list(matrix.index)[0]
        start = (self.num_peptides_1 * self.num_residues_1) + 1
        lis = list([i for i in range(start, stop+1)])
        df = matrix.loc[lis, lis]
        new = self.getPeptides(matrix=df, num=self.num_residues_2)
        for df in new:
            _all.append(df)

        # both peptide groups together
        start = (self.num_peptides_1 * self.num_residues_1) + 1
        stop = list(matrix.index)[0] + 1
        ndx_lis = [i for i in range(start, stop)]
        start = 1
        stop = (self.num_peptides_1 * self.num_residues_1) + 1
        column_lis = [i for i in range(start, stop)]
        split_matrix = matrix.loc[ndx_lis, column_lis]

        # split vertically
        smaller = []
        shape = split_matrix.shape[1]
        for i in range(0, shape, self.num_residues_1):
            df = split_matrix.iloc[:, i:i+self.num_residues_1]
            smaller.append(df)

        # split horizontally
        dfs = []
        shape = split_matrix.shape[0]
        for item in smaller:
            for i in range(0, shape, self.num_residues_2):
                df = item.iloc[i:i+self.num_residues_2, :]
                dfs.append(df)
        for df in dfs:
            _all.append(df)
        return _all

    def averagePeptidesMixed(self, values=None):
        all_split = []
        if values is None:
            filename = self.input[0]
            values = self.getValueAssignment(_input=filename)
        order = []
        for i in range(1, self.num_residues_1*self.num_peptides_1, self.num_residues_1):
            order.append(i)
        start = (self.num_residues_1*self.num_peptides_1) + 1
        stop = start + (self.num_peptides_2*self.num_residues_2) - 1
        for i in range(start, stop, self.num_residues_2):
            order.append(i)
        en_order = []
        i = 1
        for item in order:
            en_order.append((i, item))
            i += 1
        for filename in self.input:
            matrix = self.getMatrixMixed(values, _input=filename)
            split = self.getPeptidesMixed(matrix)
            for item in split:
                ndx = item.index[0]
                col = item.columns[0]
                group = [0,0]
                for _ord in en_order:
                    if ndx == _ord[1]:
                        group[1] = _ord[0]
                    if col == _ord[1]:
                        group[0] = _ord[0]
                self.order.append(group)
            all_split.append(split)
        averages = self.averagePeptides(all_split=all_split, set_order=False)
        averages = self.reindexAverages(averages)
        return averages
    
    def averageResiduesMixed(self):
        matrix = self.averageAllMixed()
        averages = []
        # get peptide group 1 split off
        indeces = self.num_peptides_1 * self.num_residues_1
        lis = [i for i in range(1, indeces + 1)]
        df = matrix.loc[lis, lis]
        peptide1 = self.getPeptides(matrix=df, num=self.num_residues_1)
        reset_index = [df.reset_index(drop=True) for df in peptide1]
        reset_columns = [df.T.reset_index(drop=True).T for df in reset_index]
        combined = pd.concat(reset_columns)
        average = combined.groupby(level=0).mean()
        averages.append(average)

        # get peptide group2 split off
        stop  = list(matrix.index)[0]
        start = (self.num_peptides_1 * self.num_residues_1) + 1
        lis = list([i for i in range(start, stop+1)])
        df = matrix.loc[lis, lis]
        peptide2 = self.getPeptides(matrix=df, num=self.num_residues_2)
        reset_index = [df.reset_index(drop=True) for df in peptide2]
        reset_columns = [df.T.reset_index(drop=True).T for df in reset_index]
        combined = pd.concat(reset_columns)
        average = combined.groupby(level=0).mean()
        averages.append(average)

        # get mixed split off
        start = (self.num_peptides_1 * self.num_residues_1) + 1
        stop = list(matrix.index)[0] + 1
        ndx_lis = [i for i in range(start, stop)]
        start = 1
        stop = (self.num_peptides_1 * self.num_residues_1) + 1
        column_lis = [i for i in range(start, stop)]
        split_matrix = matrix.loc[ndx_lis, column_lis]

        # split vertically
        smaller = []
        shape = split_matrix.shape[1]
        for i in range(0, shape, self.num_residues_1):
            df = split_matrix.iloc[:, i:i+self.num_residues_1]
            smaller.append(df)

        # split horizontally
        dfs = []
        shape = split_matrix.shape[0]
        for item in smaller:
            for i in range(0, shape, self.num_residues_2):
                df = item.iloc[i:i+self.num_residues_2, :]
                dfs.append(df)
        reset_index = [df.reset_index(drop=True) for df in dfs]
        reset_columns = [df.T.reset_index(drop=True).T for df in reset_index]
        combined = pd.concat(reset_columns)
        average = combined.groupby(level=0).mean()
        averages.append(average)
        averages = self.reindexAverages(averages)
        self.order = [['group1', 'group1_residue'],['group2','group2_residue'],['group1', 'group2_residue']]
        return average

    
    def reindexAverages(self, averages):
        if isinstance(averages, list):
            fixed = []
            for average in averages:
                df = average.reindex(index=average.index[::-1])
                fixed.append(df)
        else:
            fixed = averages.reindex(index=averages.index[::-1])
        return fixed
    
    def averageAllMixed(self, values=None):
        matrices = []
        filename = self.input[0]
        if values is None:
            values = self.getValueAssignment(_input=filename)
        for filename in self.input:
            matrix = self.getMatrixMixed(values=values, _input=filename)
            matrices.append(matrix)
        combined = pd.concat(matrices)
        average = combined.groupby(level=0).mean()
        fixed = self.reindexAverages(averages=average)
        return fixed

        





            

        

# xpm = MixedParser(['Mix1_MDmat_mean.xpm', 'Mix1_MDmat_mean.xpm'], num_residues_1=7, num_residues_2=10, num_peptides_1=3, num_peptides_2=3, igncaps=True)
# averages = xpm.averageResiduesMixed()
