import os
import argparse
from matplotlib import colors
from matplotlib import pyplot as plt
import pandas as pd


class XpmParser:
    def __init__(self, _input, num_residues, num_peptides=None, igncaps=True):
        if len(_input) == 1:
            self.input = _input[0]
        elif isinstance(_input, str):
            self.input = os.path.abspath(_input)
        else:
            self.input = []
            for item in _input:
                self.input.append(os.path.abspath(item))
        self.igncaps = igncaps
        self.num_residues = num_residues
        self.num_peptides = num_peptides
        self.vmin = None
        self.vmax = None
        self.order = []

    def getValueAssignment(self, _input=None):
        values = {}
        if _input is not None:
            filename = _input
        else:
            filename = self.input
        with open(filename) as f:
            contents = f.readlines()[9:]
        f.close()
        for line in contents:
            line_parts = line.strip().split()
            if '/*' not in line_parts[0]:
                values[line_parts[0][1:]] = float(line_parts[5][1:-1])
            else:
                break
        keys = list(values.keys())
        self.vmin = values[keys[0]]
        self.vmax = values[keys[-1]]
        return values
    
    def getAxes(self, _input=None):
        x = []
        y = []
        if _input is not None:
            filename = _input
        else:
            filename = self.input
        with open(filename) as f:
            contents = f.readlines()[9:]
        for line in contents:
            if line.startswith('/* x-axis:'):
                values = line.split(':')[1].strip().split()
                for value in values:
                    if (value != '*/'):
                        x.append(float(value))
            if line.startswith('/* y-axis:'):
                values = line.split(':')[1].strip().split()
                for value in values:
                    if (value != '*/'):
                        y.append(float(value))
        return x, y
            

    def getMatrix(self, values=None, _input=None):
        if values is None:
            if _input is not None:
                values = self.getValueAssignment(_input)
            else:
                values = self.getValueAssignment(self.input)
        if _input is not None:
            if isinstance(_input, str):
                f = open(_input, 'r')
                contents = f.readlines()
                f.close()
            if isinstance(_input, list):
                contents = _input
        else:
            f = open(self.input, 'r')
            contents = f.readlines()
            f.close()
        data = []
        matrix = {}
        skip_next = False
        for i in range(len(contents)):
            line = contents[i]
            line_parts = line.split()
            if 'static' in line_parts[0]:
                skip_next = True
                continue
            elif skip_next == True:
                skip_next = False
                continue
            elif ('/*' in line_parts[0]) or  ('*/' in line_parts[-1]):
                continue
            else:
                data.append(line)
        if self.igncaps == True:
            num = self.num_residues + 2
            overall_counter = 1
            line_counter = 1
            for line in data:
                if line_counter == 1:
                    line_counter += 1
                    continue
                if line_counter == num:
                    line_counter = 1
                    continue
                else:
                    _values = []
                    char_counter = 1
                    # stripped = line[1:-3]
                    if line == data[-1]:
                        stripped = line[1:-2]
                    else:
                        stripped = line[1:-3]
                    for char in stripped:
                        if char_counter == 1:
                            char_counter += 1
                            continue
                        elif char_counter == num:
                            char_counter = 1
                            continue
                        else:
                            value = values[char]
                            _values.append(value)
                            char_counter += 1
                    line_counter += 1
                    matrix[overall_counter] = _values
                overall_counter += 1
        else:
            line_counter = 1
            for line in data:
                _values = []
                if line != data[-1]:
                    stripped = line[1:-3]
                else:
                    stripped = line[1:-2]
                for char in stripped:
                    value = values[char]
                    _values.append(value)
                matrix[line_counter] = _values
                line_counter += 1
        df = pd.DataFrame.from_dict(matrix)
        df.reset_index(inplace=True, drop=True)
        df = df.T.reset_index(drop=True).T
        return df

    def averageMatrix(self, values=None, _input=None, matrix=None):
        if values is None:
            values = self.getValueAssignment()
        if matrix is None:
            matrix = self.getMatrix(values,_input)
        shape = int(matrix.shape[0])
        smaller = []
        for i in range(0, shape, self.num_residues):
            df = matrix.iloc[:,i:i+self.num_residues]
            smaller.append(df)
        dfs = []
        for item in smaller:
            for i in range(0, shape, self.num_residues):
                df = item.iloc[i:i+self.num_residues,:]
                dfs.append(df)
        reset_index = [df.reset_index(drop=True) for df in dfs]
        reset_columns = [df.T.reset_index(drop=True).T for df in reset_index]
        combined = pd.concat(reset_columns)
        averaged = combined.groupby(level=0).mean()
        return averaged
    
    def getPeptides(self, values=None, _input=None, matrix=None, num=None):
        if matrix is None:
            if values is None:
                if _input is None:
                    _input = self.input[0]
                values = self.getValueAssignment(_input)
            matrix = self.getMatrix(values,_input)
        shape = int(matrix.shape[0])
        smaller = []
        if num is None:
            num = self.num_residues
        else:
            num = num
        for i in range(0, shape, num):
            df = matrix.iloc[:,i:i+num]
            smaller.append(df)
        dfs = []
        for item in smaller:
            for i in range(0, shape, num):
                df = item.iloc[i:i+num,:]
                dfs.append(df)
        return dfs
    
    def averagePeptides(self, values=None, all_split=None, set_order=True):
        if values is None:
            _input = self.input[0]
            values = self.getValueAssignment(_input=_input)
        if set_order == True:
            order = []
            for i in range(1, self.num_residues*self.num_peptides, self.num_residues):
                order.append(i)
            en_order = []
            i = 1
            for item in order:
                en_order.append((i, item))
                i += 1
            print(order)
            print(en_order)
        if all_split is None:
            all_split = []
            for filename in self.input:
                matrix = self.getMatrix(values=values, _input=filename)
                # matrix.index = [i for i in range(index_start+1, index_end+2)]
                ndx_list = [i for i in range(1, len(matrix) + 1)]
                ndx_list.sort(reverse=True)
                matrix.index = ndx_list
                split = self.getPeptides(matrix=matrix)
                if set_order == True:
                    for item in split:
                        index_start = item.index[0]
                        index_end = item.index[-1]
                        ndx = item.index[-1]
                        col = item.columns[0]
                        group = [0,0]
                        for _ord in en_order:
                            if ndx == _ord[1]:
                                group[1] = _ord[0]
                            if col == _ord[1]:
                                group[0] = _ord[0]
                        self.order.append(group)
                all_split.append(split)
        num_groups = len(all_split[0])
        num_files = len(self.input)
        averages = []
        for i in range(num_groups):
            to_average = []
            for k in range(num_files):
                to_average.append(all_split[k][i])
            reset_index = [df.reset_index(drop=True) for df in to_average]
            reset_columns = [df.T.reset_index(drop=True).T for df in reset_index]
            combined = pd.concat(reset_columns)
            average = combined.groupby(level=0).mean()
            averages.append(average)
        return averages
    
    def averageReplicates(self, matrices=None):
        if matrices is None:
            matrices = []
            _input = self.input[0]
            values = self.getValueAssignment(_input=_input)
            for filename in self.input:
                matrix = self.getMatrix(values=values, _input=filename)
                matrices.append(matrix)
        reset_index = [df.reset_index(drop=True) for df in matrices]
        reset_columns = [df.T.reset_index(drop=True).T for df in reset_index]
        combined = pd.concat(reset_columns)
        average = combined.groupby(level=0).mean()
        return average
        
    def averageResidues(self):
        matrix = self.averageAll()
        shape = matrix.shape[0]
        num = self.num_residues
        smaller = []
        for i in range(0, shape, num):
            df = matrix.iloc[:,i:i+num]
            smaller.append(df)
        dfs = []
        for item in smaller:
            for i in range(0, shape, num):
                df = item.iloc[i:i+num,:]
                dfs.append(df)
        reset_index = [df.reset_index(drop=True) for df in dfs]
        reset_columns = [df.T.reset_index(drop=True).T for df in reset_index]
        combined = pd.concat(reset_columns)
        average = combined.groupby(level=0).mean()
        return average
    
    def averageAll(self, values=None):
        averages = []
        if values is None:
            filename = self.input[0]
            values = self.getValueAssignment(filename)
        for filename in self.input:
            averaged = self.averageMatrix(values, filename)
            averages.append(averaged)
        combined = pd.concat(averages)
        average = combined.groupby(level=0).mean()
        return average
    


# xpm = XpmParser(_input='IAPP1_MDmat_mean.xpm', num_peptides=6, num_residues=10, igncaps=True)
# print(xpm.getMatrix())
# xpm = XpmParser(['IAPP1_MDmat_mean.xpm', 'IAPP1_MDmat_mean.xpm'], num_peptides=6, num_residues=10, igncaps=True)
# dfs = xpm.averageAllResidues()
# fig, ax = plt.subplots()
# ax.imshow(dfs)
# plt.show()
