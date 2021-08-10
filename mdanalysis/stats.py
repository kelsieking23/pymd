import os
from itertools import groupby
import numpy as np
import pandas as pd
from scipy import stats

import sys
if os.name == 'nt':
    sys.path.append('D:/Work/iapp/')
    from pymd.structure.protein import Protein
    from pymd.mdanalysis.postprocess import PostProcess
if os.name == 'posix':
    sys.path.append('/work/cascades/kelsieking23/iapp_analysis/scripts/python/')
    from structure.protein import Protein
    from mdanalysis.postprocess import PostProcess

class StatAnalysis:
    '''
    Performs statistical analysis.
    '''

    def __init__(self, systems=None, control=None):
        '''
        Initialize StatAnalysis.
        Arguments:
        *** systems (iterable [pymd.mdanalysis.system.System], optional): iterable containing
            System objects (pymd.mdanalysis.system.System)
        *** control (pymd.mdanalysis.system.System, optional): the System object to use as control.
            if not specified, will infer. 
        '''
        if systems is not None:
            self.systems = systems
            self.control = control
        else:
            self.systems = None
            self.control = None

    

    def mindist(self, residue_ids, group='sidechain', saveto=None, csv_files=None, control=None):
        '''
        Performs t-test on average number of contacts between a given residue across systems.
        Arguments:
        *** residue_ids (iterable [str], or str): residue ID(s) to compare distances (in format RES###, ex: PHE23).
            Must be of length 2 (for now)
        *** group (str, default='sidechain'): the mindist group for analysis. 
        *** csv_files (iterable [str], optional): list of csv files to pull data from. csv should be
            in format of output from pymd.postprocess.mindist(). if not specified, will use
            files from self.system
        *** control (str, optional): name of CSV file containing control. must be specified if csv_files is
            not None. if not specified, will infer. 
        Returns:
        ** (DataFrame) returns dataframe consisting of average number of contacts per system, and p-value
            compared to control. 
        '''
        # TODO: add use case for csv_files: will need topology if no self.system specified. 
        ctrl = False
        data = {}
        for system in self.systems:
            if system is self.control:
                ctrl = True
            print(system.name)
            # get residue indexing 
            if isinstance(residue_ids, str):
                indeces = self.getResidueIndeces(residue_ids, system)
            else:
                if not self.allEqual(residue_ids):
                    if len(residue_ids) > 2:
                        raise ValueError('Must be a list of length 2')
                        # TODO: fix this
                    indeces = {}
                    for residue_id in residue_ids:
                        indeces[residue_id] = self.getResidueIndeces(residue_id, system)
                else:
                    indeces = self.getResidueIndeces(residue_ids[0], system)

            # get system average from distance.csv files
            sys_dfs = []
            for rep in range(1, system.reps+1):
                csv = os.path.join(system.directory[rep]['mindist'][group]['root'], 'distances.csv')
                if not os.path.isfile(csv):
                    post = PostProcess(system=system)
                    post.mindist(dtype=group, rep=rep)
                df = pd.read_csv(csv, index_col=0, header=0)
                sys_dfs.append(df)
            df = pd.concat(sys_dfs).groupby(level=0).mean()
            df = df.reset_index(drop=True)
            df = df.T.reset_index(drop=True).T
            
            # get distances
            distances = []
            if isinstance(indeces, list):
                for i in range(0, len(indeces)):
                    r1 = indeces[i]
                    for j in range(0, len(indeces)):
                        r2 = indeces[j]
                        print(r2)
                        if r1 == r2:
                            continue
                        else:
                            print(r1,r2)
                            distances.append(df.iloc[r1,r2])
            if isinstance(indeces, dict):
                indeces1 = []
                indeces1 = []
                i = 0
                for key in indeces.keys():
                    if i == 0:
                        indeces1 = indeces[key]
                    if i == 1:
                        indeces2 = indeces[key]
                    i += 1
                for i in range(0, len(indeces1)):
                    r1 = indeces1[i]
                    for j in range(0, len(indeces2)):
                        r2 = indeces2[j]
                        distances.append(df.iloc[r1,r2])
            if ctrl is True:
                data['control'] = distances
                ctrl = False
            else:
                data[system.name] = distances

        d = {}  
        # do t tests
        for key in data.keys():
            avg = sum(data[key])/len(data[key])
            std = np.std(np.array(data[key]))
            d[key] = {}
            d[key]['average'] = avg
            d[key]['std'] = std
            if key == 'control':
                d[key]['tstat'] = 0
                d[key]['pval'] = 0
            else:
                tstat, pval = stats.ttest_ind(data['control'], data[key])
                d[key]['tstat'] = tstat
                d[key]['pval'] = pval
        df = pd.DataFrame(d)

        if saveto is not None:
            df.to_csv(saveto)

        return df

    def dsspOverTime(self, xvg, time=None, block_average=None, saveto=None):
        dfs = {}
        for system in self.systems:
            _xvg = xvg
            if system.name == 'TAX':
                _xvg = 'dssp.xvg'
            df = system.post.dsspOverTime(xvg=_xvg, block_average=block_average, average_reps=True)
            if system.name == self.control.name:
                if time is None:
                    dfs['control'] = df
                else:
                    dfs['control'] = df.loc[time[0]:time[1],:]
            else:
                if time is None:
                    dfs[system.name] = df
                else:
                    dfs[system.name] = df.loc[time[0]:time[1],:]
        
        d = {
            'coil':{},
            'bsheet':{},
            'helix':{}
        }
        for key in dfs.keys():
            if key == 'control':
                continue
            print(key)
            tsat, pval = stats.ttest_ind(dfs['control']['coil_percent'], dfs[key]['coil_percent'])
            d['coil'][key] = pval
            print(dfs[key]['coil_percent'])
            print(dfs['control']['coil_percent'])
            tstat, pval = stats.ttest_ind(dfs['control']['bsheet_percent'], dfs[key]['bsheet_percent'])
            d['bsheet'][key] = pval
            tstat, pval = stats.ttest_ind(dfs['control']['helix_percent'], dfs[key]['helix_percent'])
            d['helix'][key] = pval
        df = pd.DataFrame(d)

        if saveto is not None:
            df.to_csv(saveto)

        return df
            







    
    @staticmethod
    def getResidueIndeces(residue_id, system):
        res_index = system.protein.ids.index(residue_id) 
        if system.peptides > 1:
            res_indeces = [res_index]
            for i in range(1, system.peptides):
                res_index += len(system.protein.ids)
                res_indeces.append(res_index)
        else:
            res_indeces = [res_index]
        return res_indeces


    @staticmethod
    def allEqual(iterable):
        g = groupby(iterable)
        return next(g, True) and not next(g, False)

