import os
import sys
import mdtraj
import numpy as np
import pandas as pd
from pymd.mdanalysis.analysis import Analysis
from pymd.mdanalysis.postprocess import PostProcess
from typing import Optional, Any
import warnings

class DSSP(Analysis):

    def __init__(self, inp: str, top: str, parent: Optional[Any] = None, block_average: Optional[int] = None, **kwargs: Any):
        self.parent = parent
        self._inp = inp
        self._topfile = top
        self.traj = None
        self._traj = None
        if 'output' in kwargs.keys():
            self._output = kwargs['output']
        else:
            self._output = None
        self.job_name = 'analysis'
        self.job_params = {}
        self.parent = parent
        self.df = pd.DataFrame()
        self.block_average = block_average
        self.__dict__.update(kwargs)
        self.test = 0
    
    def dsspOverTime(self) -> pd.DataFrame:
        if self.traj is None:
            raise ValueError('Trajectory not loaded')
        assignments = mdtraj.compute_dssp(self.traj, simplified=True)
        arr = np.array(list(map(self.calcPercentages, assignments)))
        df = pd.DataFrame(arr, columns=['helix', 'bsheet', 'coil'])
        df.index = self.traj.time
        if self.block_average is not None:
            df = self.blockAverage(df)
        # df.index = [i*0.001 for i in df.index] # type: ignore
        # if self.block_average is not None:
        #     df.index = [i*10 for i in df.index] # type: ignore
        self.df = df
        
        if self.output is not None:
            self.writeXVG(method='time')
            df = PostProcess.metadata(self.output, df=self.df)
        self.df = df
        self.parent.df = self.df
        self.save()

        return self.df
    def calcPercentages(self, assignments: pd.DataFrame) -> np.ndarray:
        a = assignments[assignments != 'NA']
        # alpha
        h = len(a[a == 'H'])
        g = len(a[a == 'G'])
        i = len(a[a == 'I'])
        alpha = ((h + g + i) / len(a)) * 100
        # beta
        b = len(a[a == 'E'])
        e = len(a[a == 'B'])
        beta = ((b + e) / len(a)) * 100
        # coil
        c = len(a[a == 'C'])
        t = len(a[a == 'T'])
        s = len(a[a == 'S'])
        irr = len(a[a == ' '])
        coil = ((c + t + s + irr) / len(a)) * 100
        return np.array([alpha, beta, coil])

    def blockAverage(self, df):
        ss = {
            'coil':{},
            'bsheet':{},
            'helix':{}
        }
        i = 0
        cp = 0
        bp = 0
        hp = 0
        if (self.block_average is None) or (self.block_average == 0):
            warnings.warn('Block average is None - returning original DataFrame')
            return df
        for index in df.index:
            t = index
            if i != self.block_average:
                cp = cp + df.loc[index, 'coil']
                bp = bp + df.loc[index, 'bsheet']
                hp = hp + df.loc[index, 'helix']
            else:
                cp = cp / (self.block_average+1)
                bp = bp / (self.block_average+1)
                hp = hp / (self.block_average+1)
                ss['coil'][t] = cp
                ss['bsheet'][t] = bp
                ss['helix'][t] = hp
                cp = cp + df.loc[index, 'coil']
                bp = bp + df.loc[index, 'bsheet']
                hp = hp + df.loc[index, 'helix']
                i = 0
            i += 1    
        df = pd.DataFrame(ss)
        return df  

    def writeXVG(self, method='time'):
        f = open(self.output, 'w') # type:ignore
        f.write('# This file was created {}\n'.format(self.now()))
        f.write(f'# Created by: pymd.mdanalysis.dssp\n')
        f.write(f'@    title "DSSP"\n')
        f.write(f'@    xaxis label "Time (ps)"\n')
        f.write(f'@    yaxis label "Secondary Structure (%)"\n')
        f.write(f'@TYPE xy\n')
        if method == 'time':
            f.write(f'@PTYPE dssp\n')
        i = 0
        for column in self.df.columns:
            f.write('@ s{} legend "{}"\n'.format(i, column))
            i += 1
        i = 0
        for index in self.df.index:
            data = [self.traj.time[i]] + list(self.df.loc[index,:])
            fmt = '{:>10.2f}   {:>6.3f}   {:>6.3f}   {:>6.3f}\n'.format(*data)
            f.write(fmt)
            i += 1
        f.close()

    # def loadTrajectory(self, stride=100, selection='backbone', b=0, e=-1):
    #     if stride != 0:
    #         traj = mdtraj.load(self.inp, top=self.topfile, stride=stride)
    #     else:
    #         traj = mdtraj.load(self.inp, top=self.topfile)
    #     traj = traj.superpose(traj)

    #     self._traj = traj
    #     sele = traj.top.select(selection)
    #     traj = traj.atom_slice(sele)
    #     if (e == -1):
    #         self.traj = traj.center_coordinates()[b:]
    #     else:
    #         self.traj = traj.center_coordinates()[b:e]
    #     return self

    @staticmethod
    def averageDSSP(dfs):
        df = pd.concat(dfs).groupby(level=0).mean()
        new_cols = ['helix_std', 'bsheet_std', 'coil_std']
        cols = ['Helix', r'$\beta$-Strand', 'Coil']
        k = 0
        for col in cols:
            stdevs = pd.DataFrame()
            i = 1
            for df in dfs:
                print(i)
                stdevs[i] = df[col]
                i += 1
            df[new_cols[k]] = stdevs.std(axis=1)
            k += 1
        return df
# d = DSSP('/Users/kelsieking/Desktop/bblab/QUR_CHARMM/1/cat.pbc.nowat.xtc', '/Users/kelsieking/Desktop/bblab/QUR_CHARMM/1/nowat.top.gro')
# x = d.dsspOverTime()
# x = pd.DataFrame(x, columns=['coil', 'beta', 'helix'])
# print(x)
