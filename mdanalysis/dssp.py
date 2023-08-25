import os
import sys
import mdtraj
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from pymd.mdanalysis.analysis import Analysis
from pymd.mdanalysis.postprocess import PostProcess
from pymd.plot.utils import create_colormap
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
            self._output = 'dssp.csv'
        self.job_name = 'dssp'
        self.job_params = {}
        self.parent = parent
        self.df = pd.DataFrame()
        self.block_average = block_average
        self.verbose = True
        self.__dict__.update(kwargs)
        self.test = 0
    
    def dsspOverTime(self) -> pd.DataFrame:
        if self.parent is not None:
            print('Calculating DSSP over time for {}'.format(self.parent.id))
        else:
            print('Calculating DSSP over time')
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
            if self.output.endswith('xvg'):
                self.writeXVG(method='time')
                df = PostProcess.metadata(self.output, df=self.df)
            else:
                df.to_csv(self.output)
        self.df = df
        # self.parent.df = self.df
        self.save()

        return self.df

    def dsspPerResidue(self, replicate_average=True, interval=1000, system_average=False):
        if self.traj is None:
            raise ValueError('Trajectory not loaded')
        self.save()
        assignments = mdtraj.compute_dssp(self.traj, simplified=True)
        df = pd.DataFrame(assignments)
        encoded = self.ssEncoder(df)
        self.df = encoded.T
        if replicate_average:
            dfs = []
            for i in range(0, self.traj.n_frames, interval):
                _df = self.df.loc[:,i:i+interval]
                _df.columns = [i for i in range(len(_df.columns))]
                _df.index = [i for i in range(len(_df.index))]
                dfs.append(_df)
            self.df = pd.concat(dfs).groupby(level=0).mean()
        labels = []
        for residue in self.top.residues:
            label = '{}{}_{}'.format(residue.name, residue.resSeq, residue.chain.index)
            labels.append(label)
        self.df.index = labels
        self.df.to_csv(self.output)
        return self.df
    
    def plotPerResidue(self, df=None, chain_average=False, output='dssp.perresidue.png', colorbar=True, nrows=1, ncols=1):
        fig, axes = plt.subplots(nrows, ncols, constrained_layout=True)
        normal_w = 8
        normal_h = 6
        fig_h = (normal_h * nrows) 
        fig_w = (normal_w * ncols) + 2
        fig.set_size_inches(fig_w*.7, fig_h*.7)
        if df is None:
            df = self.df.T
        else:
            df = df.T
        cmap = create_colormap(((255,255,255),(33,52,104),(246,140,62)), bit=True)
        k = 0
        ims = []
        for chain, ax in zip(self.top.chains, axes.flat):
            chain = [col for col in df.columns if col.endswith(str(chain.index))]
            x = df[chain].T

            im = ax.imshow(x, aspect='auto', cmap = cmap, vmin=0, vmax=1)
            ax.yaxis.set_major_locator(MultipleLocator(4))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
            ax.set_yticklabels([i for i in range(-3,42,4)])
            # ax.xaxis.set_major_locator(MultipleLocator(len(df.columns)*0.2))
            # ax.xaxis.set_minor_locator(MultipleLocator(len(df.columns)*0.02))
            xlabels = []
            for n in range(-200,1001,200):
                xlabels.append(n)
            for m in range(0,2):
                for n in range(0,1001,200):
                    xlabels.append(n)
            ax.set_xticklabels(xlabels)
            ax.set_title('Peptide {}'.format(k+1))
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('Residue')
            k += 1
        if colorbar:
            fig.colorbar(im, ax=axes.ravel().tolist())
        # fig.subplots_adjust(right=0.8)
        # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        # fig.colorbar(im, cax=cbar_ax)
        # plt.tight_layout()
        
        o = os.path.join(self.root, output)
        plt.savefig(o, dpi=300)
        if not sys.platform == 'linux':
            plt.show()
        return fig, axes

    def ssEncoder(self, df):
        df = df.replace(['C','E','H'], [0,2,1])
        return df
    
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
        print('Wrote {}'.format(self.output))

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
