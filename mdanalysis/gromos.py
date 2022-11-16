import os
import sys
import multiprocessing as mp
import pandas as pd
import time
from typing import Optional, Any
import random

from pymd.mdanalysis.analysis import Analysis
from pymd.mdanalysis.rmsd import RMSD
from pymd.utilities.rewritepdb import writePDB

class GROMOS(Analysis):

    def __init__(self, inp: str, top: str, parent: Optional[Any] = None, cutoff: Optional[float] = 0.3, **kwargs: Any):
        self.parent = parent
        self._inp = inp
        self._topfile = top
        self.traj = None
        self._traj = None
        if 'output' in kwargs.keys():
            self._output = kwargs['output']
        else:
            self._output = None
        self.job_name = 'gromos'
        self.job_params = {}
        self.parent = parent
        self.df = pd.DataFrame()
        self.X = pd.DataFrame()
        self.size = pd.DataFrame()
        self.cutoff = cutoff
        self.central_clusters = []
        self.n_clusters = []
        self.__dict__.update(kwargs)
        self.test = 0

    def getPartitions(self, nprocs):
        if nprocs == 'auto':
            nprocs = mp.cpu_count() / 2
        nframes, _, _ = self.traj._xyz.shape
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
                print(i+interval, nframes)
                if (i + interval) != (nframes - 1):
                    data = {
                        'b':i+interval,
                        'e':nframes,
                        'procid':procid
                    }
                    partitions.append(data)
                    break
                else:
                    b = partitions[-1]['b']
                    p = partitions[-1]['procid']
                    partitions[-1] = data = {
                            'b':b,
                            'e':nframes,
                            'procid':p
                    }
                    break
        return partitions

    def run(self, nprocs, output='cluster.csv', stride=100, selection='backbone', b=0, e=-1):
        if self.traj is None:
            self.loadTrajectory(stride=stride, selection=selection, b=b, e=e)
        partitions = self.getPartitions(nprocs)
        for item in partitions:
            print(item)
        pool = mp.Pool(processes=nprocs)
        results = pool.map(self.matrix, partitions)
        df = pd.concat(results, axis=1)
        df = df.fillna(df.T)
        self.df = df
        output = os.path.join(self.root, output)
        df.to_csv(output)
        return df

    def matrix(self, data):
        start_index = data['b']
        stop_index = data['e']
        procid = data['procid']
        df = pd.DataFrame()
        i = 1
        rms = RMSD(traj = self.traj, cutoff=self.cutoff, binary=False)
        nframes, _, _ = self.traj._xyz.shape
        for reference_index in range(start_index, stop_index): 
            reference = self.frames[reference_index]
            rms.reference = reference
            rmsd = rms.rmsdOverTime(frames=self.frames[reference_index:])
            fill = [None]*(reference_index)
            if fill != []:
                df[reference_index] = fill + rmsd
            else:
                df[reference_index] = rmsd
            i += 1
        stop = time.time()
        df.index = self.traj.time
        df.columns = self.traj.time[start_index:stop_index]
        self.df = df
        return df

    def gromos(self, filepath=None, prefix=None):
        if self.df.empty:
            if prefix is not None:
                filename = os.path.join(self.root, '{}.cluster.csv'.format(prefix))
            elif filepath is not None:
                filename = filepath
            else:
                filename = os.path.join(self.root, 'cluster.csv')
            self.df = pd.read_csv(filename, index_col=0, header=0)
        df = self.df
        df.index = [i for i in range(len(df.index))]
        df.columns = [i for i in range(len(df.columns))]
        central_clusters = []
        n_clusters = []
        last_df = None
        i = -1
        while not df.empty:
            i += 1
            df, _max_index, _max = self.neighbor_search(df, self.cutoff)
            central_clusters.append(_max_index)
            n_clusters.append(_max)
            last_df = df
        self.central_clusters = central_clusters
        self.n_clusters = n_clusters
        self.clustsize()
        return central_clusters, n_clusters

    def neighbor_search(self, df, cutoff=None):
        if cutoff is None:
            cutoff = self.cutoff
        _max = None
        _max_index = None
        n_neighbors = {}
        neighbors = {}
        df.columns = df.index
        for column in df.columns:
            d = df[column]
           
            n_neighbors[column] = len((d[d <= cutoff]))
            neighbors[column] = d[d<=cutoff]
            if _max is None:
                _max = len((d[d <= cutoff]))
                _max_index = column
            else:
                if _max < len((d[d <= cutoff])):
                    _max = len((d[d <= cutoff]))
                    _max_index = column
                elif _max == len((d[d <= cutoff])):
                    if random.choice([0,1]) == 0:
                        _max = len((d[d <= cutoff]))
                        _max_index = column
                else:
                    pass
        df = df.drop(neighbors[_max_index].index, axis=0)
        df = df.drop(neighbors[_max_index].index, axis=1)
        return df, _max_index, _max

    def clustsize(self):
        df = pd.DataFrame()
        df['cluster'] = [i for i in range(len(self.n_clusters))]
        df['n_clusters'] = self.n_clusters
        df['percentage'] = [(i/len(self.frames))*100 for i in self.n_clusters]
        output = os.path.join(self.root, 'size.csv')
        df.to_csv(output)
        self.size = df
    
    def writePDB(self, wcl=5):
        for i in range(0, wcl):
            # t = int(self.central_clusters[i] / self.traj.timestep)
            t = self.central_clusters[i]
            frame = self._traj._xyz[t]
            chain_index = 0
            chain_id = 'A'
            contents = ['REMARK t={}\n'.format(self.central_clusters[i])]
            contents.append('REMARK cluster size = {}\n'.format(self.n_clusters[i]))
            contents.append('REMARK clust % = {}\n'.format((self.n_clusters[i]/self._traj.n_frames)*100))
            for z in range(0, len(frame)):
                atom = self._traj.topology._atoms[z]
                if atom.residue.chain.index > chain_index:
                    chain_index = atom.residue.chain.index
                    chain_id = chr(ord(chain_id) + 1)
                x, y, z = map(self.fixCoordinates, frame[z])
                line = ['ATOM', str(atom.index), atom.name, atom.residue.name, chain_id, str(atom.residue.resSeq), x, y, z, '1.00', '0.00', atom.element.symbol]
                contents.append(line)
            output = os.path.join(self.root, 'clust{}.pdb'.format(str(i).zfill(2)))
            writePDB(contents, output)
            print('Wrote {}'.format(output))
