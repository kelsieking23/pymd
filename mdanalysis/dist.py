import multiprocessing as mp
import os
import psutil
import sys

import mdtraj
import numpy as np
import pandas as pd


class Distance:

    def __init__(self, parent=None, **kwargs):
        '''
        For now, takes a dict of parameters from run.py call
        '''
        self.parent = parent
        self._inp = None
        self._topfile = None
        self.traj = None
        self._traj = None
        self.stride = None
        self.b = None
        self.e = None
        self.res = None
        self._output = None
        self.job_name = None
        self.nprocs = 'auto'
        self.compress = False
        self.__dict__.update(kwargs)

    @property
    def inp(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self._inp)
        return self._inp

    @property
    def topfile(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self._topfile)
        return self._topfile
    
    @property
    def output(self):
        if self.parent is not None:
            if not os.path.isdir(os.path.join(self.parent.root, self.job_name)):
                os.mkdir(os.path.join(self.parent.root, self.job_name))
            if self.compress:
                self._output = os.path.splitext(self._output)[0] + '.npy'
            return os.path.join(self.parent.root, self.job_name, self._output)
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
                traj = mdtraj.load(self.inp, top=self.topfile,)
        else:
            if self.stride != 0:
                traj = mdtraj.load(self.inp, top=self.topfile, stride=self.stride)
            else:
                traj = mdtraj.load(self.inp, top=self.topfile,)
        traj = traj.superpose(traj)
        self._traj = traj
        sele = traj.top.select(self.selection)
        traj = traj.atom_slice(sele)
        if (self.e == -1):
            self.traj = traj.center_coordinates()[self.b:]
        else:
            self.traj = traj.center_coordinates()[self.b:self.e]
        self.frames = self.traj._xyz
        self.top = self.traj.top
        return self

    def getPartitions(self):
        if self.nprocs == 'auto':
            nprocs = int(mp.cpu_count() // 2)
        else:
            nprocs = self.nprocs
        if self.traj is None:
            self.load()
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
                data = {
                    'b':i+interval,
                    'e':nframes,
                    'procid':procid
                }
                partitions.append(data)
                break
        return partitions, nprocs

    
    @staticmethod
    def pairs(sele_indices, excluded_neighbors = 2):
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


    def run(self, method='atom', output='dist.csv'):
        partitions, nprocs = self.getPartitions()
        pool = mp.Pool(processes=nprocs)
        self.load()
        # try:
        #     pool = mp.Pool(processes=nprocs)
        #     results = pool.map(self.by_atom, partitions)
        # except Exception as e:
        #     print(e)
        #     sys.exit(1)
        results = pool.map(self.by_atom, partitions)
        df = pd.concat(results)
        print(df)
        # df.to_csv(output)
        return df
    
    def by_atom(self):
        sele = self.top.select(self.selection)
        pairs = self.pairs(sele, excluded_neighbors=2)
        distances = mdtraj.compute_distances(self.traj, atom_pairs=pairs)
        labels = ['_'.join(list(map(str, x))) for x in pairs]
        self.df = pd.DataFrame(distances, columns=labels)
        if self.output is not None:
            print('Writing {}'.format(self.output))
            print('Shape: {}'.format(self.df.shape))
            np.save(self.output, distances)
        return self.df





    
