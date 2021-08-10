import multiprocessing
import networkx as nx
import numpy as np
import pandas as pd
import os
from scipy.spatial import distance_matrix
from multiprocessing import Queue, Process

from numba import cuda
from numba import jit
from numba import vectorize

import sys
if os.name == 'nt':
    sys.path.append('D:/Work/iapp/')
if os.name == 'posix':
    sys.path.append('/mnt/d/Work/iapp')
from pymd.structure.protein import Protein

class Multiprocessor:
    '''
    does the multiprocessing 
    '''
    
    def __init__(self):
        self.processes = []
        self.queue = Queue()
        self.rets = None
        
    @staticmethod
    def _wrapper(func, queue, args, kwargs):
        ret = func(*args, **kwargs)
        queue.put(ret)
    
    def run(self, func, *args, **kwargs):
        p = Process(target=self._wrapper, args=(func, self.queue, args, kwargs,))
        self.processes.append(p)
        p.start()
    
    def wait(self):
        for p in self.processes:
            if self.rets is None:
                self.rets = self.queue.get()
            else:
                self.rets = np.add(self.rets, self.queue.get())
        for p in self.processes:
            p.join()
        processed = len(self.processes)
        self.processes = []
        return processed


def getPositionMatrices():
    '''
    Get xyz for all frames into numpy arrays. returns list. 
    '''
    ### also need to get the time at some point for implementation (skipping frames and such)
    frame = []
    f = open(os.path.join('ctrl_cat_pbc_550_600_dumped.xtc'), 'r')
    for line in f:
        if line.strip().startswith('natoms'):
            if frame != []:
                arr = np.array(frame, dtype='float32')
                yield arr
                frame = []
        if line.strip().startswith('x['):
            coordinates_split = line.strip().split('=')[-1].split(',')
            coordinates = []
            for split_coord in coordinates_split:
                string = ''
                for char in split_coord:
                    if (char == '{') or (char == ',') or (char == '}') or (char == ' '):
                        continue
                    else:
                        string = string + char
                coordinate = float(string)*10
                coordinates.append(coordinate)
            frame.append(coordinates)
    f.close()

def getDistanceMatrices(x):
    return distance_matrix(x,x).astype('float32')

@vectorize(['int64(float32)'], target='cuda')
def isInteraction(x):
    if x <= 3:
        return 1
    else:
        return 0

def getBinaryMatrices(x):
    d = getDistanceMatrices(x)
    return isInteraction(d)

def getResidueTerimanls(protein):
    indeces = []
    for residue in protein.residues:
        indeces.append(residue.atoms[-1].index)
    return indeces


def transformToResidueMatrix(indeces, x):
    arr = []
    last_i = -1
    last_j = -1
    for n in range(len(indeces)):
        i = indeces[n]
        row = []
        for m in range(len(indeces)):
            j = indeces[m]
            if i == j:
                row.append(0)
            elif not np.any(x[last_i+1:i+1,last_j+1:j+1]):
                row.append(0)
            else:
                row.append(1)
            last_j = j
        last_i = i
        arr.append(row)
    res = np.array(arr)
    return res

def dump():
    script_dir = os.getcwd()
    path = '/mnt/d/Work/iapp/systems/L500_PDB/Trajectories/IAPP'
    os.chdir(path)
    for folder in os.listdir(os.getcwd()):
        cwd = os.path.join(os.getcwd(), folder)
        os.chdir(cwd)
        for file_ in os.listdir(cwd):
            if file_.startswith('centered'):
                dumped = 'dump_' + file_
                os.system('/mnt/d/Work/iapp/pymd/mdanalysis/dump.sh "{}" "{}"'.format(file_, dumped))
        os.chdir(path)

if __name__=='__main__':
    import time
    start = time.time()
    tot = None
    # dump()
    i = 0
    j = 0
    path = 'D:/Work/iapp/systems/L500_PDB/Trajectories/IAPP/1'
    for file_ in os.listdir(path):
        if file_.startswith('dump'):
            print(file_)
            protein = Protein(structure='D:/Work/iapp/systems/L500_PDB/Trajectories/IAPP/1/md_1900_2000.part0022.gro')
            indeces = getResidueTerimanls(protein)
            for frame in getPositionMatrices():
                m = getBinaryMatrices(frame)
                m = transformToResidueMatrix(indeces, m)
                if tot is None:
                    tot = m
                else:
                    tot = np.add(tot, m)
                i += 1
            df = pd.DataFrame(tot)
            df.to_csv('iapp{}_finished_res.csv'.format(j))
            x = tot / i
            df = pd.DataFrame(x)
            df.to_csv('iapp{}_probs_res.csv'.format(j))
            print(time.time() - start)
            j += 1


