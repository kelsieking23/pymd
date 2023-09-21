import math
from operator import index
import mdtraj
import numpy as np
import os
import pandas as pd
import sys



class RMSD:

    def __init__(self, traj, reference=None, cutoff=0, binary=False, time=None):
        self.traj = traj._xyz
        if reference is None:
            self.reference = self.first_frame
        else:
            self.reference = reference
        # self.time = traj.time
        self.cutoff = cutoff
        self.binary = binary
        self.time = time
        self.df = pd.DataFrame()
        # rmsd = self.rmsdOverTime()
        # self.rmsd = self.makeDataFrame(rmsd)
        # print(self.rmsd)
        # self.writeFile()

    def run(self, by='protein', selection='backbone'):
        if not self
    @property
    def first_frame(self):
        return self.traj[0]

    def squareDist(self, a1, a0):
        sqd = ((a1[0] - a0[0]) ** 2) + ((a1[1] - a0[1]) ** 2) + ((a1[2] - a0[2]) ** 2)
        return sqd

    def calcRMSD(self, frames, ref=None):
        if ref is None:
            if self.reference is not None:
                ref = self.reference
        rms = []
        for frame in frames:
            rms.append(np.sqrt((np.linalg.norm(frame-ref) / frame.shape[0])))
        return np.array(rms)
    
    def rmsdOverTimeOLD(self, frames = None, reference = None):
        if frames is None:
            rmsd = list(map(self.calcRMSD, (self.traj)))
        else:
            rmsd = list(map(self.calcRMSD, frames))
        # self.rmsd = self.makeDataFrame(rmsd)
        return rmsd

    def matrix(self, to_csv=None):
        df = pd.DataFrame(index=[i for i in range(0, len(self.traj.frames))], columns=[i for i in range(0, len(self.traj.frames))])
        i = 0
        frames = [np.array(frame) for frame in self.traj.frames]
        frames = np.array(frames)
        for frame in frames:
            self.reference = frame
            rms = np.array(self.rmsdOverTimeOLD())
            df.loc[i,:] = rms
            i += 1
        X = df.to_numpy(dtype='float64')
        if to_csv is not None:
            np.savetxt(to_csv, X, delimiter=',')
        return X


    def makeDataFrame(self, rmsd):
        df = pd.DataFrame(rmsd, index=[i for i in range(len(rmsd))])
        return df
        
    def writeFile(self, output='rmsd.csv'):
        self.df.to_csv(output, header='RMSD', index_label='Time')
        

    