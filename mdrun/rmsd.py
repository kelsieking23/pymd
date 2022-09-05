import math
from operator import index
import mdtraj
import numpy as np
import os
import pandas as pd
import sys
# sys.path.append('/work/cascades/kelsieking23/software/')
sys.path.append('D:/Work')
from pymd.mdanalysis.trajectory import Trajectory


class RMSD:

    def __init__(self, traj, reference=None, cutoff=0, binary=False):
        if traj is not None:
            self.traj = traj
        if reference is None:
            self.reference = self.first_frame
        else:
            self.reference = reference
        self.time = self.traj.time
        self.cutoff = cutoff
        self.binary = binary
        # rmsd = self.rmsdOverTime()
        # self.rmsd = self.makeDataFrame(rmsd)
        # print(self.rmsd)
        # self.writeFile()

    @property
    def first_frame(self):
        return self.traj.frames[0]

    def squareDist(self, a1, a0):
        sqd = ((a1[0] - a0[0]) ** 2) + ((a1[1] - a0[1]) ** 2) + ((a1[2] - a0[2]) ** 2)
        return sqd

    def calcRMSD(self, frame, reference=None):
        total_sqd = 0
        if reference is None:
            reference = self.reference
        rmsd = np.linalg.norm(frame-reference) / frame.shape[0]
        if self.binary:
            if rmsd > self.cutoff:
                return rmsd
            else:
                return 0
        return rmsd
    
    def rmsdOverTime(self, frames = None, reference = None):
        if frames is None:
            rmsd = list(map(self.calcRMSD, self.traj.coordinates()))
        else:
            rmsd = list(map(self.calcRMSD, frames))
        self.rmsd = self.makeDataFrame(rmsd)
        return self.rmsd

    def matrix(self, to_csv=None):
        df = pd.DataFrame(index=[i for i in range(0, len(self.traj.frames))], columns=[i for i in range(0, len(self.traj.frames))])
        i = 0
        frames = [np.array(frame) for frame in self.traj.frames]
        frames = np.array(frames)
        for frame in frames:
            self.reference = frame
            rms = np.array(self.rmsdOverTime())
            df.loc[i,:] = rms
            i += 1
        X = df.to_numpy(dtype='float64')
        if to_csv is not None:
            np.savetxt(to_csv, X, delimiter=',')
        return X


    def makeDataFrame(self, rmsd):
        df = pd.DataFrame(rmsd, index=self.time)
        return df
        
    def writeFile(self, output='rmsd.csv'):
        self.rmsd.to_csv(output, header='RMSD', index_label='Time')
        

    