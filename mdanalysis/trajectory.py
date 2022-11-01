import multiprocessing as mp
import mdtraj
from mdtraj.core.topology import Atom
import numpy as np
import os
import pandas as pd
import sys
sys.path.append('/work/cascades/kelsieking23/software/')
# sys.path.append("D:/Work/")
from pymd.mdanalysis.ndx import Ndx
from scipy.spatial import distance_matrix

class Trajectory:

    '''
    XtcParser takes a dumped XTC and performs various operations. 
    '''
    def __init__(self, xtc, gro, interval=(0,-1), skip=0, slice=tuple(), selection=None, peptides=None):
        self.xtc = xtc
        self.gro = gro
        if skip <= 1:
            traj = mdtraj.load(self.xtc, top=self.gro)
        else:
            traj = mdtraj.load(self.xtc, top=self.gro, stride=skip)
        self.full = traj
        top = traj.topology
        self.top = top
        if selection is not None:
            sele = top.select(selection)
            self.frames = traj._xyz[:,sele,:]
        else:
            self.frames = traj._xyz
        self.frames = self.frames[interval[0]:interval[1], :, :]
        self.time = traj._time[interval[0]:interval[1]]
        # if slice != tuple():
        #     self.frames = self.frames[slice[0]:slice[1], :, :]
        #     self.time = self.time[slice[0]:slice[1]]
        # self._traj = self.traj
        # self.top = self._traj._topology
        # self._top = self.top

        self.sele = selection

        if peptides is None:
            self.peptides = len(self.top._chains)
        else:
            self.peptides=peptides


    def coordinates(self):
        for frame in self.frames:
            yield frame
    
    def sample(self, skip=10):
        self.traj._xyz = self.frames[::skip]
        return self

    # @property
    # def frames(self):
    #     return self.traj._xyz
    
    def select(self, selection, ndxt=None):
        self.sele=selection
        sele = self._getSelection(selection, ndxt)
        newtraj = Trajectory(self.xtc, self.gro, sele)
        newtraj._selectCoords()
        newtraj._selectAtoms()
        newtraj._selectResidues()
        return newtraj

    def _getSelection(self, selection, ndxt=None):
        '''
        get selection
        '''
        ndx = Ndx(self.gro, peptides=self.peptides, ndxt=ndxt)
        sele = list(map(int, ndx.select(selection)))
        return [(i-1) for i in sele]

    def _selectCoords(self):
        x = list(map(self.mapSele, self.frames))
        self.traj._xyz = x


    def _selectAtoms(self):
        atoms = []
        i = 0
        for atom in self.top._atoms:
            if self.sele is not None:
                if i in self.sele:
                    atoms.append(atom)
            else:
                atoms.append(atom)
            i += 1
        self.top._atoms = atoms

    def _selectResidues(self):
        residues = []
        for residue in self.top._residues:
            for atom in self.top._atoms:
                if atom in residue._atoms:
                    residues.append(residue)
                    break
        self.top._residues = residues

    def mapSele(self, frame):
        return frame[self.sele]

    def getFrame(self, t):
        i = 0
        for time in self.time:
            if time == t:
                return self.frames[i]
            i += 1


# print(array.shape)
# nsamples, nx, ny = array.shape
# X = array.reshape((nsamples,nx*ny))
# df = pd.DataFrame(X)
# print(df.head())