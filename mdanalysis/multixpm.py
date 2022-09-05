import os
import argparse
from matplotlib import colors
from matplotlib import pyplot as plt
import pandas as pd

import sys
sys.path.append(os.getcwd())
from pymd.mdanalysis.xpm import XpmParser
from pymd.mdanalysis.mixedxpm import MixedParser


class MultiXPM(XpmParser):

    def __init__(self, _input, num_peptides, num_residues, num_residues_2=None, num_peptides_2=None, igncaps=True):
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
        if num_residues_2 is not None:
            self.num_residues_2 = num_residues_2
        if num_peptides_2 is not None:
            self.num_peptides_2 = num_peptides_2
        self.vmin = None
        self.vmax = None
        self.order = []

    def reader(self):
        f = open(self.input, 'r')
        for line in f:
            yield line
        f.close()

    def getFrame(self):
        contents = []
        first_line  = True
        for line in self.reader():
            if ('/* XPM */' in line) and (first_line is True):
                first_line  = False
            elif ('/* XPM */' in line) and (first_line is False):
                yield contents
                contents = []
            else:
                contents.append(line)

    def averageFrames(self):
        i = 0
        df = None
        values = self.getValueAssignment()
        m = MixedParser(_input=None, num_residues_1=self.num_residues, num_residues_2=self.num_residues_2, num_peptides_1=self.num_peptides, num_peptides_2=self.num_peptides_2, igncaps=self.igncaps)
        for frame in self.getFrame():
            if self.num_peptides_2 is None:
                if df is None:
                    df = self.getMatrix(_input=frame, values=values)
                else:
                    df = df + self.getMatrix(_input=frame, values=values)
                i += 1
            else:
                if df is None:
                    df = m.getMatrixMixed(_input=frame, values=values)
                else:
                    df = df + m.getMatrixMixed(_input=frame, values=values)
                i += 1
        u = df / i
        return  df / i

    def getDataFrames(self):
        values = self.getValueAssignment()
        for frame in self.getFrame():
            yield self.getMatrix(_input=frame, values=values) 




