import os
import sys

import numpy as np
import pandas as pd

if os.name == 'nt':
    sys.path.append('D:/Work/iapp/')
if os.name == 'posix':
    sys.path.append('/mnt/d/Work/iapp')

from pymd.structure.protein import Protein, Residue, Atom
from pymd.mdanalysis.xpm import XpmParser
from pymd.mdanalysis.multixpm import MultiXPM


class Contacts(XpmParser, MultiXPM):

    def __init__(self, xpm, structure, num_peptides, num_residues, igncaps=True, multi=True):
        
        if igncaps is True:
            self.protein = Protein(structure=structure, ignore=['ACE', 'NH2'])
        else:
            self.protein = Protein(structure=structure)
        self.multi = multi
        self.input = xpm
        self.num_peptides = num_peptides
        self.num_residues = num_residues
        self.igncaps = igncaps

    def calculateContacts(self):
        contacts = {}
        if self.multi is False:
            return self.calculateContactsFrame()
        values = self.getValueAssignment()
        i = 0
        for frame in self.getFrame():
            df = self.getMatrix(_input=frame, values=values)
            for index in df.index:
                if index not in contacts.keys():
                    contacts[index] = {}
                for column in df.columns:
                        distance = df.iloc[index, column]
                        calc = self.getContactType(self.protein.residues[index], self.protein.residues[column], distance)
                        contacts[index][column] = calc
                        # if calc is True:
                        #     print(self.protein.residues[index].id, self.protein.residues[index].chain, self.protein.residues[index].index, '||', self.protein.residues[column].id, self.protein.residues[column].chain, self.protein.residues[column].index)
            i += 1
            if i == 1:
                break

    def getContactType(self, res1, res2):
        hbond = False
        ca = False
        saltbridge = False
        pi = False

        if res1 is res2:
            return 0
        
        for a in res1.atoms:
            for b in res2.atoms:
                if HBond(a, b).isHBond is True:
                    return True
        return 0


        
    @staticmethod
    def distance(x,y):
        return np.sqrt((x.coordinates[0]-y.coordinates[0])**2 + (x.coordinates[1] - y.coordinates[1])**2 + (x.coordinates[2] - y.coordinates[2])**2)
    
    @staticmethod 
    def angle(x,y,z):

        a = np.array(x.coordinates) 
        b = np.array(y.coordinates)
        c = np.array(z.coordinates) 

        ba = a - b
        bc = c - b

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)

        return np.degrees(angle)



class HBond(Contacts):

    def __init__(self, atom1, atom2):
        if atom1.hbond is True:
            if ('N' in atom2.name) or ('O' in atom2.name):
                self.donor = atom1
                self.acceptor = atom2
            else:
                self.donor = None
                self.acceptor = None
        elif atom2.hbond is True:
            if ('N' in atom2.name) or ('O' in atom2.name):
                self.donor = atom2
                self.acceptor = atom1
            else:
                self.donor = None
                self.acceptor = None
        else:
            self.donor = None
            self.acceptor = None

    @property
    def isHBond(self):
        distance = False
        theta = False
        if (self.donor is None) or (self.acceptor is None):
            return False
        if self.distance(self.donor, self.acceptor) <= 3:
            distance = True
        if self.angle(self.donor.electronegative, self.donor, self.acceptor) >= 120:
            theta = True
        if (distance is True) and (theta is True):
            return True
        return False
        
        
        
        




c = Contacts(xpm='D:/Work/grant/MDmat_Frames/IAPP/1/IAPP1_MDmat_Frames.xpm', structure='D:/Work/iapp/systems/L500_PDB/Trajectories/IAPP/1/md_1500.gro', num_peptides=6, num_residues=10)
c.calculateContacts()
def angle(x,y,z):

    a = np.array(x) 
    b = np.array(y)
    c = np.array(z) 
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)

