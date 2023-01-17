import os
import sys

from collections.abc import Iterable
import numpy as np
import mdtraj
import pandas as pd

class Peptide:
    '''
    class name for now might be misnomer, but designing this to be a replacement for Protein class
    for use with a trajectory. 
    '''

    def __init__(self, structure, xtc=None, stride=0, selection='all', ignore=[], name=''):
        self.structure = structure
        self.selection = selection
        self.xtc = xtc
        self.ignore = ignore
        self.topology = self._topology()
        self.top = self.topology
        self._traj = self._iterload(stride)
        self.properties = Properties(self) # a place to store all the random stuff calculated so it doesn't get cluttered here
        self.data = np.array([])
        self.name = name
        self.connections = {}
    # Initilization functions
    def _topology(self):
        top = mdtraj.load(self.structure).topology
        if isinstance(self.selection, str):
            sele = top.select(self.selection)
        elif isinstance(self.selection, Iterable):
            sele = self.selection
        else:
            raise ValueError('selection must be string or iterable')
        return self._chain_topology_from_subset(top, sele)
    
    def _chain_topology_from_subset(self, top, sele):
        '''
        this isa modified mdtraj function that doesn't reset verything to index 0 if you just
        pick out a single chain. 
        '''
        new = mdtraj.Topology()
        old_atom_to_new_atom = {}
        selection = self.selection
        chain_index = int(selection.split()[-1])
        newChain = new.add_chain_custom(chain_index)
        for chain in top._chains:
            if chain.index != chain_index:
                continue
            for residue in chain._residues:
                resSeq = getattr(residue, 'resSeq', None) or residue.index
                newResidue = None
                old_index = residue.index
                for atom in residue._atoms:
                    if atom.index in sele:
                        try:  # OpenMM Topology objects don't have serial attributes, so we have to check first.
                            serial = atom.serial
                        except AttributeError:
                            serial = None
                        if newResidue is None:
                            newResidue = new.add_residue_custom_index(residue.name, newChain,
                                                                    old_index, resSeq, residue.segment_id)
                        newAtom = new.add_atom_custom_index(atom.name, atom.element, atom.index,
                                    newResidue, serial=serial)
                        old_atom_to_new_atom[atom] = newAtom
        
        bondsiter = top.bonds
        if not hasattr(bondsiter, '__iter__'):
            bondsiter = bondsiter()
        for bond in bondsiter:
            try:
                atom1, atom2 = bond
                new.add_bond(old_atom_to_new_atom[atom1],
                                        old_atom_to_new_atom[atom2],
                                        type=bond.type,
                                        order=bond.order)
            except KeyError:
                pass
                # we only put bonds into the new topology if both of their partners
                # were indexed and thus HAVE a new atom
        return new

    @property
    def trajectory(self):
        return self._traj
    
    @property
    def traj(self):
        '''
        alias for trajectory
        '''
        return self.trajectory
    
    def _iterload(self, stride):
        if self.xtc is None:
            return None
        return mdtraj.iterload(self.xtc, chunk=100, atom_indices=self.atom_indeces, stride=stride, top=self.topology)
    
    @property
    def sele(self):
        return self.top.select(self.selection)
    @property
    def atom_indeces(self):
        return np.array([atom.index for atom in self.topology.atoms if atom.index in self.sele])
    
    @property
    def residue_indeces(self):
        res_indeces = []
        for atom in self.topology.atoms:
            if atom.index in self.sele:
                if atom.residue.index not in res_indeces:
                    res_indeces.append(atom.residue.index)
        return np.array(res_indeces)
        
    
    @property
    def residue_ids(self):
        return np.array(['{}{}'.format(self.top.residue(i).name, self.top.residue(i).resSeq) for i in self.residue_indeces])
    
    # Property functions
    def com(self):
        if self.traj is None:
            raise ValueError('No trajectory loaded!')
        _coms = []
        for chunk in self.traj:
            _coms.append(mdtraj.compute_center_of_mass(chunk))
        _coms = np.concatenate(_coms)
        setattr(self.properties, 'com', _coms)
        self.data = _coms
        return _coms

class Properties:

    def __init__(self, peptide):
        self.parent = peptide
    


# i = 2
# PATH = 'D:/Work/Projects/Abeta-barto/concat/{}-mer'.format(i)
# xtc = os.path.join(PATH, 'concat.pbc.nowat.stride100.xtc')
# pdb = os.path.join(PATH, 'top.nowat.pdb')

# p = Peptide(pdb, xtc, selection='chainid 0')
# print(p.residue_ids)
# atoms = [atom.index for atom in p.top.atoms]
# print(atoms)
# print(p.top.to_dataframe())
# p2 = Peptide(pdb, xtc, selection='chainid 1')
# atoms = [atom.index for atom in p2.top.atoms]
# print(p2.top.to_dataframe())

# all = Peptide(pdb, xtc, selection='all')
# atoms = [atom.index for atom in all.top.atoms]
# p1 = all.top.select('chainid 1')
# p2 = all.top.subset(p1)
# print([atom.index for atom in p2.atoms])
# top = mdtraj.load(pdb).topology
# sele = top.select('chainid 1')

