import multiprocessing as mp
import mdtraj
import numpy as np
import os
import pandas as pd

import sys
if sys.platform == 'win32':
    sys.path.append('D:/Work/iapp/')
if sys.platform == 'linux':
    sys.path.append('/mnt/d/Work/iapp/')
# from pymd.structure.protein import Protein

class XtcParser:

    '''
    XtcParser takes a dumped XTC and performs various operations. 
    '''
    def __init__(self, xtc, gro=None):
        self.xtc = xtc
        self.gro = gro

        self.traj = mdtraj.load(self.xtc, top=self.gro)


    def coordinates(self):
        for frame in self.frames:
            yield frame

    @property
    def frames(self):
        return self.traj._xyz


    def framesOld(self):
        f = open(self.xtc, 'r')
        i = 0
        for line in f:
            if line.strip().startswith('natoms'):
                line_parts = line.strip().split()
                step = int(line_parts[3])
                time = step / 500000
            if line.strip().startswith('x['):
                atom_index = int(line.strip()[2:7].strip())
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
                data = {atom_index:coordinates}
                yield data, time
        f.close()

    def rawXTC(self):
        f = open(self.xtc, 'rb')
        i = 0
        for line in f:
            print(line)
            i += 1
            if i == 20:
                break

    
    def xtcReader(self):
        f = open(self.xtc, 'r')
        i = 0
        for line in f:
            if line.strip().startswith('natoms'):
                line_parts = line.strip().split()
                step = int(line_parts[3])
                time = step / 500000
            if line.strip().startswith('x['):
                atom_index = int(line.strip()[2:7].strip())
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
                data = {atom_index:coordinates}
                yield data, time
        f.close()

    def getPositionMatrices(self):
        '''
        Returns a list of dataframes. The index of each element in the list corresponds with the atom index. Each dataframe contains the x, y, and z coordinates (indeces) for each frame in the .xtc (columns).
        '''
        f = open(self.xtc, 'r')
        i = 0
        frames = 0
        positions = [{} for i in range(len(self.atoms))]
        for line in f:
            if line.strip().startswith('natoms'):
                line_parts = line.strip().split()
                step = int(line_parts[3])
                time = step / 500000
                frames += 1
            if line.strip().startswith('x['):
                atom_index = int(line.strip()[2:7].strip())
                positions[atom_index][time] = {}
                coordinates_split = line.strip().split('=')[-1].split(',')
                coordinates = []
                dimensions = ['x', 'y', 'z']
                k = 0
                for split_coord in coordinates_split:
                    string = ''
                    for char in split_coord:
                        if (char == '{') or (char == ',') or (char == '}') or (char == ' '):
                            continue
                        else:
                            string = string + char
                    coordinate = float(string)*10
                    positions[atom_index][time][dimensions[k]] = coordinate
                    k += 1
                i +=1
        f.close()
        self.frames = frames
        dfs = []
        for pos in positions:
            df = pd.DataFrame(pos)
            dfs.append(df)
        return positions

    def getFrames(self):
        f = open(self.xtc, 'r')
        frames = 0
        for line in f:
            if line.strip().startswith('natoms'):
                frames += 1
        f.close()
        return frames

    def getAtomInteractions(self):
        '''
        Get total number of interactions over time from .xtc. Utilizes multiprocessing; will run processes equal to number of cpu cores at a given time. Returns a dataframe.
        '''
        import time
        start = time.time()
        matrices = self.getPositionMatrices()
        output = mp.Queue()


        def calculate(ref, index):
            interactions = {}
            if index not in interactions.keys():
                interactions[index] = {}
            for column in ref.columns:
                current = ref.loc[:,column]
                i = 0
                for matrix in matrices:
                    if i not in interactions[index].keys():
                        interactions[index][i] = 0
                    if i == index:
                        interactions[index][i] = 0
                    else:
                        m = matrix.loc[:,column]
                        d = np.sqrt(np.sum([(a-b)*(a-b) for a, b in zip(current, m)])) 
                        if d <= 3:
                            interactions[index][i] += 1
                    i += 1
            output.put(interactions)
            return interactions

        interactions_master = {}
        cores = mp.cpu_count()
        jobs = []
        order = {}
        atoms = []
        for i in range(0, len(self.atoms)):
            jobs.append(mp.Process(target=calculate, args=(matrices[i], i)))
            atoms.append(i)
            if i % cores == 0:
                for job in jobs:
                    job.start()
                for job in jobs:
                    job.join()
                for job in jobs:
                    ints = output.get()
                    for key, value in ints.items():
                        order[key] = value
                for index in atoms:
                    dic = order[index]
                    interactions_master[index] = {}
                    for key, value in dic.items():
                        interactions_master[index][key] = value
                order = {}
                jobs = []
                atoms = []
        if jobs != []:
            for job in jobs:
                job.start()
            for job in jobs:
                job.join()
            for job in jobs:
                ints = output.get()
                for key, value in ints.items():
                    order[key] = value
            for index in atoms:
                dic = order[index]
                interactions_master[index] = {}
                for key, value in dic.items():
                    interactions_master[index][key] = value
        end = time.time()
        elapsed = (end - start)/60
        print(elapsed)
        df = pd.DataFrame(interactions_master)
        return df

    def getResidueInteractionProbabilities(self, csv=None):
        '''
        Transforms atom interaction matrix (from self.getAtomInteractions()) into a residue interaction matrix. Returns normalized probabilities - each dataframe value equals: ((total # interactions over time)/(number of atoms in both residues))/(# frames). 
        Returns dataframe.
        Arguments:
        ** csv (string): optional. if you have a csv written of the output from self.getAtomInteractions(), you can pass this to save time and bypass this calculation. 
        '''
        # get interaction data
        if csv is None:
            matrix = self.getAtomInteractions()
        else:
            matrix = pd.read_csv(csv, header=0, index_col=0)
        
        # get frames
        if self.frames is None:
            self.frames = self.getFrames()

        residue_interactions = {}

        # turn atom data into residue data
        # residue a
        for i in range(0, len(self.residues)):
            resA = self.residues[i]
            if i not in residue_interactions.keys():
                residue_interactions[i] = {}
            resA_indeces = []
            for atom in resA.atoms:
                resA_indeces.append(atom.index)
            
            # residue b
            for k in range(0, len(self.residues)):
                if k not in residue_interactions[i].keys():
                    residue_interactions[i][k] = 0
                resB = self.residues[k]
                resB_indeces = []
                for atom in resB.atoms:
                    resB_indeces.append(atom.index)

                if i != k:
                    num_interactions = 0
                    for resA_atom_index in resA_indeces:
                        for resB_atom_index in resB_indeces:
                            num_interactions += matrix.iloc[resA_atom_index, resB_atom_index]
                    residue_interactions[i][k] = (num_interactions / (len(resB_indeces) + len(resA_indeces))) / self.frames
                else:
                    residue_interactions[i][k] = 0

        df = pd.DataFrame(residue_interactions)
        return df

    def getResidueInteractions(self, csv=None):
        # get interaction data
        if csv is None:
            matrix = self.getAtomInteractions()
        else:
            matrix = pd.read_csv(csv, header=0, index_col=0)
        
        residue_interactions = {}

        for i in range(0, len(self.residues)):
            resA = self.residues[i]
            if i not in residue_interactions.keys():
                residue_interactions[i] = 0
            resA_indeces = []
            for atom in resA.atoms:
                resA_indeces.append(atom.index)

            # residue b
            for k in range(0, len(self.residues)):
                resB = self.residues[k]
                resB_indeces = []
                for atom in resB.atoms:
                    resB_indeces.append(atom.index)
 
                _break = False
                if i != k:
                    for resA_atom_index in resA_indeces:
                        for resB_atom_index in resB_indeces:
                            if matrix.iloc[resA_atom_index, resB_atom_index] != 0:
                                residue_interactions[i] += 1
                                _break = True
                                break
                        if _break is True:
                            break
        return residue_interactions
            

# x = XtcParser('cat_pbc_pro.xtc', gro='md_500_600_nowat.gro')
# i = 0
# for coord in x.coordinates():
#     print(coord[0])
#     i += 1
#     if i == 5:
#         break
# x = XtcParser(xtc='D:/Work/iapp/tests/cat_pbc_pro.xtc')
# x.rawXTC()
# x = XtcParser(xtc='D:/Work/iapp/systems/L500_PDB/Trajectories/IAPP/1/dump_centered_IAPP1_1500_2000ns.xtc', gro='D:/Work/iapp/systems/L500_PDB/Trajectories/IAPP/1/md_1500.gro')
# i = 0
# for frame, time in x.xtcReader():
#     print(frame)
#     print(time)
#     i += 1
#     break
