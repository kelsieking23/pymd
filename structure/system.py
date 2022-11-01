####################################################################################################################
# system.py                                                                                                        #
# Written by Kelsie M. King                                                                                        #
# Last Edited 06/03/2020                                                                                           #
# System contains the System class, which is the main workhorse of the module. This class contains all functions   #
# related to data analysis of the protein and ligand.                                                              #
####################################################################################################################

import numpy as np
import pandas as pd
import statistics
import sys
import os 

# sys.path.append('D:/Work')
# from grapher import Grapher
from pymd.structure.protein import Protein
from pymd.structure.ligand import Ligand
from pymd import utilities

class System:
    '''
    Performs operations on Protein and Ligand classes.
    Arguments:
    ** protein: string. Path to protein structure file. Creates instance of Protein class.
    ** ligands (optional): list or string. Path(s) to ligand structure file(s). Creates instance(s) of Ligand class. ### needs to be modified to make true - only works for list ATM
    ** residues (optional): list. Residue IDs to get coordinates of & perform operations with. ID syntax is <residue_name><residue_number> (ex: SER144, HIS41, CYS145)
    ** covalent_ids (optional): list or string. If the ligand(s) is/are covalent, or otherwise in the protein structure file, this indicates which lines of data are associated with the 
        ligand(s) in the PDB file. ID syntax is <residue_name><residue_number> (ex: FHR405 (see PDBID: 6LZE for reference))
    ** model (optional): list or string. Model coordinates to use as ligand coordinates if multiple models are in the ligand structure file. Default is MODEL 1. If list, must be passed 
        in the same order as the ligands passed. For example, models[0] must be associated with ligands[0]. 
    ** names (optional): string or list or dict. Name(s) to associate with ligand(s). If a list, names must be passed in the same order as the ligands passed. For example, names[0] must be
        associated with ligands[0]. If dict, key must = path to ligand file, value must = name. ## need to modify to allow dict.
    ** alt (optional): list. To map ligand names to alternative names.
    ** ignh (optional): boolean. Specifies if hydrogen should be ignored when collecting coordinates. Default is false
    '''
    def __init__(self, protein, ligands=None, residues=None, covalent_ids=None, models=None, names=None, alt=None, ignh=False):
        self.protein = Protein(structure=protein, residues=residues, covalent_ids=covalent_ids)
        if residues is not None:
            self.residues = residues
        else:
            self.residues = self.protein.ids
        if isinstance(models, str):
            models = [models] * len(ligands)
        elif isinstance(models, list):
            pass
        else:
            models = ['MODEL 1'] * len(ligands)
        if ligands is not None:
            self.ligands = []
            if covalent_ids is None:
                if names is not None:
                    if isinstance(names, list):
                        for i in range(len(ligands)):
                            if alt is not None:
                                self.ligands.append(Ligand(structure=ligands[i], alt=alt[i], model=models[i], name=names[i], ignh=ignh))
                            else:
                                self.ligands.append(Ligand(structure=ligands[i], model=models[i], name=names[i], ignh=ignh))
                    if isinstance(names, dict):
                        for key in names.keys():
                            name = names[key]
                            self.ligands.append(Ligand(structure=key, model=models[key], name=name, ignh=ignh))
                        loaded = []
                        for ligand in self.ligands:
                            loaded.append(ligand.structure)
                        i = 0
                        for ligand in ligands:
                            if ligand not in loaded:
                                self.ligands.append(Ligand(structure=ligand, model=models[i], ignh=ignh))
                                i += 1
                else:
                    if isinstance(ligands, list):
                        i = 0
                        for ligand in ligands:
                            self.ligands.append(Ligand(structure=ligand, model=models[i], name=os.path.splitext(os.path.basename(ligand))[0], ignh=ignh))
                            i += 1
                    if isinstance(ligands, str):
                        self.ligands.append(Ligand(structure=ligands, model=models[0], ignh=ignh))
            else:
                self.ligands = self.protein.ligands
        
        self.names = []
        if ligands is not None:
            for ligand in self.ligands:
                if ligand.name is None:
                    ligand.name = ''.join(os.path.basename(ligand.structure).split('.')[:-1])
                self.names.append(ligand.name)
        else:
            if names is not None:
                for name in names:
                    self.names.append(name)
        self.alts = []
        if alt is not None:
            for altname in alt:
                self.alts.append(altname)
                
    ######################################################################################
    ######################## Interaction Matrices Functions ##############################
    ######################################################################################
    @staticmethod
    def euclideanDistance(x, y):
        try:
            assert len(x) == len(y)
        except:
            raise ValueError('X and Y must be of equal length')
        return np.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2)

    @staticmethod
    def rmsd(x, y):
        try:
            assert len(x) == len(y)
        except:
            raise ValueError('X and Y must be of equal length')
        total_squared_distance = 0
        for i in range(0, len(x)):
            total_squared_distance += (x[i][0] - y[i][0])**2 + (x[i][1] - y[i][1])**2 + (x[i][2] - y[i][2])**2
        return np.sqrt(total_squared_distance / len(x))   

    def interactions(self, distance=5, norm=False, heavy=False, to_csv=None):
        '''
        Create basic interaction matrix. Interaction for a given ligand atom is defined as an atom within 5 Angstroms of the residue center of mass. 
        Arguments:
        ** distance (optional): int or float. Interaction cutoff distance (in Angstroms). Default is 5 Angstroms. 
        ** norm (optional): boolean. Specify if interactions should be normalized by number of atoms in ligand.
           If true, will normalize by compound size. Default is False. 
        ** heavy (optional): boolean. If normalizing by compound size (if norm = True), specify if interactions should be normalized by heavy atom count or all atom count.
           If true, will normalize by heavy atom count. Default is False. 
        ** to_csv (optional): string. Path to .csv output file. Will only output a .csv if an argument is passed. 
        Returns: DataFrame of interaction matrix. 
        '''
        interactions = {}
        for ligand in self.ligands:
            interactions[ligand.name] = {}
            num_interactions = 0
            for res_id in self.protein.ids:
                res_com = self.protein.coms[res_id]
                for atom in ligand.atoms:
                    d = np.sqrt((atom.coordinates[0] - res_com[0])**2 + (atom.coordinates[1] - res_com[1])**2 + (atom.coordinates[2] - res_com[2])**2)
                    if d <= distance:
                        num_interactions += 1
                interactions[ligand.name][res_id] = num_interactions
                num_interactions = 0
            if norm == True:
                atom_count = ligand.atomCount(heavy=heavy)
                for res_id in self.protein.ids:
                    interactions[ligand.name][res_id] = interactions[ligand.name][res_id] / atom_count
        df = pd.DataFrame(interactions)
        if to_csv is not None:
            df.to_csv(to_csv)
        return df

    def interactionsTotal(self, distance=5, norm=False, heavy=False, to_csv=None):
        '''
        Create a pandas Series with totaled interactions. Interaction for a given ligand atom is defined as an atom within 5 Angstroms of the residue center of mass. 
        Arguments:
        ** distance (optional): int or float. Interaction cutoff distance (in Angstroms). Default is 5 Angstroms. 
        ** norm (optional): boolean. Specify if interactions should be normalized by number of atoms in ligand.
           If true, will normalize by compound size. Default is False. 
        ** heavy (optional): boolean. If normalizing by compound size (if norm = True), specify if interactions should be normalized by heavy atom count or all atom count.
           If true, will normalize by heavy atom count. Default is False. 
        ** to_csv (optional): string. Path to .csv output file. Will only output a .csv if an argument is passed. 
        Returns: DataFrame of interaction matrix. 
        '''
        interactions = {}
        for ligand in self.ligands:
            interactions[ligand.name] = {}
            interactions[ligand.name]['num_interactions'] = 0
            num_interactions = 0
            for res_id in self.protein.ids:
                res_com = self.protein.coms[res_id]
                for coord in ligand.coordinates:
                    d = np.sqrt((coord[0] - res_com[0])**2 + (coord[1] - res_com[1])**2 + (coord[2] - res_com[2])**2)
                    if d <= distance:
                        num_interactions += 1
            interactions[ligand.name]['num_interactions'] = num_interactions
            num_interactions = 0
            if norm == True:
                atom_count = ligand.atomCount(heavy=heavy)
                interactions[ligand.name]['num_interactions'] = interactions[ligand.name]['num_interactions'] / atom_count
        df = pd.DataFrame(interactions)
        if to_csv is not None:
            df.to_csv(to_csv)
        return df
    
    def minimumDistance(self, to_csv=None):
        '''
        Creates minimum distance matrix. Distance is measured as the distance between a given ligand atom and the center of mass of a given residue. 
        Arguments:
        ** to_csv (optional): string. Path to .csv file output. Will only output if argument is passed. 
        Returns: DataFrame of minimum distance matrix. 
        '''
        mindist = {}
        for ligand in self.ligands:
            mindist[ligand.name] = {}
            for res_id in self.protein.ids:
                mindist[ligand.name][res_id] = None
                res_com = self.protein.coms[res_id]
                for coord in ligand.coordinates:
                    d = np.sqrt((coord[0] - res_com[0])**2 + (coord[1] - res_com[1])**2 + (coord[2] - res_com[2])**2)
                    if mindist[ligand.name][res_id] == None:
                        mindist[ligand.name][res_id] = d
                    if d < mindist[ligand.name][res_id]:
                        mindist[ligand.name][res_id] = d
                    else:
                        continue
        df = pd.DataFrame(mindist)
        if to_csv is not None:
            df.to_csv(to_csv)
        return df

    def minimumDistanceByAtom(self):
        '''
        Gets minimum distance by atom. Creates an attribute of ligand containing minimum distance matrix. 
        '''
        for ligand in self.ligands:
            mindist = {}
            for residue in self.protein.residues:
                mindist[residue.id] = {}
                d = None
                protein_atom_name = None
                ligand_atom_name = None
                for protein_atom in residue.atoms:
                    for ligand_atom in ligand.atoms:
                        _d = self.euclideanDistance(protein_atom.coordinates, ligand_atom.coordinates)
                        if d is None:
                            d = _d
                            protein_atom_name = protein_atom.name
                            ligand_atom_name = ligand_atom.name
                        if (d > _d):
                            d = _d
                            protein_atom_name = protein_atom.name
                            ligand_atom_name = ligand_atom.name
                mindist[residue.id]['distance'] = d 
                mindist[residue.id]['protein_atom'] = protein_atom_name
                mindist[residue.id]['ligand_atom'] = ligand_atom_name
            df = pd.DataFrame(mindist).T
            ligand.minimum_distance = df

    def averageDistance(self, to_csv=None):
        '''
        Creates average distance matrix. Distance is measured as the distance between a given ligand atom and the center of mass of a given residue. 
        Arguments:
        ** to_csv (optional): string. Path to .csv file output. Will only output if argument is passed. 
        Returns: DataFrame of minimum distance matrix. 
        '''
        average_dist = {}
        for ligand in self.ligands:
            average_dist[ligand.name] = {}
            for res_id in self.protein.ids:
                average_dist[ligand.name][res_id] = []
                res_com = self.protein.coms[res_id]
                for coord in ligand.coordinates:
                    d = np.sqrt((coord[0] - res_com[0])**2 + (coord[1] - res_com[1])**2 + (coord[2] - res_com[2])**2)
                    average_dist[ligand.name][res_id].append(d)
                average_dist[ligand.name][res_id] = sum(average_dist[ligand.name][res_id]) / len(average_dist[ligand.name][res_id])
        df = pd.DataFrame(average_dist)
        if to_csv is not None:
            df.to_csv(to_csv)
        return df

    def interactionDistanceRatio(self, distance=5, method='average', to_csv=None):
        '''
        Calculates the ratio between number of interactions to a given residue and the average distance to a given residue. Creates an interaction matrix.  
        Distance is measured as the distance between a given ligand atom and the center of mass of a given residue. 
        Arguments:
        ** distance (optional): int or float. Specifies the interaciton distance cutoff. Default is 5 Angstroms. 
        ** method (optional): string. Can be 'average' or 'minimum'. Specifies if the ratio should be taken by minimum distance or average distance to a given residue.
           Default is 'average'. 
        ** to_csv (optional): string. Path to .csv file output. Will only output if argument is passed. 
        Returns: DataFrame of minimum distance matrix. 
        '''
        interactions = self.interactions(distance=distance, norm=False, heavy=False, to_csv=None)
        if method == 'average':
            distances = self.averageDistance(to_csv=None)
        if method == 'minimum':
            distances = self.minimumDistance(to_csv=None)
        df = interactions / distances
        if to_csv is not None:
            df.to_csv(to_csv)  
        return df

    def residueInteractions(self, distance=5, to_csv=None):
        interactions = {}
        for res_id in self.protein.ids:
            residue_coordinates = self.protein.coordinates[res_id]
            interactions[res_id] = {}
            for ligand in self.ligands:
                interactions[res_id][ligand.name] = 0
                ligand_coordinates = ligand.coordinates
                for residue_coordinate in residue_coordinates:
                    if isinstance(ligand_coordinates, list):
                        for ligand_coordinate in ligand_coordinates:
                            d = np.sqrt((ligand_coordinate[0] - residue_coordinate[0])**2 + (ligand_coordinate[1] - residue_coordinate[1])**2 + (ligand_coordinate[2] - residue_coordinate[2])**2)
                            if d <= 5:
                                # if res_id == 'PHE140':
                                #     print(ligand.name)
                                #     print(ligand_coordinate)
                                #     print(residue_coordinate)
                                interactions[res_id][ligand.name] = 1
                    if isinstance(ligand_coordinates, dict):
                        for atom_id in ligand_coordinates:
                            ligand_coordinate = ligand_coordinates[atom_id]
                            d = np.sqrt((ligand_coordinate[0] - residue_coordinate[0])**2 + (ligand_coordinate[1] - residue_coordinate[1])**2 + (ligand_coordinate[2] - residue_coordinate[2])**2)
                            if d <= 5:
                                # if res_id == 'PHE140':
                                #     print(ligand.name)
                                #     print(ligand_coordinate)
                                #     print(residue_coordinate)
                                interactions[res_id][ligand.name] = 1
        df = pd.DataFrame(interactions)
        df['sum'] = df.sum(axis=1)
        num_residues = len(self.protein.ids)
        num_ligands = len(self.ligands)
        percents = []
        for i in range(num_ligands):
            value = df.iloc[i, -1]
            percent = (value / num_residues) * 100
            percents.append(percent)
        df['%'] = percents
        if to_csv is not None:
            df.to_csv(to_csv)
        return df

    def interactionFrequency(self, distance=5, heavy=True, to_csv=None):
        interactions = {}
        for res_id in self.protein.ids:
            residue_coordinates = self.protein.coordinates[res_id]
            interactions[res_id] = {}
            for ligand in self.ligands:
                interactions[res_id][ligand.name] = 0
                for residue_coordinate in residue_coordinates:
                    for atom in ligand.atoms:
                        if (heavy is True) and (atom.type == 'H'):
                            continue
                        d = np.sqrt((atom.coordinates[0] - residue_coordinate[0])**2 + (atom.coordinates[1] - residue_coordinate[1])**2 + (atom.coordinates[2] - residue_coordinate[2])**2)
                        if d <= 5:
                            interactions[res_id][ligand.name] +=1
                atom_count = ligand.atomCount(heavy=heavy)
                interactions[res_id][ligand.name] = interactions[res_id][ligand.name] / atom_count
        df = pd.DataFrame(interactions)
        if to_csv is not None:
            df.to_csv(to_csv)
        return df


    ######################################################################################
    ################################# Vina Energy & LE ###################################
    ###################################################################################### 

    def vinaEnergy(self):
        energies = {}
        for ligand in self.ligands:
            energies[ligand.name] = {}
            f = open(ligand.log, 'r')
            contents = []
            for line in f:
                contents.append(line)
            energy_lines = contents[-10:-1]
            f.close()
            vina_energies = []
            for line in energy_lines:
                line_parts = line.split()
                energy = float(line_parts[1])
                vina_energies.append(energy)
            average = sum(vina_energies) / len(vina_energies)
            energies[ligand.name]['average_vina'] = average
            stdev = statistics.stdev(vina_energies)
            energies[ligand.name]['stdev_vina'] = stdev
            lowest_energy = vina_energies[0]
            energies[ligand.name]['lowest_energy_vina'] = lowest_energy
            atom_count = ligand.atomCount(heavy=True)
            le = lowest_energy / atom_count
            energies[ligand.name]['LE_vina'] = le
        df = pd.DataFrame(energies)
        return df
            
    
    ######################################################################################
    ################################# Schrod Functions ###################################
    ###################################################################################### 
    def mmgbsa(self, csv, to_csv=None):
        '''
        Gets MM/GBSA and ligand efficiency from Schrodinger .csv output. 
        Arguments:
        ** csv: string. Path to .csv file containing MM/GBSA information for ligands. 
        Returns: DataFrame
        ** to_csv (optional): string. Path to .csv file output. Will only output if argument is passed. 
        '''
        mmgbsa = {}
        f = open(csv, 'r')
        i = 0
        for line in f:
            if i == 0:
                i += 1
                continue
            else:
                line_parts = line.split(',')
                name = line_parts[0]
                dG = float(line_parts[1])
                for ligand in self.ligands:
                    if name == ligand.name:
                        mmgbsa[ligand.name] = {}
                        atoms = ligand.atomCount(heavy=True)
                        ligand_efficiency = dG / atoms
                        mmgbsa[ligand.name]['mmgbsa'] = dG
                        mmgbsa[ligand.name]['LE_mmgbsa'] = ligand_efficiency
        f.close()
        df = pd.DataFrame(mmgbsa)
        if to_csv is not None:
            df.to_csv(to_csv)  
        return df

    @staticmethod
    def getQikpropData(out):
        all_contents = []
        info = []
        f = open(out, 'r')
        i = 0
        for line in f:
            line_parts = line.split()
            if 'QikProp' in line_parts:
                if i == 0:
                    i += 1
                    continue
                else:
                    all_contents.append(info)
                    info = []
            else:
                info.append(line)
        all_contents.append(info)
        f.close()
        properties = {}
        drug_standards = {}
        for report in all_contents:
            drug_id = report[0].strip()
            print(drug_id)
            if drug_id == '':
                drug_id = report[2].strip()
            # if drug_id not in self.names:
            #     continue
            properties[drug_id] = {}
            for line in report:
                i = 0
                line_parts_unstripped = line.split()
                line_parts = list(map(str.strip, line_parts_unstripped))
                if line_parts != []:
                    if 'Solute' in line_parts[0]:
                        if 'Molecular' in line_parts[1]:
                            if 'Weight' in line_parts[2]:
                                mw = line_parts[4]
                                mw = float(line_parts[4])
                                properties[drug_id]['MW'] = mw
                                if 'MW' not in drug_standards.keys():
                                    low_limit = float(line_parts[6])
                                    if '*' not in line_parts[8]:
                                        up_limit = float(line_parts[8][:-1])
                                    else:
                                        up_limit = float(line_parts[8][:-2])
                                    limits = (low_limit, up_limit)
                                    drug_standards['MW'] = limits
                            elif 'Volume' in line_parts[2]:
                                vol = float(line_parts[4])
                                properties[drug_id]['MV'] = vol
                                if 'MV' not in drug_standards.keys():
                                        low_limit = float(line_parts[6])
                                        up_limit = float(line_parts[7][1:-2])
                                        limits = (low_limit, up_limit)
                                        drug_standards['MV'] = limits
                            else:
                                continue
                        elif 'Electron' in line_parts[1]:
                            electron_affinity = float(line_parts[5])
                            properties[drug_id]['eV'] = electron_affinity
                            if 'eV' not in drug_standards.keys():
                                low_limit = float(line_parts[7])
                                up_limit = float(line_parts[9][:-1])
                                drug_standards['eV'] = (low_limit, up_limit)
                        elif 'SASA' in line_parts[2]:
                            if 'Total' in line_parts[1]:
                                total_sasa = float(line_parts[4])
                                properties[drug_id]['TotalSASA'] = total_sasa
                                if 'TotalSASA' not in drug_standards.keys():
                                    low_limit = float(line_parts[6])
                                    up_limit = line_parts[7]
                                    if '*' in up_limit:
                                        up_limit = float(line_parts[7][1:-2])
                                    else:
                                        up_limit = float(line_parts[7][1:-1])
                                    drug_standards['TotalSASA'] = (low_limit, up_limit)
                            elif 'Hydrophobic' in line_parts[1]:
                                hydrophobic_sasa = float(line_parts[4])
                                properties[drug_id]['HydrophobicSASA'] = hydrophobic_sasa
                                if 'HydrophobicSASA' not in drug_standards.keys():
                                    low_limit = float(line_parts[6])
                                    up_limit = line_parts[8]
                                    if '*' in up_limit:
                                        up_limit = float(line_parts[8][:-2])
                                    else:
                                        up_limit = float(line_parts[8][:-1])
                                    drug_standards['HydrophobicSASA'] = (low_limit, up_limit)
                            elif 'Hydrophilic' in line_parts[1]:
                                hydrophilic_sasa = float(line_parts[4])
                                properties[drug_id]['HydrophilicSASA'] = hydrophilic_sasa
                                if 'HydrophilicSASA' not in drug_standards.keys():
                                    low_limit = float(line_parts[6])
                                    up_limit = line_parts[8]
                                    if '*' in up_limit:
                                        up_limit = float(line_parts[8][:-2])
                                    else:
                                        up_limit = float(line_parts[8][:-1])
                                    drug_standards['HydrophilicSASA'] = (low_limit, up_limit)
                            elif 'Carbon' in line_parts[1]:
                                carbon_pi_sasa = float(line_parts[5])
                                properties[drug_id]['CarbonPiSASA'] = carbon_pi_sasa
                                if 'CarbonPiSASA' not in drug_standards.keys():
                                    low_limit = float(line_parts[7])
                                    up_limit = line_parts[9]
                                    if '*' in up_limit:
                                        up_limit = float(line_parts[9][:-2])
                                    else:
                                        up_limit = float(line_parts[9][:-1])
                                    drug_standards['CarbonPiSASA'] = (low_limit, up_limit)
                            elif 'Weakly' in line_parts[1]:
                                weakly_polar_sasa = float(line_parts[5])
                                properties[drug_id]['WeakPolarSASA'] = weakly_polar_sasa
                                if 'WeakPolarSASA' not in drug_standards.keys():
                                    low_limit = float(line_parts[7])
                                    up_limit = line_parts[9]
                                    if '*' in up_limit:
                                        up_limit = float(line_parts[9][:-2])
                                    else:
                                        up_limit = float(line_parts[9][:-1])
                                    drug_standards['WeakPolarSASA'] = (low_limit, up_limit)
                            else:
                                continue
                        elif 'Donor' in line_parts[2]:
                            donors = float(line_parts[7])
                            properties[drug_id]['Donors'] = donors
                            if 'Donors' not in drug_standards.keys():
                                low_limit = float(line_parts[9])
                                up_limit = float(line_parts[11][:-1])
                                drug_standards['Donors'] = (low_limit, up_limit)
                        elif 'Acceptor' in line_parts[2]:
                            acceptors = float(line_parts[7])
                            properties[drug_id]['Acceptors'] = acceptors
                            if 'Acceptors' not in drug_standards.keys():
                                low_limit = float(line_parts[9])
                                up_limit = float(line_parts[11][:-1])
                                drug_standards['Acceptors'] = (low_limit, up_limit)
                        else:
                            continue
                             
                    elif 'Lipinski' in line_parts[0]:
                        violations = line_parts[6]
                        if 'M' in violations:
                            violations = float(line_parts[6][:-1])
                            properties[drug_id]['LipinskiViolations'] = violations
                            if 'MWOutsideTrainingRange' not in properties.keys():
                                properties[drug_id]['MWOutsideTrainRange'] = True
                        else:
                            violations = float(line_parts[6])
                            properties[drug_id]['LipinskiViolations'] = violations
                            if 'MWOutsideTrainingRange' not in properties.keys():
                                properties[drug_id]['MWOutsideTrainingRange'] = False
                        if 'LipinskiViolations' not in drug_standards.keys():
                            drug_standards['LipinskiViolations'] = (0, 5)
                    elif '%' in line_parts[0]:
                        absorption_percent = line_parts[8]
                        if 'M' in absorption_percent:
                            absorption_percent = float(line_parts[8][:-1])
                            properties[drug_id]['PercentAbsorbGI'] = absorption_percent
                            properties[drug_id]['MWOutsideTrainingRange'] = True
                        else:
                            absorption_percent = float(line_parts[8])
                            properties[drug_id]['PercentAbsorbGI'] = absorption_percent
                        if 'PercentAbsorbGI' not in drug_standards.keys():
                            drug_standards['PercentAbsorbGI'] = (25, 100)
                    elif 'Qual.' in line_parts[0]:
                        qual = line_parts[7]
                        if 'M' in qual[-1]:
                            qual = qual[:-1]
                        properties[drug_id]['AbsorbQualModel'] = qual
                    elif 'QP' in line_parts[0]:
                        if 'octanol/water' in line_parts[4]:
                            logp = line_parts[6]
                            if 'M' == logp[0]:
                                logp = float(logp[:-1])
                            else:
                                logp = float(logp)
                            properties[drug_id]['logP'] = logp
                            if 'logP' not in drug_standards.keys():
                                low_limit = float(line_parts[8])
                                up_limit = line_parts[10]
                                if '*' in up_limit:
                                    up_limit = float(up_limit[:-2])
                                else:
                                    up_limit = float(up_limit[:-1])
                                drug_standards['logP'] = (low_limit, up_limit)
                                
                    else:
                        continue

        properties['Standard'] = drug_standards
        # self.properties = properties
        df = pd.DataFrame.from_dict(properties)

        return df

    def qikpropScreen(self, csv=None):
        ignore = ['MWOutsideTrainingRange', 'AbsorbQualModel']
        if csv is not None:
            df = self.getQikpropData(csv=csv)
        else:
            try:
                df = pd.DataFrame(self.properties)
            except:
                print('Please input an .out file or run getQikpropData first.')
                return 0
        data = {}
        for ligand in self.ligands:
            name = ligand.name 
            data[name] = {}
            for prop in df.index:
                if prop in ignore:
                    continue
                data[name][prop] = 0
                standard = df.loc[prop, 'Standard']
                value = df.loc[prop, name]
                if isinstance(standard, tuple):
                    low_lim = standard[0]
                    up_lim = standard[1]
                    if low_lim < value:
                        if value < up_lim:
                            data[name][prop] = 1
                elif prop == 'LipinskiViolations':
                    max_ = 4
                    if value <= max_:
                        data[name][prop] = 1
                elif prop == 'PercentAbsorbGI':
                    min_ = 25
                    if value > min_:
                        data[name][prop] = 1
                else:
                    continue
        df = pd.DataFrame(data)
        df['sum'] = df.sum(axis=1)
        column_sum = df.sum(axis=0)
        column_sum.name = 'sum'
        df = df.append(column_sum)
        return df

                


    def fingerprint(self, csv, output_residues=None, output_contact_types=None):
        '''
        Gets fingerprint data from Schrodinger .csv output. 
        Arguments: 
        ** csv: string. Path to .csv file.
        ** output_residues (optional): list. List of residue IDs to output residue-specific information. 
        '''
        # initialize
        self.contact_types = ['contact', 'backbone', 'sidechain', 'polar', 'hydrophobic', 'acceptor', 'donor', 'aromatic', 'charged']
        self.fingerprints = {} # holds contact data
        self.fingerprints['contact'] = {}
        self.fingerprints['backbone'] = {}
        self.fingerprints['sidechain'] = {}
        self.fingerprints['polar'] = {}
        self.fingerprints['hydrophobic'] = {}
        self.fingerprints['acceptor'] = {}
        self.fingerprints['donor'] = {}
        self.fingerprints['aromatic'] = {}
        self.fingerprints['charged'] = {}
        contact_indeces = {} # holds index position of contact type for residue in fingerprint csv
        backbone_indeces = {}
        sidechain_indeces = {}
        polar_indeces = {}
        hydrophobic_indeces = {}
        acceptor_indeces = {}
        donor_indeces = {}
        aromatic_indeces = {}
        charged_indeces = {}
        for residue in self.residues:
            # get residue number
            residue_number = int(residue[3:])
            # get index position in csv for each interaction type
            contact_index = 9*(residue_number) - 8
            backbone_index = 9*(residue_number) - 7
            sidechain_index = 9*(residue_number) - 6
            polar_index = 9*(residue_number) - 5
            hydrophobic_index = 9*(residue_number) - 4
            acceptor_index = 9*(residue_number) - 3
            donor_index = 9*(residue_number) - 2
            aromatic_index = 9*(residue_number) - 1
            charged_index = 9*(residue_number)
            # add to dictionary
            contact_indeces[residue] = contact_index
            backbone_indeces[residue] = backbone_index
            sidechain_indeces[residue] = sidechain_index
            polar_indeces[residue] = polar_index
            hydrophobic_indeces[residue] = hydrophobic_index
            acceptor_indeces[residue] = acceptor_index
            donor_indeces[residue] = donor_index
            aromatic_indeces[residue] = aromatic_index
            charged_indeces[residue] = charged_index
        # parse csv
        f = open(csv, 'r')
        i = 0
        for line in f:
            if i == 0: # skip first line
                i += 1
                continue
            else:
                line_parts = line.split(',')
                ligand = line_parts[0]
                print(ligand)
                if ligand not in self.names:
                    continue
                else:
                    self.fingerprints['contact'][ligand] = {}
                    for key in contact_indeces.keys():
                        index = contact_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['contact'][ligand][key] = num_contacts

                    self.fingerprints['backbone'][ligand] = {}
                    for key in backbone_indeces.keys():
                        index = backbone_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['backbone'][ligand][key] = num_contacts

                    self.fingerprints['sidechain'][ligand] = {}
                    for key in sidechain_indeces.keys():
                        index = sidechain_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['sidechain'][ligand][key] = num_contacts

                    self.fingerprints['polar'][ligand] = {}
                    for key in polar_indeces.keys():
                        index = polar_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['polar'][ligand][key] = num_contacts

                    self.fingerprints['hydrophobic'][ligand] = {}
                    for key in hydrophobic_indeces.keys():
                        index = hydrophobic_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['hydrophobic'][ligand][key] = num_contacts

                    self.fingerprints['acceptor'][ligand] = {}
                    for key in acceptor_indeces.keys():
                        index = acceptor_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['acceptor'][ligand][key] = num_contacts

                    self.fingerprints['donor'][ligand] = {}
                    for key in donor_indeces.keys():
                        index = donor_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['donor'][ligand][key] = num_contacts

                    self.fingerprints['aromatic'][ligand] = {}
                    for key in aromatic_indeces.keys():
                        index = aromatic_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['aromatic'][ligand][key] = num_contacts

                    self.fingerprints['charged'][ligand] = {}
                    for key in charged_indeces.keys():
                        index = charged_indeces[key]
                        num_contacts = int(line_parts[index])
                        self.fingerprints['charged'][ligand][key] = num_contacts
        f.close()
        # output dataframe fingerprint data about specific residues
        if output_residues is not None:
            output_residue_data = {}
            for residue in output_residues:
                temp = {}
                for contact_type in self.fingerprints.keys():
                    temp[contact_type] = {}
                    for ligand in self.fingerprints[contact_type].keys():
                        contact = self.fingerprints[contact_type][ligand][residue]
                        temp[contact_type][ligand] = contact
                df = pd.DataFrame(temp)
                output_residue_data[residue] = df
            return output_residue_data
        elif output_contact_types is not None:
            if isinstance(output_contact_types, list):
                output_contact_data = {}
                for output_contact_type in output_contact_types:
                    df = pd.DataFrame(self.fingerprints[output_contact_type])
                    output_contact_data[output_contact_type] = df
                return output_contact_data
            elif isinstance(output_contact_types, str):
                if output_contact_types != 'all':
                    output_contact_data = self.fingerprints[output_contact_types]
                    df = pd.DataFrame(output_contact_data)
                    return df
                else:
                    output_contact_data = {}
                    for contact_type in self.contact_types:
                        df = pd.DataFrame(self.fingerprints[contact_type])
                        output_contact_data[contact_type] = df
                    return output_contact_data
        else:
            return self.fingerprints
    ######################################################################################
    ############################ Stats & Other Calculations ##############################
    ######################################################################################
    
    def rmsd(self, ligand, reference, remap=False):
        '''
        Calculates RMSD between ligand and reference. 
        Arguments:
        ** ligand: Ligand object. Specifies which ligand to calculate RMSD. 
        ** reference: string or Ligand object. Specifies which reference to use. Can pass reference structure file, or an instance of Ligand class. 
        ** remap (optional): boolean. Specifies if ligand must be remapped to reference. Default is false. 
        '''
        if isinstance(reference, str):
            reference = Ligand(structure=reference)
        if remap == True:
            ligand.coordinates = utilities.Remap(ligand=ligand, reference=reference).remapped
        total_squared_distance = 0
        for i in range(len(ligand.coordinates)):
                ligand_coord = ligand.coordinates[i]
                ref_coord = reference.coordinates[i]
                squared_distance = ((ligand_coord[0] - ref_coord[0])**2 + (ligand_coord[1] - ref_coord[1])**2 + (ligand_coord[2] - ref_coord[2])**2)
                total_squared_distance += squared_distance
        rmsd = (total_squared_distance / len(ligand.coordinates)) ** 0.5
        return rmsd

#### personal implementation functions
def qikprop_graph():
    alts=['DB01988', 'DB07456', 'DB08450', 'DB11791', 'DB12983','DB13014', 'DB14773', 'DB15637']
    residues=['ASN142', 'CYS145', 'HIS164', 'PHE140', 'GLU166', 'MET165', 'HIS41', 'HIS163', 'LEU141', 'GLY143', 'LEU27', 'MET49', 'ARG188', 'GLN189', 'THR25', 'SER144', 'TYR54', 'ASP187']
    # unfunc = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../TopDrugs/6lu7/1069/ligand_out.pdbqt','../TopDrugs/6lu7/1516/ligand_out.pdbqt','../TopDrugs/6lu7/4877/ligand_out.pdbqt','../TopDrugs/6lu7/5822/ligand_out.pdbqt','../TopDrugs/6lu7/6908/ligand_out.pdbqt','../TopDrugs/6lu7/7585/ligand_out.pdbqt','../TopDrugs/6lu7/7791/ligand_out.pdbqt','../TopDrugs/6lu7/7813/ligand_out.pdbqt','../TopDrugs/6lu7/8778/ligand_out.pdbqt','../TopDrugs/6lu7/8887/ligand_out.pdbqt','../TopDrugs/6lu7/9201/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
    # rank0 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1069/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_1516/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7585/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8887/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/0/0/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
    # rank1 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1069/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_1516/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7585/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8887/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/1/0/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
    # rank2 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1069/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_1516/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7585/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8887/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/2/0/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
    # rank3 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1069/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_1516/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7585/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8887/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/3/0/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
    # rank4 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1069/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_1516/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7585/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8887/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/4/0/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
    # systems = [unfunc, rank0, rank1, rank2, rank3, rank4]
    unfunc = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../TopDrugs/6lu7/1516/ligand_out.pdbqt','../TopDrugs/6lu7/4877/ligand_out.pdbqt','../TopDrugs/6lu7/5822/ligand_out.pdbqt','../TopDrugs/6lu7/6908/ligand_out.pdbqt','../TopDrugs/6lu7/7791/ligand_out.pdbqt','../TopDrugs/6lu7/7813/ligand_out.pdbqt','../TopDrugs/6lu7/8778/ligand_out.pdbqt','../TopDrugs/6lu7/9201/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778', '9201'], alt=alts)
    rank0 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1516/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/0/0/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778','9201'], alt=alts)
    systems = [unfunc, rank0]
    qikprop_all = []
    i = -1
    names=['1516', '4877', '5822', '6908','7791', '7813', '8778','9201']
    for system in systems:
        if i == -1:
            qikprop_path = '/Users/kelsieking/Documents/Schrodinger/unfunc_qikprop/unfunc_qikprop.out'
            qikprop = system.getQikpropData(csv=qikprop_path)
        else:
            prefix = 'func_' + str(i)
            qikprop_path = '/Users/kelsieking/Documents/Schrodinger/' + prefix + '_qikprop/' + prefix + '_qikprop.out'
            qikprop = system.getQikpropData(csv=qikprop_path)
        i += 1
        qikprop_all.append(qikprop)

    grapher = Grapher(screener=rank0)
    angstrom_symbol = r'$\AA$'
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_mw.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='MW', xlabel='Compound', ylabel='Molecular Weight (g/mol)')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_mv.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='MV', xlabel='Compound', ylabel= 'Molecular Volume (' + angstrom_symbol + '$^{3}$' + ')')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_totalsasa.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='TotalSASA', xlabel='Compound', ylabel= 'Total SASA (' + angstrom_symbol + '$^{2}$' + ')')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_hydrophobicsasa.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='HydrophobicSASA', xlabel='Compound', ylabel= 'Hydrophobic SASA (' + angstrom_symbol + '$^{2}$' + ')')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_hydrophilicsasa.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='HydrophilicSASA', xlabel='Compound', ylabel= 'Hydrophilic SASA (' + angstrom_symbol + '$^{2}$' + ')')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_donors.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='Donors', xlabel='Compound', ylabel='Number of Hydrogen Bond Donors')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_acceptors.png', alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='Acceptors', xlabel='Compound', ylabel='Number of Hydrogen Bond Acceptors')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_ev.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='eV', xlabel='Compound', ylabel='Electron Affinity (kJ/mol)')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_logp.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='logP', xlabel='Compound', ylabel='logP (octonol/water)')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_lipinski.png',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='LipinskiViolations', xlabel='Compound', ylabel='Number of Lipinski Violations')
    grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_percentabsgi',alt=True, dfz=qikprop_all[0], dfn=qikprop_all[1], labelz='Bare', labeln='Functionalized', index='PercentAbsorbGI', xlabel='Compound', ylabel='Human Oral Absorption in GI (%)')
    names=['1516', '4877', '5822', '6908','7791', '7813', '8778','9201']
# residues=['ASN142', 'CYS145', 'HIS164', 'PHE140', 'GLU166', 'MET165', 'HIS41', 'HIS163', 'LEU141', 'GLY143', 'LEU27', 'MET49', 'ARG188', 'GLN189', 'THR25', 'SER144', 'TYR54', 'ASP187']
# unfunc = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../TopDrugs/6lu7/1069/ligand_out.pdbqt','../TopDrugs/6lu7/1516/ligand_out.pdbqt','../TopDrugs/6lu7/4877/ligand_out.pdbqt','../TopDrugs/6lu7/5822/ligand_out.pdbqt','../TopDrugs/6lu7/6908/ligand_out.pdbqt','../TopDrugs/6lu7/7585/ligand_out.pdbqt','../TopDrugs/6lu7/7791/ligand_out.pdbqt','../TopDrugs/6lu7/7813/ligand_out.pdbqt','../TopDrugs/6lu7/8778/ligand_out.pdbqt','../TopDrugs/6lu7/8887/ligand_out.pdbqt','../TopDrugs/6lu7/9201/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
# print(unfunc.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_qikprop/unfunc_qikprop.out'))
# rank0 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1069/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_1516/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7585/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8887/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/0/0/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])
# df = unfunc.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_mmgbsa/unfunc_mmgbsa-out.csv')
# df = rank0.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_0_mmgbsa/func_0_mmgbsa-out.csv')
# dic = unfunc.fingerprint(csv='../fingerprint/unfunc_finger.csv')
# df = pd.DataFrame(dic['polar'])
# print(df.loc[:, '9201'])
# df = pd.DataFrame(dic['aromatic'])
# print(df.loc[:, '9201'])
# print('******************')
# print('SB')
# print('******************')
# df = unfunc.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_mmgbsa/unfunc_mmgbsa-out.csv')
# print(df)
# print('******************')
# print('BB')
# print('******************')
# df = unfunc.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_bigbox_mmgbsa/unfunc_bigbox_mmgbsa-out.csv')
# print(df)
# print('******************')
# print('******************')
# print('SB')
# print('******************')
# df = rank0.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_0_mmgbsa/func_0_mmgbsa-out.csv')
# print(df)
# print('******************')
# print('BB')
# print('******************')
# df = rank0.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/functionalized_0_bb_mmgbsa/functionalized_0_bb_mmgbsa-out.csv')
# print(df)
# # print('******************')
# data = rank0.residueInteractions()
# print(data)
# data = rank0.vinaEnergy()
# print(data)
# df = rank0.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_0_mmgbsa/func_0_mmgbsa-out.csv')
# print(df)
# df = rank0.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_0_qikprop/func_0_qikprop.out')
# print(df)
# shit = unfunc.fingerprint(csv='../fingerprint/unfunc_finger.csv', output_contact_types=['contact'])
# shit = shit['contact']
# shit2 = rank0.fingerprint(csv='../fingerprint/func_0_finger.csv', output_contact_types=['contact'])
# shit2 = shit2['contact']
# new = shit2 - shit
# print(shit)
# print(shit2)
# print('\n')
# df = unfunc.minimumDistance()
# df2 = rank0.minimumDistance()
# new = df2 - df
# print(new)
# unfunc = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../TopDrugs/6lu7/1069/ligand_out.pdbqt','../TopDrugs/6lu7/1516/ligand_out.pdbqt','../TopDrugs/6lu7/4877/ligand_out.pdbqt','../TopDrugs/6lu7/5822/ligand_out.pdbqt','../TopDrugs/6lu7/6908/ligand_out.pdbqt','../TopDrugs/6lu7/7585/ligand_out.pdbqt','../TopDrugs/6lu7/7791/ligand_out.pdbqt','../TopDrugs/6lu7/7813/ligand_out.pdbqt','../TopDrugs/6lu7/8778/ligand_out.pdbqt','../TopDrugs/6lu7/8887/ligand_out.pdbqt','../TopDrugs/6lu7/9201/ligand_out.pdbqt'], residues=residues, names=['1069', '1516', '4877', '5822', '6908', '7585', '7791', '7813', '8778', '8887','9201'])

# unfunc = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../TopDrugs/6lu7/1516/ligand_out.pdbqt','../TopDrugs/6lu7/4877/ligand_out.pdbqt','../TopDrugs/6lu7/5822/ligand_out.pdbqt','../TopDrugs/6lu7/6908/ligand_out.pdbqt','../TopDrugs/6lu7/7791/ligand_out.pdbqt','../TopDrugs/6lu7/7813/ligand_out.pdbqt','../TopDrugs/6lu7/8778/ligand_out.pdbqt','../TopDrugs/6lu7/9201/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908','7791', '7813', '8778','9201'])
# rank0 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1516/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/0/0/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778','9201'])
# rank1 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1516/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/1/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/1/0/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778','9201'])
# rank2 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1516/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/2/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/2/0/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778','9201'])
# rank3 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1516/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/3/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/3/0/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778','9201'])
# rank4 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1516/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/4/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/4/0/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778','9201'])
# dfu = unfunc.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_qikprop/unfunc_qikprop.out')
# dfm = unfunc.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_mmgbsa/unfunc_mmgbsa-out.csv')
# dff = unfunc.fingerprint(csv='../fingerprint/unfunc_finger.csv', output_contact_types='all')
# dff0 = rank0.fingerprint(csv='../fingerprint/func_0_finger.csv', output_contact_types='all')
# for contact_type in rank0.contact_types:
#     print(contact_type)
#     tf = dff0[contact_type]
#     uf = dff[contact_type]
#     if contact_type == 'charged':
#         print(tf.loc['GLU166', '4877'])
#         print(uf.loc['GLU166', '4877'])
#         # print(uf)
#     print(tf-uf)
#     print('\n')
# print(rank0.interactions())
# print(rank0.residueInteractions())
# print(rank0.interactions() - unfunc.interactions())
# df0 = rank0.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_0_qikprop/func_0_qikprop.out')
# df0m = rank0.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_0_mmgbsa/func_0_mmgbsa-out.csv')
# print(df0m)
# print(df0m.loc['mmgbsa',:].min(), df0m.loc['mmgbsa',:].max())
# print(df0m.loc['LE_mmgbsa', :].min(), df0m.loc['LE_mmgbsa', :].max())
# for ligand in rank0.ligands:
#     if ligand.name == '9201':
#         atom_count = ligand.atomCount(heavy=True)
#         print(atom_count)
#         df = rank0.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_0_mmgbsa/func_0_mmgbsa-out.csv')
#         mmgbsa = df.loc['mmgbsa', '9201']
#         print(mmgbsa)
#         print(mmgbsa / atom_count)
# print('\n')
# for ligand in unfunc.ligands:
#     if ligand.name == '9201':
#         atom_count = ligand.atomCount(heavy=True)
#         print(atom_count)
#         df = unfunc.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_mmgbsa/unfunc_mmgbsa-out.csv')
#         mmgbsa = df.loc['mmgbsa', '9201']
#         print(mmgbsa)
#         print(mmgbsa / atom_count)
# df1 = rank1.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_1_qikprop/func_1_qikprop.out')
# df1m = rank1.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_1_mmgbsa/func_1_mmgbsa-out.csv')
# df2 = rank1.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_2_qikprop/func_2_qikprop.out')
# df2m = rank2.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_2_mmgbsa/func_2_mmgbsa-out.csv')
# df3 = rank3.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_3_qikprop/func_3_qikprop.out')
# df3m = rank3.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_3_mmgbsa/func_3_mmgbsa-out.csv')
# df4 = rank4.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_4_qikprop/func_4_qikprop.out')
# df4m = rank4.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_4_mmgbsa/func_4_mmgbsa-out.csv')
# x = dfu.iloc[4,:-1]
# y = dfm.iloc[1,:]
# df = pd.concat([x, y], axis=1).T
# x = df0.iloc[4,:-1]
# y = df0m.iloc[1,:]
# dfn = pd.concat([x,y], axis=1).T
# df = df.join(dfn, rsuffix='_0')
# x = df1.iloc[4,:-1]
# y = df1m.iloc[1,:]
# dfn = pd.concat([x,y], axis=1).T
# df = df.join(dfn, rsuffix='_1')
# x = df2.iloc[4,:-1]
# y = df2m.iloc[1,:]
# dfn = pd.concat([x,y], axis=1).T
# df = df.join(dfn, rsuffix='_2')
# x = df3.iloc[4,:-1]
# y = df3m.iloc[1,:]
# dfn = pd.concat([x,y], axis=1).T
# df = df.join(dfn, rsuffix='_3')
# x = df4.iloc[4,:-1]
# y = df4m.iloc[1,:]
# dfn = pd.concat([x,y], axis=1).T
# df = df.join(dfn, rsuffix='_4')
# utilities.linearRegression(df, title='MMGBSA vs VOLUME')

def qikpropEnergyLR(qiloc, eiloc):
    qikprop_dfs = []
    mmgbsa_dfs = []
    df = unfunc.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_qikprop/unfunc_qikprop.out')
    qikprop_dfs.append(df)
    df = unfunc.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_mmgbsa/unfunc_mmgbsa-out.csv')
    mmgbsa_dfs.append(df)
    df = rank0.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_0_qikprop/func_0_qikprop.out')
    qikprop_dfs.append(df)
    df = rank0.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_0_mmgbsa/func_0_mmgbsa-out.csv')
    mmgbsa_dfs.append(df)
    df = rank1.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_1_qikprop/func_1_qikprop.out')
    qikprop_dfs.append(df)
    df = rank1.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_1_mmgbsa/func_1_mmgbsa-out.csv')
    mmgbsa_dfs.append(df)
    df = rank1.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_2_qikprop/func_2_qikprop.out')
    qikprop_dfs.append(df)
    df = rank2.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_2_mmgbsa/func_2_mmgbsa-out.csv')
    mmgbsa_dfs.append(df)
    df = rank3.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_3_qikprop/func_3_qikprop.out')
    qikprop_dfs.append(df)
    df = rank3.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_3_mmgbsa/func_3_mmgbsa-out.csv')
    mmgbsa_dfs.append(df)
    df = rank4.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_4_qikprop/func_4_qikprop.out')
    qikprop_dfs.append(df)
    df = rank4.mmgbsa(csv='/Users/kelsieking/Documents/Schrodinger/func_4_mmgbsa/func_4_mmgbsa-out.csv')
    mmgbsa_dfs.append(df)
    for i in range(len(qikprop_dfs)):
        x = qikprop_dfs[i].iloc[qiloc,:-1]
        y = mmgbsa_dfs[i].iloc[eiloc,:]
        dfn = pd.concat([x, y], axis=1).T
        if i == 0:
            df = dfn
            continue
        else:
            suffix = '_' + str(i - 1)
            df = df.join(dfn, rsuffix=suffix)
    shit = 'shit_' + str(qiloc) + '.csv'
    print(shit)
    df.to_csv(shit)

# for i in range(0,10):
#     qikpropEnergyLR(qiloc=i, eiloc=1)
# modified = System(protein='../TopDrugs/6lu7_protein.pdbqt',ligands=['../scipaper_compare/11a/6lu7/11a.pdbqt', '../scipaper_compare/11a_benzene/pdbqts/11a_benzene_no2_meta/11a_benzene_no2_meta_out.pdbqt', '../scipaper_compare/11a_benzene_no2_glactam/11a_benzene_no2_glactam_out.pdbqt'], names=['11a_out_1', '11a_benzene_no2_meta_out_1', '11a_benzene_no2_glactam_out_1'], residues=residues)
# dfs = modified.fingerprint(csv='../fingerprint/11a_modify_finger.csv', output_contact_types='all')
# print(modified.mmgbsa(csv='../11a_modify_energies.csv'))
# for ctype in modified.contact_types:
#     print(ctype)
#     print(dfs[ctype])
#     print('\n')
    # print('****************')
    # df = dfs[ctype]
    # shits = []
    # names = []
    # for name in modified.names:
    #     names.append(name)
    #     shits.append(df.loc[:, name])
    # print(names[1], names[0])
    # print(shits[1] - shits[0])
    # print('\n')
    # print(names[2], names[0])
    # print(shits[2] - shits[0])
    # print('\n')
    
# sadness = System(protein='../TopDrugs/6lu7_protein.pdbqt',ligands=['../Top5_for_Top_11_1/Top5_9201/0/0/ligand_out.pdbqt', '../Top5_for_Top_11_1/Top5_9201/1/0/ligand_out.pdbqt', '../Top5_for_Top_11_1/Top5_9201/2/0/ligand_out.pdbqt', '../Top5_for_Top_11_1/Top5_9201/3/0/ligand_out.pdbqt'], names=['9201_0', '9201_1', '9201_2', '9201_3'], residues=residues)
# df = sadness.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/9201_all_qikprop/9201_all_qikprop.out')
# print(df)
# alts=['DB01988', 'DB07456', 'DB08450', 'DB11791', 'DB12983','DB13014', 'DB14773', 'DB15637']

# unfunc = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../TopDrugs/6lu7/1516/ligand_out.pdbqt','../TopDrugs/6lu7/4877/ligand_out.pdbqt','../TopDrugs/6lu7/5822/ligand_out.pdbqt','../TopDrugs/6lu7/6908/ligand_out.pdbqt','../TopDrugs/6lu7/7791/ligand_out.pdbqt','../TopDrugs/6lu7/7813/ligand_out.pdbqt','../TopDrugs/6lu7/8778/ligand_out.pdbqt','../TopDrugs/6lu7/9201/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778', '9201'], alt=alts)
# rank0 = System(protein='../TopDrugs/6lu7_protein.pdbqt', ligands=['../Top5_For_Top_11_All/Top5_1516/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_4877/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_5822/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_6908/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7791/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_7813/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_8778/0/0/ligand_out.pdbqt','../Top5_For_Top_11_All/Top5_9201/0/0/ligand_out.pdbqt'], residues=residues, names=['1516', '4877', '5822', '6908', '7791', '7813', '8778','9201'], alt=alts)
# df1 = unfunc.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/unfunc_qikprop/unfunc_qikprop.out')
# df2 = rank0.getQikpropData(csv='/Users/kelsieking/Documents/Schrodinger/func_0_qikprop/func_0_qikprop_changed2fit9201_allout.out')   
# grapher = Grapher(screener=rank0)        
# grapher.doubleBar(saveto='../figures/prop_compare/prop_compare_percentabsgi_test',alt=True, dfz=df1, dfn=df2, labelz='Bare', labeln='Functionalized', index='PercentAbsorbGI', xlabel='Compound', ylabel='Human Oral Absorption in GI (%)')
