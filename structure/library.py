import os
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
from pymd import dir_path


def charged():
    return ['ARG', 'LYS', 'ASP', 'GLU']

def pos_charged():
    return ['ARG', 'LYS']

def neg_charged():
    return ['ASP', 'GLU']

def polar():
    return ['SER', 'THR', 'ASN', 'GLN', 'CYS', 'HIS', 'HSD', 'HSE']

def hydrophobic():
    return ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'GLY', 'PRO']

def _canonical():
    canon = charged() + polar() + hydrophobic()
    unique = []
    for item in canon:
        if item not in unique:
            unique.append(item)
    return unique

def propdict():
    return {
        'ILE':'h',
        'VAL':'h',
        'LEU':'h',
        'PHE':'h,a',
        'CYS':'h',
        'MET':'h',
        'ALA':'h',
        'GLY':'h',
        'THR':'p',
        'SER':'p',
        'TRP':'h,a',
        'TYR':'h,a',
        'PRO':'h',
        'HIS':'p',
        'ASN':'p',
        'GLN':'p',
        'ASP':'c,nc',
        'GLU':'c,nc',
        'LYS':'c,pc',
        'ARG':'c,pc'
    }

def kyte_doolittle_ranks():
    return {
        'ILE':0,
        'VAL':0,
        'LEU':1,
        'PHE':2,
        'CYS':2,
        'MET':2,
        'ALA':2,
        'GLY':3,
        'THR':3,
        'SER':3,
        'TRP':3,
        'TYR':3,
        'PRO':3,
        'HIS':4,
        'ASN':4,
        'GLN':4,
        'ASP':4,
        'GLU':4,
        'LYS':4,
        'ARG':5
    }

def oneletter_to_threeletter():
    return {
        'A':'ALA',
        'R':'ARG',
        'N':'ASN',
        'D':'ASP',
        'C':'CYS',
        'E':'GLU',
        'Q':'GLN',
        'G':'GLY',
        'H':'HIS',
        'I':'ILE',
        'L':'LEU',
        'K':'LYS',
        'M':'MET',
        'F':'PHE',
        'P':'PRO',
        'S':'SER',
        'T':'THR',
        'W':'TRP',
        'Y':'TYR',
        'V':'VAL'
    }

def threeletter_to_oneletter():
    dic = {v:k for k,v in oneletter_to_threeletter().items()}
    dic['HSD'] = 'H'
    return dic


def mjhw_zscores():
    '''
    Scores taken from Table S1:
    https://www.pnas.org/doi/full/10.1073/pnas.2003773117#sec-3
    '''
    oneletters = {
        'A':-0.0645,
        'C':0.502,
        'D':-1.34,
        'E':-1.35,
        'F':1.48,
        'G':-0.382,
        'H':-0.0737,
        'I':1.14,
        'K':-1.47,
        'L':1.33,
        'M':0.767,
        'N':-0.551,
        'P':-0.424,
        'Q':-0.518,
        'R':-1.21,
        'S':-0.540,
        'T':-0.216,
        'V':0.811,
        'W':1.29,
        'Y':0.812
    }
    dic = {}
    for k, v in oneletters.items():
        dic[oneletter_to_threeletter()[k]] = v
    return dic

cannonical = _canonical()


def get_scale(scale, get_as='dict'):
    scale_csv = os.path.join(dir_path, 'structure', 'hydrophobicity_scores.csv')
    df = pd.read_csv(scale_csv)
    df.columns = [col.lower().strip() for col in df.columns]
    if get_as == 'dict':
        _scale = {}
        for code, value in zip(df['3letter'], df[scale]):
            _scale[code] = float(value)
        return _scale
    if get_as == 'list':
        _scale = []
        for code, value in zip(df['3letter'], df[scale]):
            _scale.append([code, float(value)])
        return _scale
