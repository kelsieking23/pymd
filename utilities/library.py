
def residues():
    return ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD']

def lipids():
    return ['POPC', 'CHL1', 'SDPE', 'POPE', 'PSM', 'SOPS', 'POPE', 'POPS', 'SM', 'CHOL', 'DLPG', 'DDPC', 'POPI', 'POPG',
            'CER1']

def ions():
    return ['k'.upper(), 'cl'.upper(), 'na'.upper(), 'sod'.upper(), 'cla'.upper()]

def solvent():
        solvent = ['sol', 'tip3p', 'tip3']
        return [sol.upper() for sol in solvent]

def caps():
    caps = ['ace', 'nh2', 'nme']
    return [cap.upper() for cap in caps]

def backbone():
    backbone = ['ca', 'c', 'n']
    return [atom.upper() for atom in backbone]

def mainchain(self):
    mainchain = ['ca', 'c', 'o', 'n', 'hn', 'h', 'ha']
    return [atom.upper() for atom in mainchain]

def gangliosides():
     gangliosides = ['cer', 'bglc', 'bgal', 'ane5']
     return [gang.upper() for gang in gangliosides]

def headgroup_ids():
     headgroup_ids = ['o7', 'p8', 'p9', 'o10', 'o11']
     return [hg.upper() for hg in headgroup_ids]

def dna():
     return ['DT', 'DC', 'DA', 'DG', 'DU']

def lipid():
     return ['POPC', 'CHL1', 'SDPE', 'POPE', 'PSM', 'SOPS', 'POPE', 'POPS', 'SM', 'CHOL', 'DLPG', 'DDPC']

def code_conversions():
    d = {
         'ala':'a',
         'arg':'r',
         'asn':'n',
         'asp':'d',
         'cys':'c',
         'glu':'e',
         'gln':'q',
         'gly':'g',
         'his':'h',
         'hsd':'h',
         'ile':'i',
         'leu':'l',
         'lys':'k',
         'met':'m',
         'phe':'f',
         'pro':'p',
         'ser':'s',
         'thr':'t',
         'trp':'w',
         'tyr':'y',
         'val':'v'
    }
    return d