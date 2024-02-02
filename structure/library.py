def charged():
    return ['ARG', 'LYS', 'ASP', 'GLU']

def pos_charged():
    return ['ARG', 'LYS']

def neg_charged():
    return ['ASP', 'GLU']

def polar():
    return ['SER', 'THR', 'ASN', 'GLN', 'CYS', 'HIS']

def hydrophobic():
    return ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'GLY', 'PRO']

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