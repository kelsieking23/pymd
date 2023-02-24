from pymd.utilities.rewritepdb import writePDB, editChainIDResidue
import MDAnalysis as mda

def convertCoordinates(coordinates):
    xyz = list(map(float, coordinates))
    transformed = []
    for coord in xyz:
        c = str(round(coord*10,3))
        while len(c.split('.')[-1]) < 3:
            c = c + '0'
        transformed.append(c)
    return transformed
    
def convertGro(structure):
    universe = mda.Universe(structure)
    newfilename = structure[:-3] + 'pdb'
    with mda.Writer(newfilename) as pdb:
        pdb.write(universe)
    return newfilename
def convertGroOld(structure, ligands=None):
    valid_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD',
                'ACE', 'NH2', 'NME']
    f = open(structure, 'r')
    contents = f.readlines()
    f.close()
    data = []
    i = 0
    k = 0
    residues = []
    sm = None
    atom = 1
    residue = None
    using_solv_res = False
    for line in contents:
        if i < 2:
            i += 1
            continue
        elif line == contents[-1]:
            break
        else:
            line_parts = line.split()
            gro_res_id = line[0:8]
            res_name = ''.join(char for char in gro_res_id if char.isalpha())
            res_num = ''
            chars = 0
            for char in gro_res_id:
                if (char.isnumeric()) or (char == ' '):
                    res_num = res_num + char
                else:
                    break
            res_num = res_num.strip()
            res_name = gro_res_id[chars:]
            res_id = res_name + res_num
            if residue is None:
                print(res_num)
                residue = int(res_num)
            # if (line_parts[0].endswith('NA')) or (line_parts[0].endswith('CL')):
            #     res_name = line_parts[0][-2:]
            # else:
            #     res_name = line_parts[0][-3:]
            # if residue is None:
            #     residue = int(line_parts[0][:-3])
            # res_num = str(residue)
            # res_id = res_name + res_num
            if ligands is not None:
                if res_name in ligands:
                    if k == 0:
                        k = 1
                        sm = res_name
                    else:
                        pass
                else: residues.append(res_name)
            elif (res_name == 'SOL') or (res_name == 'NA') or (res_name == 'CL') or (res_name == 'TIP3'):
                pass
            else:
                residues.append(res_name)
            atom_num = str(atom)
            atom_type = line[8:15].strip()
            x = line[20:28].strip()
            y = line[28:36].strip()
            z = line[36:45].strip()
            x,y,z = convertCoordinates([x,y,z])
            # if len(line_parts) >= 6:
            #     atom_type = line_parts[1]
            #     x, y, z = convertCoordinates(line_parts[3:6])
            # if len(line_parts) == 5:
            #     x, y, z = convertCoordinates(line_parts[2:])
            #     if not (line_parts[1].startswith('HW')):
            #         atom_type = line_parts[1][:2]
            #     else:
            #         atom_type = line_parts[1][:3]
            if res_name not in valid_residues:
                newline = ['ATOM', atom_num, atom_type, res_name, res_num, x, y, z, '1.00', '0.00']
            else:
                newline = ['ATOM', atom_num, atom_type, res_name, 'X', res_num, x, y, z, '1.00', '0.00']
            # if res_name == 'SOL':
            #     newline = ['ATOM', atom_num, atom_type, res_name, res_num, x, y, z, '1.00', '0.00']
            # elif (res_name == 'NA') or (res_name == 'CL'):
            #     newline = ['ATOM', atom_num, atom_type, res_name, res_num, x, y, z, '1.00', '0.00']
            # else:
            #     newline = ['ATOM', atom_num, atom_type, res_name, 'X', res_num, x, y, z, '1.00', '0.00']
            data.append(newline)
            atom += 1
            if atom == 100000:
                atom = 1
            try:
                next_line = contents[i+1]
            except:
                continue
            if next_line != contents[-1]:
                next_line_parts = next_line.split()
                next_res_num = next_line_parts[0][:-3]
                _res_num = line_parts[0][:-3]
                if next_res_num != _res_num:
                    if (next_line_parts[0].endswith('NA')) or (next_line_parts[0].endswith('CL')):
                        next_res_name = next_line_parts[0][-2:]
                    else:
                        next_res_name = next_line_parts[0][-3:]
                    if next_res_name == residues[0]:
                        # print(next_res_name)
                        residue = int(next_res_num)
                    else:
                        residue += 1 
                if residue > 9999:
                    residue = 1
            i += 1
    newfilename = structure[:-3] + 'pdb'
    writePDB(data, newfilename)
    nterm = residues[0]
    cterm = residues[-1]
    editChainIDResidue(newfilename, newfilename, nterm, cterm, sm)
    return newfilename

