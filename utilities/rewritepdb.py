'''
***
Writes a new .pdb file with chain information.
Returns a list containing edited file contents. 
Assumes the .pdb file to be edited does not currently have chain information.
***
* filename: the path to the .pdb file to be edited
* newfilename: the path to the resulting .pdb filename
* cterm: the c-terminal residue number 
***
'''
def addChainID(filename, newfilename, nterm, cterm):
    
    # initialize variables
    contents = [] # hold contents of .pdb file
    new_contents = [] # hold contents of new .pdb file
    next_resnum = None # hold the next residue in the .pdb file
    switch = False # tells when next atom should switch chain id
    chain_id = 'A' # hold current chain ID
    i = 0 # iterator
    
    f = open(filename, 'r') # open .pdb file
    for line in f:
        contents.append(line) # append contents of .pdb file to a list
    f.close()
    
    for line in contents:
        if ('ATOM' in line) and ('SOL' not in line): # check for ATOM line
            line_parts = line.split()
            res_num = line_parts[4] # current residue number
            res_name = line_parts[3]
            # if int(res_num) % cterm == 0: # check if the residue number is the last in the chain
            if res_name == cterm:
                try:
                    parts = contents[i+1].split() # check next atom info
                except:
                    parts = ['TER']
                if 'TER' not in parts: # check that next line is an atom
                    next_name = parts[3]
                    if next_name == nterm:
                        switch = True # increment chain ID if next res name is nterm
                    # format the new line with chain ID, append to newcontents

                if len(line_parts[2]) < 4: 
                    line_format = '{:<4s}{:>7s}  {:<3} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if len(line_parts[2]) == 4:
                    line_format = '{:<4s}{:>7s} {:<4} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if switch is True: # check if the next atom is on the next chain
                    chain_id = chr(ord(chain_id) + 1) # increment chain ID
                    switch = False
            else:
                # format the new line with chain ID, append to new_contents
                if len(line_parts[2]) < 4:
                    line_format = '{:<4s}{:>7s}  {:<3} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if len(line_parts[2]) == 4:
                    line_format = '{:<4s}{:>7s} {:<4} {:>3s} ' + chain_id + '{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
        else:
            # append remark/comment lines to new_contents
            new_contents.append(line)
        i += 1 # incrementor

    # write new .pdb file
    f = open(newfilename, 'w') 
    for item in new_contents:
        f.write(item)
        f.write('          \n') 
    f.close()

    # return list of edited contents
    return newfilename


'''
***
Renumbers several peptides to uniform residue numbering. 
Returns a list containing edited file contents. 
***
* filename: the path to the .pdb file to be edited
* newfilename: the path to the resulting .pdb filename
* cterm: the c-terminal residue number 
* nterm: the n-terminal residue number
***
'''
def renumberPDB(filename, newfilename, nterm, cterm):
    # initialize variables
    contents = [] # hold .pdb file contents
    new_contents = [] # hold edited contents
    num = 1 # hold current residue number
    i = 0 # iterator
    
    # open .pdb, append contents to contents list
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()


    for line in contents:
        if 'ATOM' in line: # check for atom line
            line_parts = line.split() # hold line contents
            res_num = line_parts[5] # hold current residue number
            if int(res_num) % cterm == 0: # check if current residue is the last in the chain
                parts = contents[i+1].split() # check next atom 
                if 'TER' not in parts: # check if next atom is atom
                    next_resnum = int(parts[5]) # hold next residue number
                # format the new line with appropriate residue number, append to newcontents
                if len(line_parts[2]) < 4:
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s}  {:<3} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if len(line_parts[2]) == 4:
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s} {:<4} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                # check if next atom is on a new peptide
                if int(next_resnum) != int(res_num):
                    num = nterm
            else:
                parts = contents[i+1].split() # check next atom
                if 'TER' not in parts: # check if next atom is atom
                    next_resnum = parts[5] # hold next residue number
                # format the new line with appropriate residue number, append to newcontents
                if len(line_parts[2]) < 4: 
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s}  {:<3} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if len(line_parts[2]) == 4:
                    line_parts[5] = str(num)
                    line_format = '{:<4s}{:>7s} {:<4} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}'
                    formatted_list = line_format.format(*line_parts)
                    new_contents.append(formatted_list)
                if int(next_resnum) != int(res_num):
                    num += 1
        else:
            # append comment lines to new_contents
            new_contents.append(line)
        i += 1 # iterator

    # write new .pdb file
    f = open(newfilename, 'w') 
    for line in new_contents:
        f.write(line)
        f.write('      \n')
    f.close()

    # return list of edited contents
    return newfilename

def editChainID(filename, newfilename, nterm, cterm, sm=None):
    first_nterm = False
    chain_id = 'A'
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()
    newcontents = []
    nterm = str(nterm)
    cterm = str(cterm)
    first_sm = False
    sol = False
    i = 0
    for line in contents:
        line_parts = line.split()
        if 'SOL' not in line_parts[3]:
            try:
                atom_num = line_parts[5]
                res_name = line_parts[3]
                # chain_id = line_parts[4]
            except:
                i += 1
                continue
            if atom_num == nterm:
                next_line = contents[i+1]
                next_atom_num = next_line.split()[5]
                if next_atom_num == atom_num:
                    if first_nterm == False:

                        newcontents.append(line)
                    else:
                        line_parts[4] = chain_id
                    atom_name = line_parts[2]
                    if len(atom_name) > 3:
                        string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    else:
                        string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                        newcontents.append(string)
                else:
                    first_nterm = True
                    line_parts[4] = chain_id
                    atom_name = line_parts[2]
                    if len(atom_name) > 3:
                        string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    else:
                        string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    newcontents.append(string)
            elif atom_num == cterm:
                next_line = contents[i+1]
                next_atom_num = next_line.split()[5]
                if next_atom_num == nterm:
                    line_parts[4] = chain_id
                    atom_name = line_parts[2]
                    if len(atom_name) > 3:
                        string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    else:
                        string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    newcontents.append(string)
                    chain_id = chr(ord(chain_id) + 1)
                else:
                    line_parts[4] = chain_id
                    atom_name = line_parts[2]
                    if len(atom_name) > 3:
                        string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    else:
                        string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    newcontents.append(string)
            elif res_name == sm:
                if first_sm == False:
                    chain_id = chr(ord(chain_id) + 1)
                    line_parts[4] = chain_id                
                    atom_name = line_parts[2]
                    if len(atom_name) > 3:
                        string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    else:
                        string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                    newcontents.append(string)
                    first_sm = True
                else:
                    next_line = contents[i+1]
                    try:
                        next_atom_num = next_line.split()[5]
                    except:
                        line_parts[4] = chain_id                
                        atom_name = line_parts[2]
                        if len(atom_name) > 3:
                            string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                        else:
                            string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                        newcontents.append(string)
                    if atom_num != next_atom_num:
                        line_parts[4] = chain_id                
                        atom_name = line_parts[2]
                        if len(atom_name) > 3:
                            string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                        else:
                            string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                        newcontents.append(string)
                        chain_id = chr(ord(chain_id) + 1)
                    else:
                        line_parts[4] = chain_id                
                        atom_name = line_parts[2]
                        if len(atom_name) > 3:
                            string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                        else:
                            string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)

                        newcontents.append(string)
            else:
                line_parts[4] = chain_id
                atom_name = line_parts[2]
                if len(atom_name) > 3:
                    string = '{:<4s}{:>7s} {:^4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)
                else:
                    string = '{:<4s}{:>7s} {:>4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}\n'.format(*line_parts)

                newcontents.append(string)
            i += 1
        else:
            break
    for line in contents[i:]:
        newcontents.append(line)
    f = open(newfilename, 'w')
    for line in newcontents:
        f.write(line)
    f.close()

def formatString(line_parts, chain_id):
    atom_name = line_parts[2]
    line_parts[4] = chain_id
    if len(atom_name) > 3:
        string = '{:<4s}{:>7s} {:<4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}{:>12s}\n'.format(*line_parts)
    else:
        string = '{:<4s}{:>7s}  {:<4s}{:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}  {:>3s}  {:>3s}{:>12s}\n'.format(*line_parts)
    return string

def editChainIDResidue(filename, newfilename, nterm, cterm, sm=None):
    first_nterm = False
    chain_id = 'A'
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()
    newcontents = []
    nterm = str(nterm)
    cterm = str(cterm)
    first_sm = False
    i = 0
    sol = False
    for line in contents:
        line_parts = line.split()
        if 'SOL' not in line:
            try:
                atom_num = line_parts[5]
                res_name = line_parts[3]
                # chain_id = line_parts[4]
            except:
                i += 1
                continue
            if 'ATOM' in line_parts[0]:
                if len(line_parts) < 12:
                    atom_name = line_parts[2]
                    if not (atom_name == 'NA') or (atom_name=='CL'):
                        line_parts.append(atom_name[0])
                    else:
                        add = line_parts[2][0] + line_parts[2][1].lower()
                        line_parts.append(add)
                if res_name == nterm:
                    next_line = contents[i+1]
                    next_res_name = next_line.split()[3]
                    if next_res_name != res_name:
                        first_nterm = True
                    newcontents.append(formatString(line_parts, chain_id))
                elif res_name == cterm:
                    try:
                        next_line = contents[i+1]
                        next_res_name = next_line.split()[3]
                    except:
                        next_res_name = 'TER'
                    if (next_res_name == nterm) or (next_res_name == sm):
                        newcontents.append(formatString(line_parts, chain_id))
                        chain_id = chr(ord(chain_id) + 1)
                    elif next_res_name == 'TER':
                        newcontents.append(formatString(line_parts, chain_id))
                        chain_id = chr(ord(chain_id) + 1)  
                    else:
                        newcontents.append(formatString(line_parts, chain_id))
                elif res_name == sm:
                    res_num = line_parts[5]
                    try:
                        next_line = contents[i+1]
                        next_res_num = next_line.split()[5]
                    except:
                        next_res_num = 'NONE'
                    newcontents.append(formatString(line_parts, chain_id))
                    if next_res_num != res_num:
                        chain_id = chr(ord(chain_id) + 1)
                else:
                    newcontents.append(formatString(line_parts, chain_id))
            else:
                newcontents.append(line)
            i += 1
        else:
            sol = True
            break
    if sol is True:
        for line in contents[i:]:
            newcontents.append(line)
    f = open(newfilename, 'w')
    for line in newcontents:
        f.write(line)
    f.close()

def writePDB(data, newfilename):
    f = open(newfilename, 'w')
    ter = False
    for line in data:
        if ('TER' in line):
            # if ([data][-1] == ['TER']) or ([data][-1] == 'TER'):
            #     ter = True
            #     string = '{}'.format(*line)
            #     f.write(string)
            #     break
            ter = True
            string = '{}'.format(line)
            f.write(string)
            break
        atom_name = line[2]
        atom_symbol = line[2][0]
        res_name = line[3]
        if line[-1] != atom_symbol:
            line.append(atom_symbol)
        if line[4].isalpha():
            if ('SOL' not in line) and ('NA' not in line) and ('CL' not in line):
                if isinstance(line[6], str):
                    xyz = list(map(float, line[6:9]))
                    line[6] = xyz[0]
                    line[7] = xyz[1]
                    line[8] = xyz[2]
            else:
                if isinstance(line[5], str):
                    xyz = list(map(float, line[5:8]))
                    line[5] = xyz[0]
                    line[6] = xyz[1]
                    line[7] = xyz[2]
            if len(atom_name) > 3:
                if len(line[-1]) == 1:
                    string = '{:<4s}{:>7s} {:<4s} {:>3s} {:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>3s}  {:>3s}{:>12s}\n'.format(*line)
                else:
                    string = '{:<4s}{:>7s} {:<4s} {:>3s} {:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>3s}  {:>3s}{:>10s}\n'.format(*line)
            else:
                if len(line[-1]) == 1:
                    string = '{:<4s}{:>7s}  {:<4s}{:>3s} {:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>3s}  {:>3s}{:>12s}\n'.format(*line)
                else:
                    string = '{:<4s}{:>7s}  {:<4s}{:>3s} {:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>3s}  {:>3s}{:>10s}\n'.format(*line)
        else:
            if isinstance(line[5], str):
                xyz = list(map(float, line[5:8]))
                line[5] = xyz[0]
                line[6] = xyz[1]
                line[7] = xyz[2]
            if len(atom_name) > 3:
                string = '{:<4s}{:>7s} {:^4s} {:>3s}  {:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>3s}  {:>3s}{:>12s}\n'.format(*line)
            else:
                string = '{:<4s}{:>7s}  {:<4s}{:>3s}  {:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>3s}  {:>3s}{:>12s}\n'.format(*line)
        f.write(string)
    f.close()
    return newfilename


# editChainIDResidue(filename, newfilename, nterm, cterm)
# filename: the file you want to edit the chain ID for
# newfilename: the output file path
# nterm: the name of the nterm (example: 'ACE')
# cterm: the name of the cterm (example: 'NH2')
# import os
# base = 'C:\\Users\\KelsieKing\\Desktop\\iapp\\systems\\DRUDE\\'
# for filename in os.listdir(base):
#     path = os.path.abspath(os.path.join(base, filename))
#     filebase = filebase.split('.')[0]
#     newfilename = filebase + '_newchain.pdb'
#     newfilepath = os.path.join(base, newfilename)
#     editChainIDResidue(path, newfilepath, 'ACE', 'NH2')
