# standard python imports
import numpy as np
from numpy.core.shape_base import block
import pandas as pd
from pathlib import Path
import os
import matplotlib
import matplotlib.pyplot as plt
import re
import sys

from pymd.structure.protein import Protein
from pymd.utilities.rewritepdb import writePDB
from pymd.utilities.rewritepdb import editChainIDResidue
# from pymd.mdanalysis.run import Run
# from pymd.mdanalysis.plot_ss_xpm import parseXPM

# if sys.platform == 'linux':
#     matplotlib.use('Agg')


class PostProcess:

    def __init__(self, system=None, tu='ns'):
        if system is not None:
            self.system = system
            # self.directory = self.system.directory
            # self.protein = self.system.protein
        else:
            self.system = None
            self.directory = None
            self.protein = None
            self.peptides = None
        self.tu = tu
        self.data = {}


    @property
    def conversions(self):
        conversions = {
            'ps':0.001,
            'us':1000
        }
        return conversions
    @staticmethod
    def blockAverage(df, block_average):
        d = {}
        for column in df.columns:
            d[column] = {}
            i = 0
            t = 0
            time = 0
            for val in df[column]:
                if i == block_average:
                    t += val
                    avg = t / i
                    t = 0
                    i = 0
                    try:
                        d[column][time/10] = avg
                    except:
                        print(time)
                else:
                    t += val
                i += 1
                time += 1
        df = pd.DataFrame(d)
        return df

    def getDataFrame(self, *args, select=-1):
        df = pd.DataFrame()
        filenum = 1
        convert = False
        name = ''
        for filename in args:
            skiprows = 0
            with open(filename, 'r') as f:
                i = 0
                for line in f:
                    if (i == 0) and (not line.startswith('#')):
                        break
                    elif line.startswith('#'):
                        skiprows += 1
                    elif line.startswith('@'):
                        # check if time is in ps, if yes, convert
                        if (len(line.split()) > 3) and (line.split()[3] == '"Time'):
                            if (line.split()[4].strip('"').strip('(').strip(')') != self.tu):
                                print(line.split()[4])
                                convert = True
                        skiprows += 1
                    else:
                        break
                    i += 1
            if filename.endswith('xvg'):
                data = pd.read_csv(filename, delim_whitespace=True, header=None, skiprows=skiprows, index_col=0)
            else:
                data = pd.read_csv(filename, header=0, index_col=0)
            metadata = self.metadata(filename, df = data).attrs
            if select == -1:
                if metadata['s_list'] != []:
                    data.columns = metadata['s_list']
                for column in data.columns:
                    if column in df.columns:
                        k = 1
                        _col = str(column) + '_{}'.format(str(k))
                        while _col in df.columns:
                            _col = str(column) + '_{}'.format((str(k)))
                            k += 1
                        df[_col] = data[column]
                    else:
                        df[column] = data[column]
            else:
                if metadata['s'] == {}:
                    df[filenum] = data[1]
                else:
                    if len(args) == 1:
                        df = data
                    else:
                        df[filenum] = data[select+1]
            filenum += 1
            name = name + filename + ' '
        df.name = name
        df = self.metadata(*args, df=df)
        self.data = df.attrs
        if convert is True:
            df.index = df.index / 1000
            df.attrs['x_label'] = 'Time (ns)'
        return df

    def remapData(self, df, pdb):
        f = open(pdb, 'r')
        chains = {}
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if 'SOL' in line:
                break
            chainid = line.split()[4]
            resnum = line.split()[5]
            if not chainid in chains.keys():
                chains[chainid] = []
            if resnum not in chains[chainid]:
                chains[chainid].append(resnum)
        f.close()
        for chainid in chains.keys():
            print(chainid, chains[chainid])
            # print(df.loc[chains[chainid],:])

    def timeAverage(self, df, interval=200):
        '''
        parameters:
        df (dataframe)
        interval (int, default=200): interval for averages
        '''
        try:
            assert isinstance(interval, int)
        except:
            raise ValueError('Interval must be integer. Input {} cannot be interpreted as integer'.format(interval))
        end = int(df.index[-1])
        start = int(df.index[0])
        mean = pd.DataFrame()
        sd = pd.DataFrame()
        for i in range(start, end, interval):
            mean['{}-{} {}'.format(i,i+interval, self.tu)] = df.loc[i:i+interval, :].mean()
            sd['{}-{} {}-std'.format(i,i+interval, self.tu)] = df.loc[i:i+interval, :].std()
        data = pd.concat([mean, sd], axis=1)
        return data

    @staticmethod
    def metadata(*args, df):
        title = None
        x_label = None
        y_label = None
        s = {}
        s_list = []
        ptype = None
        _type = None
        file = args[0]
        f = open(file, 'r')
        for line in f:
            if not line.startswith('@'):
                continue
            if 'xaxis' in line.split():
                x_label = line.strip('\n').split('"')[-2]
            if 'yaxis' in line.split():
                y_label = line.strip('\n').split('"')[-2]
            if 'title' in line.split():
                title = line.strip('\n').split('"')[-2]
            if line.split()[1].startswith('s'):
                col = line.split()[1]
                name = line.strip('\n').split('"')[-2]
                s[col] = name
                s_list.append(name)
            if line.startswith('@PTYPE'):
                ptype = line.split()[-1].strip('\n').strip()
            if line.startswith('@TYPE'):
                _type = line.split()[-1].strip('\n').strip()
        if x_label is None:
            x_label = ''
        if y_label is None:
            y_label = ''
        if title is None:
            title = ''
        df.attrs = {
            'title':title,
            'x_label':x_label,
            'y_label':y_label, 
            's':s,
            'ptype':ptype,
            'type':_type,
            's_list':s_list
        }
        f.close()
        return df
        
    def rmsdPerPeptide(self, files, peptides): ## files argument will be taken away -will infer from self.system.directory[rep][rmsd][data]
        
        # find per peptide files
        regex = re.compile('rmsd_\d.xvg')
        per_peptide_files = []
        for f in files:
            filename = f.split(os.sep)[-1]
            if regex.findall(filename) != []:
                per_peptide_files.append(f)
        
        # order them by peptide number
        ordered = []
        for peptide in range(1, peptides+1): #### CHANGE TO SELF.SYSTEM.PEPTIDES AFTER TESTING
            for f in per_peptide_files:
                if 'rmsd_{}.xvg'.format(peptide) in f:
                    ordered.append(f)
                    break
        
        # get data from the xvg
        i = 0
        first = ordered[0]
        convert = False
        f = open(first, 'r')
        for line in f:
            if line.startswith('#'):
                pass
            elif line.startswith('@'):
                try:
                    line_parts = line.split()
                    if 'xaxis' in line_parts[1]:
                        if line_parts[-1] == '(ps)"':
                            convert = True
                except:
                    pass
            else:
                break
            i += 1
        f.close()
        
        # create dataframe
        skiprows = i
        df = pd.read_csv(first, delim_whitespace=True, header=None, index_col=0, skiprows=skiprows, engine='python')
        df.columns = ['peptide1']
        i = 2
        for f in ordered[1:]:
            s = pd.read_csv(f, delim_whitespace=True, header=None, index_col=0, skiprows=skiprows)
            col = 'peptide{}'.format(i)
            s.columns = [col]
            df[col] = s[col]
            i += 1

        return df

    def rmsd(self, group, rep=None, average_replicates=True):
        if rep is None:
            return self.iterate(self.rmsd, group=group, average_replicates=average_replicates)

        filename = os.path.join(self.system.directory[rep]['rmsd']['root'], '{}.xvg'.format(group))
        f = open(filename, 'r')
        time = []
        rmsd = []
        convert = False
        for line in f:
            line_parts = line.split()
            if ('#' in line_parts[0]):
                pass
            elif ('@' in line_parts[0]):
                if line_parts[1] == 'xaxis':
                    if '(ps)"' in line_parts[-1]:
                        convert = True
            else:
                x = float(line_parts[0].strip())
                if convert is True:
                    x = x / 1000
                y = float(line_parts[1].strip())
                time.append(x)
                rmsd.append(y)
        f.close()
        df = pd.DataFrame()
        df['rmsd'] = rmsd
        df.index = time
        return df

    def rmsf(self, group, rep=None, res=False, start=None, stop=None, tu='ns'):
        if rep is None:
            return self.iterate(self.rmsf, group=group, res=res, start=start, stop=stop, tu=tu)
        
        if (start is None) and (stop is None):
            rmsf_filename = os.path.join(self.system.directory[rep]['rmsf']['root'], 'rmsf_{}.xvg'.format(group))
            rmsdev_filename = os.path.join(self.system.directory[rep]['rmsf']['root'], 'rmsdev_{}.xvg'.format(group))
        else:
            if tu == 'ns':
                pass
            if tu == 'ps':
                start = start / 1000
                stop = stop / 1000
            rmsf_filename = os.path.join(self.system.directory[rep]['rmsf']['root'], 'rmsf_{}_{}_{}.xvg'.format(group, start, stop))
            rmsdev_filename = os.path.join(self.system.directory[rep]['rmsf']['root'], 'rmsdev_{}_{}_{}.xvg'.format(group, start, stop))
        
        df = pd.DataFrame()
        # parse rmsf file
        f = open(rmsf_filename, 'r')
        xdata = []
        ydata = []
        for line in f:
            if (line.startswith('#')) or (line.startswith('@')):
                continue
            line_parts = line.split()
            x = int(line_parts[0].strip())
            y = float(line_parts[1].strip())
            xdata.append(x)
            ydata.append(y)
        f.close()

        df['rmsf'] = ydata
        df.index = xdata
        
        # parse rmsdsev file
        f = open(rmsf_filename, 'r')
        ydata = []
        for line in f:
            if (line.startswith('#')) or (line.startswith('@')):
                continue
            line_parts = line.split()
            y = float(line_parts[1].strip())
            ydata.append(y)
        f.close()
        
        df['rmsdev'] = ydata

        return df

    def clusters(self, *args, models=1, sm=None,nterm=None, cterm=None, edit_chainid=True):

        def fixCluster(filepath, models, nterm, cterm):
            contents = []
            model = 1
            f_ = open(filepath, 'r')
            for line in f_:
                if not line.startswith('ENDMDL'):
                    contents.append(line)
                else:
                    # print('writing MODEL {}'.format(model))
                    contents.append(line)
                    new = filepath[:-4] + '_model{}.pdb'.format(model)
                    f = open(new, 'w')
                    for line in contents:
                        f.write(line)
                    f.close()
                    if edit_chainid is True:
                        # if self.system is None:
                        #     print(1)
                        #     protein = Protein(structure=filepath)
                        # else:
                        #     print(2)

                        #     protein = self.system.protein
                        # if nterm is None:
                        #     nterm = protein.nterm
                        # if cterm is None:
                        #     cterm = protein.cterm
                        if self.system is not None:
                            editChainIDResidue(new, new, nterm, cterm, sm=self.system.ligands)
                        else:
                            editChainIDResidue(new, new, nterm, cterm, sm=sm)
                    if model < models:
                        # print('finished model {}'.format(model))
                        model += 1
                        contents = []
                    else:
                        break
            f_.close()


        if (args == tuple()) and (self.system is not None):
            for rep in range(1, self.system.reps+1):
                files = self.system.directory[rep]['clusters']['data']
                regex = re.compile('clusters_\d+_\d+.pdb')
                for f in files:
                    filename = f.split(os.sep)[-1]
                    if (regex.findall(filename) != []) and (filename.count('.') < 2):
                        fixCluster(f, models, nterm,cterm)
        else:
            for f in args:
                fixCluster(f, models,nterm,cterm)

    def clusterPercentages(self, *args, frames=1000):
        df = pd.DataFrame()
        rep = 1
        for filepath in args:
            f = open(filepath, 'r')
            percents = []
            for line in f:
                if (line.startswith('#')) or (line.startswith('@')):
                    continue
                else:
                    percent = ((int(line.split()[-1].strip()) / frames) * 100) 
                    percent = round(percent, 2)
                    percents.append(percent)
            f.close()
            df[rep] = percents[0:6]
            rep += 1
        return df
      
    def clusterPercentagesOLD(self, period=(400,600)):
        df = pd.DataFrame()
        for rep in range(1, self.system.reps+1):
            root = self.system.directory[rep]['clusters']['root']
            if period is not None:
                filepath = os.path.join(root, 'cluster_size_{}_{}.xvg'.format(period[0], period[1]))
            else:
                filepath = os.path.join(root, 'cluster_size.xvg')
            frames = (period[1] - period[0])*100
            f = open(filepath, 'r')
            percents = []
            for line in f:
                if (line.startswith('#')) or (line.startswith('@')):
                    continue
                else:
                    percent = ((int(line.split()[-1].strip()) / frames) * 100) 
                    percent = round(percent, 2)
                    percents.append(percent)
            f.close()
            df[rep] = percents[0:5]
        return df

    def mindist(self, dtype, rep=None, normalize=True, average_peptides=True, testing=False, average_replicates=False):

        if rep is not None:
            pass
        else:
            return self.iterate(self.mindist, dtype=dtype, normalize=normalize, average_peptides=average_peptides, testing=testing, average_replicates=average_replicates)
        
        path = self.system.directory[rep]['mindist'][dtype]['root']
        if normalize is True:
            csv_filename = os.path.join(path, 'distances.csv')
        else:
            csv_filename = os.path.join(path, 'distances_raw.csv')

        # parse ndx
        if normalize is True:
            num_atoms = {}
            current = None
            ndx = os.path.join(self.system.directory[rep]['root'], 'mindist_{}.ndx'.format(dtype))
            f = open(ndx, 'r')
            for line in f:
                if line.startswith('['):
                    current = int(line.strip()[1:-1].strip().split('_')[0][2:])
                    if current not in num_atoms.keys():
                        num_atoms[current] = 0
                else:
                    num_atoms[current] += len(line.split())
            f.close()
        # for key,value in num_atoms.items():
        #     print(key,value)
        # sys.exit(0)

        if (csv_filename not in os.listdir(path)) or (testing == True):    
            print('No distance matrix found in {}'.format(path))   
            print('Creating distance matrix') 
            averages = {}
            print('Parsing files in {}...'.format(path))
            for filename in os.listdir(path):
                if not filename.endswith('xvg'):
                    continue
                filepath = os.path.join(path, filename)
                parts = filename.split('_')
                try:
                    ri_1 = parts[0]
                    ri_2 = parts[1]
                    if 'ri' in ri_1:
                        index_1 = int(ri_1[2:])
                    else:
                        index_1 = int(ri_1[1:])
                    if 'ri' in ri_2:
                        index_2 = int(ri_2[2:])
                    else:
                        index_2 = int(ri_2[1:])
                except:
                    continue
                if index_1 not in averages.keys():
                    averages[index_1] = {}
                if index_2 not in averages.keys():
                    averages[index_2] = {}
                if index_1 != index_2:
                    df = pd.read_csv(filepath, delim_whitespace=True, skiprows=24, header=None)
                    df.columns = ['time', 'ncontacts']
                    avg = df['ncontacts'].mean(axis=0)
                    if normalize is True:
                        num_interactors = num_atoms[index_1] + num_atoms[index_2]
                        avg = avg / num_interactors
                    averages[index_1][index_2] = avg
                    averages[index_2][index_1] = avg
                else:
                    averages[index_1][index_2] = 0
                    averages[index_2][index_1] = 0
            df = pd.DataFrame(averages)
            df.fillna(0, inplace=True)
            df.sort_index(ascending=True, inplace=True)
            df = df.T.sort_index(ascending=True).T
            res_ids = []
            if self.system.peptides is not None:
                for i in range(self.system.peptides):
                    for item in self.protein.ids:
                        res_ids.append(item)
            else:
                res_ids = self.protein.ids

            if dtype == 'sidechain':
                i = 1
                glycine_indeces = []
                missing_indeces = []
                for res_id in res_ids:
                    res_name = res_id[:3]
                    if res_name == 'GLY':
                        glycine_indeces.append(i)
                    if res_name not in list(df.index):
                        missing_indeces.append(i)
                    i += 1
                for missing_index in missing_indeces:
                    if missing_index in glycine_indeces:
                        temp = []
                        for i in range(0, len(df.index)):
                            temp.append(0)
                        df.insert(missing_index-1, missing_index, temp)
                        temp = []
                        for i in range(0, len(df.columns)):
                            temp.append(0)
                        df = df.T
                        df.insert(missing_index-1, missing_index, temp)
            if average_peptides is True:
                cols = []
                rows = []
                # split into columns
                prev = None
                for chain in self.protein.chains.keys():
                    if prev is None:
                        len_chain = len(self.protein.chains[chain]['indeces'])
                        cols.append(df.iloc[:, 0:len_chain])
                    else:
                        len_chain = len_chain + len(self.protein.chains[chain]['indeces'])
                        cols.append(df.iloc[:, prev:len_chain])
                    prev = len_chain
                for d in cols:
                    prev = None
                    for chain in self.protein.chains.keys():
                        if prev is None:
                            len_chain = len(self.protein.chains[chain]['indeces'])
                            rows.append(d.iloc[0:len_chain, :])
                        else:
                            len_chain = len_chain + len(self.protein.chains[chain]['indeces'])
                            rows.append(d.iloc[prev:len_chain, :])
                        prev = len_chain
                reset_index = [r.reset_index(drop=True) for r in rows]
                reset_columns = [r.T.reset_index(drop=True).T for r in reset_index]
                combined = pd.concat(reset_columns)
                df = combined.groupby(level=0).mean()
                df.index = [i for i in range(1, len(df.index)+1)]
                df.columns = [i for i in range(1, len(df.index)+1)]                
            df.to_csv(csv_filename)
        else:
            print('Reading distance matrix {}'.format(csv_filename))
            df = pd.read_csv(csv_filename, index_col=0, header=0)
        return df
    
    def mindistResidueSM(self, group, rep, normalize=True, ndx=None):
        mindist = {}
        path = os.path.join(self.system.directory[rep]['mindist'][group]['root'])
        if group == 'sm_sm':
            groups = ['ringA', 'ringB', 'ringC'] # fixing file naming fuckup for sm_sm
        
        
        for filename in os.listdir(path):
            if filename.endswith('csv'):
                # df = pd.read_csv(os.path.join(path, filename), index_col=0, header=0).T
                # return df
                continue
            if filename.startswith('#'):
                continue
            parts = filename.split('_')
            # get small molecule group and ri
            if group != 'sm_sm':
                if ('distances.csv') not in parts:
                    sm_group = '_'.join(parts[:-2])
                    ri = int(parts[-2][2:])
            # else:
                # if re.match(r'ring.3\d_ring.3\d\.xvg$', filename):
                #     print(filename)
                #     ri = line_parts[0]
                #     sm_group = line_parts[1].split('.')[0]
                # else:
                #     continue
                
            
            '''
            # else: # fixing file naming fuckup for sm_sm
            #     if len(parts) == 3:
            #         sm_group = parts[0]
            #         ri = parts[1]
            #     elif len(parts) == 5:
            #         sm_group = '_'.join(parts[0:2])
            #         ri = '_'.join(parts[2:4])
            #     else:
            #         if parts[0] == self.system.name:
            #             sm_group = parts[0]
            #             ri = '_'.join(parts[1:3])
            #         else:
            #             if (parts[1] == 'hydroxyl') or (parts[1] == 'carbonyl'):
            #                 sm_group = '_'.join(parts[0:2])
            #                 ri = parts[3]
            #             else:
            #                 sm_group = parts[0]
            #                 ri = '_'.join(parts[1:3])
            '''

                     
            if ri not in mindist.keys():
                mindist[ri] = {}
            if sm_group not in mindist[ri].keys():
                mindist[ri][sm_group] = 0
            if group == 'sm_sm':
                if sm_group not in mindist.keys():
                    mindist[sm_group] = {}
                if ri not in mindist[sm_group].keys():
                    mindist[sm_group][ri] = 0
            
            # open xvg & get info
            convert = False
            filepath = os.path.join(path, filename)
            f = open(filepath, 'r')
            lines = 0
            data_lines = 0
            mindist_total = 0
            for line in f:
                line_parts = line.split()
                if (line.startswith('@')) or (line.startswith('#')):
                    # get units, non-data lines
                    lines += 1
                    try:
                        if line_parts[1] == "xaxis":
                            timescale = line_parts[-1][1:3]
                            if timescale != 'ns':
                                convert = True
                    except:
                        continue
                else:
                    # get data lines
                    mindist_total += int(line_parts[1])
                    data_lines += 1
            f.close()

            average = mindist_total / data_lines
            mindist[ri][sm_group] = average
            if group == 'sm':
                mindist[sm_group][ri] = average
            
        df = pd.DataFrame(mindist).T.sort_index(axis=0)

        # fix glycines for sidechain
        if 'sidechain' in group:
            gly_residues = []
            # protein = Protein(structure=gro, ignore=['ACE', 'NH2'], ligands='MYR') # change to self.system.protein
            d = []
            for column in df.columns:
                d.append(0)
            for gly in self.system.protein.glycines:
                df.loc[gly] = d
            df = df.sort_index()
        
        # average peptides
        if (self.system.peptides is not None) and (group != 'sm_sm'):
            dfs = []
            len_peptide = len(self.system.protein.ids) # change to getting from self.system.protein.chains
            num_chains = len(self.system.protein.chains.keys()) # change to getting from self.system.protein.chains
            df = df.reset_index(drop=True)
            for i in range(0, (len_peptide*num_chains), len_peptide):
                dfs.append(df.iloc[i:i+(len_peptide), :])
            reset_index = [df.reset_index(drop=True) for df in dfs]
            combined = pd.concat(reset_index)
            df = combined.groupby(level=0).mean()
            df.index = [i for i in range(1, len(df.index)+1)]

        if normalize is False:
            return df.T        
        
        # normalize
        ndx = os.path.join(self.system.root, str(rep), ndx)
        f = open(ndx, 'r')
        groups = {}
        residue = {}
        current = None
        if self.system.run.types is None:
            self.system.run.getTypes(gro=self.system.gro, set_types=True)
        types = self.system.run.types
        for line in f:
            if line.startswith('['):
                current = line.strip()[1:-1].strip()
                if (current.endswith('Sidechain')) or (current.endswith('Mainchain')):
                    current = current.split('_')[0][2:]
                    residue[current] = 0
                else:
                    if current not in groups.keys():
                        groups[current] = 0
            else:
                test = line.split()[0]
                if not (test in types['protein']):
                    groups[current] += len(line.split())
                else:
                    residue[current] += len(line.split())

        f.close()
        to_del = []
        # for key,value in groups.items():
        #     if value == 0:
        #         to_del.append(key)
        # for key in to_del:
        #     del groups[key]
        for grp,atom_num in groups.items():
            org_values = df[grp]
            new_label = '{}_norm'.format(grp)
            new_values = []
            i = 1
            for value in org_values:
                if i not in self.system.protein.glycines:
                    print(residue[str(i)], atom_num)
                    interactors = residue[str(i)] + atom_num
                    df[new_label] = df[grp] / interactors
                else:
                    df[new_label] = 0
                i += 1
            # if group == 'sm_sm':
            #     df = df.T
            #     df[new_label] = df.T[grp] / atom_num
            #     df = df.T
        df = df.drop(list(groups.keys()), axis=1)
        if group == 'sm_sm':
            df = df.drop(list(groups.keys()), axis=0)
        if group == 'sm_sm':
            df = df.fillna(0)
        output = os.path.join(path, 'distances.csv')
        df.to_csv(output)
        print(df.T)
        return df.T

    def mindistOverTime(self, group, rep=None, block_average=100, normalize=True, average_replicates=False):
        if rep is not None:
            path = self.system.directory[rep]['mindist'][group]['root']
        else:
            return self.iterate(self.mindistOverTime, group=group, average_replicates=average_replicates)

        files = os.listdir(path)
        convert = False
        first = True
        already_done = []
        repeat = False

        # parse ndx
        if normalize is True:
            num_atoms = {}
            current = None
            if group == 'peptide':
                dt = 'peptide_peptide'
            else:
                dt = group
            ndx = os.path.join(self.system.directory[rep]['root'], 'mindist_{}.ndx'.format(dt))
            f = open(ndx, 'r')
            for line in f:
                if line.startswith('['):
                    current = line.strip()[1:-1].strip()
                    # current = current[0] + '_' + current[1]
                    if current not in num_atoms.keys():
                        num_atoms[current] = 0
                else:
                    num_atoms[current] += len(line.split())
            f.close()
        for filename in files:
            split = filename.split('_')
            if group == 'peptide':
                p1 = ''.join(split[:2])
                p2 = ''.join(split[2:4])
            else:
                p1 = split[0]
                p2 = split[1]
            peptides = '{}_{}'.format(p1, p2)
            reverse = '{}_{}'.format(p2, p1)
            if reverse in already_done:
                continue
            data_lines = []
            index = []
            filepath = os.path.join(path, filename)
            f = open(filepath, 'r')
            for line in f:
                line_parts = line.split()
                if (line.startswith('#')) or (line.startswith('@')):
                    try:
                        if 'xaxis' in line_parts[1]:
                            if '(ps)"' in line_parts[-1]:
                                convert = True
                    except:
                        line_parts
                        continue
                else:
                    if first is True:
                        index.append(float(line_parts[0]))
                    if normalize is True:
                        num_interactors = num_atoms[p1] + num_atoms[p2]
                        data_lines.append(float(line_parts[1]) / num_interactors)
                    else:
                        data_lines.append(float(line_parts[1]))
            f.close()
            if first is True:
                df = pd.DataFrame(data=data_lines, index=index, columns=[peptides])
                first = False
                index = list(df.index)
            else:
                df[peptides] = data_lines
            if repeat is True:
                temp = {}
                temp[peptides] = df[peptides]
                temp[reverse] = df[reverse]
                temp_df = pd.DataFrame(temp)
            already_done.append(peptides)
        if convert is True:
            df.index = [i/1000 for i in df.index]
            index = list(df.index)
        if block_average is not None:
            d = {}
            for column in df.columns:
                d[column] = {}
                i = 0
                t = 0
                time = 0
                for val in df[column]:
                    if i == block_average:
                        t += val
                        avg = t / i
                        t = 0
                        i = 0
                        try:
                            d[column][index[time]] = avg
                        except:
                            print(time)
                    else:
                        t += val
                    i += 1
                    time += 1
            df = pd.DataFrame(d)
        return df

    def dsspPerResidue(self, xpm, rep, start=0, stop='frames', peptides=None):
        PATH = 'D:/Work/pymd/'
        # PATH = os.getcwd()
        if peptides is not None:
            chains = self.system.protein.chains
        perl_path = os.path.join(PATH, 'mdanalysis', 'plot_ss_xpm.pl')
        xpm = os.path.join(self.system.directory[rep]['dssp']['root'], xpm)
        # print('perl {} {}'.format(perl_path, xpm))
        # os.system('perl {} {}'.format(perl_path, xpm))
        parseXPM(xpm, start, stop)
        if self.directory is not None:
            f = open('summary_SS.xvg', 'r')
            contents = f.readlines()
            f.close()
            newfilename = os.path.join(self.directory[rep]['dssp']['root'], 'summary_SS.xvg')
            f = open(newfilename, 'w')
            for line in contents:
                f.write(line)
            f.close()
            if os.name == 'nt':
                os.system('del {}'.format('summary_SS.xvg'))
            if os.name == 'posix':
                os.system('rm {}'.format('summary_SS.xvg'))
        else:
            newfilename = os.path.abspath('summary_SS.xvg')
        summary = {
            'alpha':{},
            'beta':{},
            'coil':{}
        }
        f = open(newfilename, 'r')
        for line in f:
            if ('@' not in line) and ('#' not in line):
                line_parts = line.strip().split()
                res_num = int(line_parts[0])
                alpha = sum(list(map(float, line_parts[1:4])))
                beta = sum(list(map(float, line_parts[4:8])))
                coil = float(line_parts[8])
                summary['alpha'][res_num] = alpha
                summary['beta'][res_num] = beta
                summary['coil'][res_num] = coil
        f.close()
        if peptides is not None:
            num_residues = 0
            num_chains = 0
            for key in chains.keys():
                num_residues += len(chains[key]['indeces'])
                num_chains += 1
        df = pd.DataFrame(summary)
        if peptides is not None:
            if len(df.index) % num_residues != 0:
                len_chain = 0
                i = 1
                for key in chains.keys():
                    if i != num_chains:
                        len_chain += len(chains[key]['indeces'])
                        drop = len_chain + i
                        df.drop(labels=drop, axis=0, inplace=True)
                    i += 1
            df.index = [i for i in range(1, num_residues + 1)]
            dfs = []
            len_chain = 0
            prev = None
            for key in chains.keys():
                if len_chain == 0:
                    len_chain = len(chains[key]['indeces'])
                    dfs.append(df.iloc[:len_chain, :])
                    prev = len_chain
                else:
                    len_chain = len_chain + len(chains[key]['indeces'])
                    print(prev, len_chain)
                    dfs.append(df.iloc[prev:len_chain, :])
                    prev += len(chains[key]['indeces'])
            reset_index = [df.reset_index(drop=True) for df in dfs]
            combined = pd.concat(reset_index)
            df = combined.groupby(level=0).mean()
            df.index = [i for i in range(1, len(df.index)+1)]
        if self.directory is None:
            if os.name == 'nt':
                os.system('del {}'.format(newfilename))
            if os.name == 'posix':
                os.system('rm {}'.format(newfilename))
        # p_coil = (df['coil'].sum() / (df['coil'].sum() + df['alpha'].sum() + df['beta'].sum()))*100
        # p_alpha = (df['alpha'].sum() / (df['coil'].sum() + df['alpha'].sum() + df['beta'].sum()))*100
        # p_beta = (df['beta'].sum() / (df['coil'].sum() + df['alpha'].sum() + df['beta'].sum()))*100
        # print('%Coil: {}'.format(p_coil))
        # print('%Alpha: {}'.format(p_alpha))
        # print('%Beta: {}'.format(p_beta))
        # print('\n\n')
        # return df
        # df = pd.DataFrame(summary)
        return df

    def dsspPerResidueAverage(self, gro, xpm, peptides=None):
        dfs = self.iterate(self.dsspPerResidue, gro=gro, xpm=xpm, peptides=peptides)
        alpha = pd.DataFrame()
        beta = pd.DataFrame()
        coil = pd.DataFrame()
        i = 1
        for df in dfs:
            alpha[i] = df['alpha']
            beta[i] = df['beta']
            coil[i] = df['coil']
            if i == 1:
                alpha.index = df.index
                beta.index = df.index
                coil.index = df.index
            i += 1
        dfs = [beta, coil, alpha]
        names = [r'$\beta$-Sheet', 'Coil', r'$\alpha$-Helix']
        i = 0
        for df in dfs:
            df['mean'] = df.mean(axis=1)
            df['stdev'] = df.std(axis=1)
            df.name = names[i]
            i += 1
        return dfs

    def dsspOverTime(self, *args, block_average=2):
        '''
        Post-processing for dssp over time data. Writes CSV.
        Positional Arguments:
        **  any number of .xvgs. if multiple, and if average is true, will return average, if not, will yeild dataframes. if single, will return a single dataframe.
        Returns:
        ** pandas DataFrame
        '''
        # if os.path.isabs(xvg):
        #     f = open(xvg, 'r')
        # elif rep is not None:
        #     xvg = os.path.join(self.system.directory[rep]['dssp']['root'], xvg)
        # else:
        #     return self.iterate(self.dsspOverTime, xvg=xvg, block_average=block_average, average=average, tu=tu, average_replicates=average_replicates)
        dfs = []
        for xvg in args:
            f = open(xvg, 'r')
            i = 0
            header = []
            convert = False
            for line in f:
                line_parts = line.split()
                if line.startswith('#'):
                    pass
                elif (line.startswith('@')):
                    try:
                        if (line_parts[2]) == 'legend':
                            header.append(line_parts[-1][1:-1])
                    except:
                        pass
                    if line_parts[1] == 'xaxis':
                        if '(ps)"' in line_parts[-1]:
                            convert = self.conversions['ps']
                        if 'xm' in line_parts[-1]:
                            convert = self.conversions['us']
                else:
                    break
                i += 1
            f.close()
            df = pd.read_csv(xvg, delim_whitespace=True, header=None, index_col=0, skiprows=i, skipfooter=2, engine='python')
            df.columns = header
            # if convert is not None:
            #     df.index = [i*convert for i in df.index]
            coil_types = ['Coil', 'Bend', 'Turn']
            bsheet_types = ['B-Sheet', 'B-Bridge']
            helix_types = ['3-Helix', 'A-Helix', '5-Helix']
            tmp = pd.DataFrame()
            for coil_type in coil_types:
                if coil_type in df.columns:
                    if 'coil_percent' not in tmp.columns:
                        tmp['coil_percent'] = df[coil_type]
                    else:
                        tmp['coil_percent'] = tmp['coil_percent'] + df[coil_type]
            if 'coil_percent' not in tmp.columns:
                tmp['coil_percent'] = [0]*len(df)
            for bsheet_type in bsheet_types:
                if bsheet_type in df.columns:
                    if 'bsheet_percent' not in tmp.columns:
                        tmp['bsheet_percent'] = df[bsheet_type]
                    else:
                        tmp['bsheet_percent'] = tmp['bsheet_percent'] + df[bsheet_type]
            if 'bsheet_percent' not in tmp.columns:
                tmp['bsheet_percent'] = [0]*len(df)
            for helix_type in helix_types:
                if helix_type in df.columns:
                    if 'helix_percent' not in tmp.columns:
                        tmp['helix_percent'] = df[helix_type]
                    else:
                        tmp['helix_percent'] = tmp['helix_percent'] + df[helix_type]
            if 'helix_percent' not in tmp.columns:
                tmp['helix_percent'] = [0]*len(df)

            tmp['total'] = tmp.sum(axis=1)
            tmp.index = df.index
            df = pd.DataFrame()
            df['coil_percent'] = round((tmp['coil_percent'] / tmp['total']) * 100, 2)
            df['bsheet_percent'] = round((tmp['bsheet_percent'] / tmp['total']) * 100, 2)
            df['helix_percent'] = round((tmp['helix_percent'] / tmp['total']) * 100, 2)
            # if block_average is not None:
            #     df = self.blockAverage(df, block_average)
            if block_average is not None:
                ss = {
                    'coil_percent':{},
                    'bsheet_percent':{},
                    'helix_percent':{}
                }
                i = 0
                cp = 0
                bp = 0
                hp = 0
                for index in df.index:
                    t = index
                    if i != block_average:
                        cp = cp + df.loc[index, 'coil_percent']
                        bp = bp + df.loc[index, 'bsheet_percent']
                        hp = hp + df.loc[index, 'helix_percent']
                    else:
                        cp = cp / (block_average+1)
                        bp = bp / (block_average+1)
                        hp = hp / (block_average+1)
                        ss['coil_percent'][t] = cp
                        ss['bsheet_percent'][t] = bp
                        ss['helix_percent'][t] = hp
                        cp = cp + df.loc[index, 'coil_percent']
                        bp = bp + df.loc[index, 'bsheet_percent']
                        hp = hp + df.loc[index, 'helix_percent']
                        i = 0
                    i += 1    
                df = pd.DataFrame(ss)

            df['total'] = round(df.sum(axis=1),2)


            if len(args) == 1:
                return df
            else:
                dfs.append(df)
        df = pd.concat(dfs).groupby(level=0).mean()
        print(len(dfs))
        new_cols = ['coil_std', 'bsheet_std', 'helix_std']
        cols = ['coil_percent', 'bsheet_percent', 'helix_percent']
        k = 0
        for col in cols:
            stdevs = pd.DataFrame()
            i = 1
            for df in dfs:
                print(i)
                stdevs[i] = df[col]
                i += 1
            df[new_cols[k]] = stdevs.std(axis=1)
            k += 1
        return df

    def hbondsHeatmap(self, dtype, period, rep=None, normalize=False):

        data = {}
        files = self.system.directory[rep]['hbonds'][dtype]['data']
        hbonds_csv = os.path.join(self.system.directory[rep]['hbonds'][dtype]['root'], 'hbonds.csv')
        if hbonds_csv in files:
            df = pd.read_csv(hbonds_csv, index_col=0, header=0)
            return df
        for filepath in files:
            filename = os.path.basename(filepath)
            if filename.endswith('xvg'):
                split = filename.split('_')
                if filename.endswith('hbond.xvg'):
                    if len(split[0]) > 3:
                        s0 = split[0][3:]
                    else:
                        continue
                else:
                    s0 = split[0][:-3]
                if s0 not in data.keys():
                    data[s0] = {}
                if filename.endswith('hbond.xvg'):
                    if len(split[1]) > 3:
                        s1 = split[1][3:]
                    else:
                        continue
                else:
                    s1 = split[-1].split('.')[0][:-3]
                if s1 not in data.keys():
                    data[s1] = {}
                convert = False
                hbond_count = 0
                frames = 0
                f = open(filepath, 'r')
                for line in f:
                    line_parts = line.split()
                    if line.startswith('#'):
                        pass
                    elif (line.startswith('@')):
                        try:
                            if (line_parts[2]) == 'legend':
                                header.append(line_parts[-1][1:-1])
                        except:
                            pass
                        if line_parts[1] == 'xaxis':
                            if '(ps)"' in line_parts[-1]:
                                convert = True
                    else:
                        x = float(line_parts[0].strip())
                        y = int(line_parts[1].strip())
                        if convert is True:
                            x = x / 1000
                        if (x >= period[0]) and (x <= period[1]):
                            hbond_count += y
                            frames += 1
                f.close()
                avg = hbond_count / frames
                data[s0][s1] = avg
                data[s1][s0] = avg
        df = pd.DataFrame(data)
        df.sort_index(axis=1, inplace=True)
        df.fillna(0, inplace=True)
        output = os.path.join(self.system.directory[rep]['hbonds'][dtype]['root'], 'hbonds.csv')
        df.to_csv(output)
        return df

    def hbondsOverTimeSave(self, dtype, rep=None, block_average=100):

        if rep is not None:
            pass
        else:
            return self.iterate(self.hbondsOverTime, dtype=dtype, block_average=block_average)

        xvgs = self.system.directory[rep]['hbonds'][dtype]['data']
        
        dfs = []
        residue_positions = {}
        pos = 0
        _block_average = block_average
        for xvg in xvgs:
            if xvg.endswith('csv'):
                continue
            print(xvg)
            f = open(xvg, 'r')
            ri = int(os.path.basename(xvg).split('_')[0][2:].strip())
            residue_positions[ri] = {}
            xdata = []
            ydata = []
            convert = False
            for line in f:
                line_parts = line.split()
                if ('#' in line_parts[0]):
                    pass
                elif ('@' in line_parts[0]):
                    if line_parts[1] == 'xaxis':
                        if '(ps)"' in line_parts[-1]:
                            convert = True
                else:
                    x = int(line_parts[0])
                    if convert is True:
                        x = x / 1000
                    y = int(line_parts[1])
                    if x < _block_average:
                        xdata.append(x)
                        ydata.append(y)
                    else:
                        xdata.append(x)
                        ydata.append(y)
                        mean = sum(ydata) / len(ydata)
                        residue_positions[ri][_block_average] = mean
                        xdata = []
                        ydata = []
                        _block_average += 100
            f.close()
            _block_average = block_average
        df = pd.DataFrame(residue_positions)
        df = df.reindex(sorted(df.columns), axis=1)
        print(df)
        
    def hbondsMeanSave(self, dtype, rep=None, period=(400,600), normalize=True):
        if rep is None:
            return self.iterate(self.hbondsMean, dtype=dtype, period=period)
        total = 0
        values = []
        frames = 0
        files = self.system.directory[rep]['hbonds'][dtype]['data']
        pairs = {}
        num_files = 0
        for filepath in files:
            filename = os.path.basename(filepath)
            if filename.endswith('xvg'):
                split = filename.split('_')
                if filename.endswith('.xvg'):
                    if len(split) == 3:
                        if (not len(split[0]) > 3) or (not len(split[1]) > 3):
                            if dtype == 'sm':
                                continue
                            else:              
                                s0 = split[0]
                                s1 = split[1]
                        else:
                            s0 = split[0]
                            s1 = split[1]
                    elif len(split) != 3:
                        if (not len(split[0]) > 3) or (not len(split[1]) > 3):
                            if dtype == 'sm':
                                continue
                            else:
                                s0 = split[0]
                                s1 = split[1].split('.')[0]
                        else:
                            s0 = split[0]
                            s1 = split[1].split('.')[0]
                if s0 not in pairs.keys():
                    pairs[s0] = {}
                if s1 not in pairs.keys():
                    pairs[s1] = {}
                if s0 not in pairs[s1].keys():
                    pairs[s1][s0] = 1
                    pairs[s0][s1] = 1
                else:
                    continue
                convert = False
                f = open(filepath, 'r')
                for line in f:
                    line_parts = line.split()
                    if line.startswith('#'):
                        pass
                    elif (line.startswith('@')):
                        try:
                            if (line_parts[2]) == 'legend':
                                header.append(line_parts[-1][1:-1])
                        except:
                            pass
                        if line_parts[1] == 'xaxis':
                            if '(ps)"' in line_parts[-1]:
                                convert = True
                    else:
                        x = float(line_parts[0].strip())
                        if convert is True:
                            x = x / 1000
                        y = int(line_parts[1].strip())
                        if (x >= period[0]) and (x <= period[1]):
                            total += y
                            values.append(y)
                            frames += 1
                f.close()
                num_files += 1
        
        mean = total / frames
        stdev = np.std(values)
        return mean, stdev

    def hbondsTotal(self, rep, dtype, period=(400,600)):
        total = 0
        values = []
        frames = 0
        files = self.system.directory[rep]['hbonds'][dtype]['data']
        csv = os.path.join(self.system.directory[rep]['hbonds'][dtype]['root'], 'total.csv')
        if csv in files:
            f = open(csv, 'r')
            line = f.readlines()[0]
            f.close()
            total = float(line.split(',')[0])
            frames = float(line.split(',')[1])
            return total, frames
        for filepath in files:
            filename = os.path.basename(filepath)
            if filename.endswith('xvg'):
                split = filename.split('_')

                if filename.endswith('hbond.xvg'):
                    if len(split[0]) > 3 and (dtype != 'protein_sm'):
                        s0 = split[0][3:]
                    elif dtype == 'protein_sm':
                        s1 = split[0]
                    else:
                        continue

                if filename.endswith('hbond.xvg'):
                    if (len(split[1]) > 3) and (dtype != 'protein_sm'):
                        s1 = split[1][3:]
                    elif dtype == 'protein_sm':
                        s1 = split[1]
                    else:
                        continue


                convert = False
                f = open(filepath, 'r')
                for line in f:
                    line_parts = line.split()
                    if line.startswith('#'):
                        pass
                    elif (line.startswith('@')):
                        try:
                            if (line_parts[2]) == 'legend':
                                header.append(line_parts[-1][1:-1])
                        except:
                            pass
                        if line_parts[1] == 'xaxis':
                            if '(ps)"' in line_parts[-1]:
                                convert = True
                    else:
                        x = float(line_parts[0].strip())
                        if convert is True:
                            x = x / 1000
                        y = int(line_parts[1].strip())
                        if (x >= period[0]) and (x <= period[1]):
                            total += y
                            values.append(y)
                            frames += 1
                f.close()
        f = open(csv, 'w')
        f.write('{},{}'.format(total, frames))
        f.close()
        return total, frames

    def hbondsOverTime(self, dtype, rep=None, block_average=100):

        if rep is not None:
            pass
        else:
            return self.iterate(self.hbondsOverTimeOLD, dtype=dtype, block_average=block_average)
        by_time = {}
        xvgs = self.system.directory[rep]['hbonds'][dtype]['data']
        # for xvg in xvgs:
        #     if xvg.endswith('csv'):
        #         continue
        #         df = pd.read_csv(xvg, index_col=0, header=0)
        #         time = os.path.basename(xvg).split('_')[-1].split('.')[0]
        #         by_time[time] = df
        # if by_time != {}:
        #     for key, value in by_time.items():
        #         print(key)
        #         print(value)
        #         print('\n')
        #     return by_time                
        

        dfs = []
        residue_positions = {}
        pos = 0
        for xvg in xvgs:
            if xvg.endswith('csv'):
                continue
            print(xvg)
            f = open(xvg, 'r')
            i = 0
            header = []
            convert = False
            for line in f:
                line_parts = line.split()
                if ('#' in line_parts[0]):
                    pass
                elif ('@' in line_parts[0]):
                    if line_parts[1] == 'xaxis':
                        if '(ps)"' in line_parts[-1]:
                            convert = True
                else:
                    break
                i += 1
            f.close()
            header = ['hydrogen bonds', 'pairs']
            df = pd.read_csv(xvg, delim_whitespace=True, header=None, index_col=0, skiprows=i, skipfooter=2, engine='python')
            print(df)
            print('\n')
            df.columns = header
            if convert is True:
                df.index = [i/1000 for i in df.index]
            if block_average is not None:
                last_t_added = None
                hbonds = {}
                i = 0
                hb = 0
                hb_block = []
                for index in df.index:
                    t = index
                    if (t % block_average != 0) or (t == 0):
                        h = df.loc[t, 'hydrogen bonds']
                        hb += h
                        hb_block.append(h)
                    else:
                        hb = hb / i
                        hbonds[t] = {}
                        hbonds[t]['hydrogen bonds'] = hb
                        hbonds[t]['stdev'] = np.std(np.array(hb_block)) 
                        hb = 0
                        hb_block = []
                        i = 0
                        last_t_added = t
                    i += 1  
                if hb_block != []:
                    if round(t) == last_t_added + 100:
                        t = last_t_added + block_average
                    hb = hb / i
                    hbonds[t] = {}
                    hbonds[t]['hydrogen bonds'] = hb
                    hbonds[t]['stdev'] = np.std(np.array(hb_block))
                df = pd.DataFrame(hbonds)
                dfs.append(df)
            else:
                dfs.append(df)

            xvg_base = os.path.basename(xvg)
            ri = xvg_base.split('_')[0]
            residue_positions[ri] = pos
            pos += 1

            blank = pd.DataFrame(0, index=[i for i in df.index], columns=[c for c in df.columns])


        # if self.system.peptides is not None:
        #     if 'sidechain' in dtype:
        #         res_ids = []
        #         for i in range(self.system.peptides):
        #             for item in self.protein.ids:
        #                 res_ids.append(item)
        #         i = 1
        #         glycine_indeces = []
        #         missing_indeces = []
        #         for res_id in res_ids:
        #             res_name = res_id[:3]
        #             if res_name == 'GLY':
        #                 glycine_indeces.append(i)
        #             if res_name not in list(df.index):
        #                 missing_indeces.append(i)
        #             i += 1
        #         print(missing_indeces)
        #         print(glycine_indeces)
        #         for missing_index in missing_indeces:
        #             print(missing_index)
        #             if missing_index in glycine_indeces:
        #                 temp = []
        #                 for i in range(0, len(df.index)):
        #                     temp.append(0)
        #                 df.insert(missing_index-1, missing_index, temp)
        #                 temp = []
        #                 for i in range(0, len(df.columns)):
        #                     temp.append(0)
        #                 df = df.T
        #                 df.insert(missing_index-1, missing_index, temp)


        peptide_dfs = {}
        for chain in self.system.protein.chains.keys():
            indeces = self.system.protein.chains[chain]['indeces']
            i = 1
            for ind in indeces:
                if chain == 'A':
                    peptide_dfs[ind] = []
                ri = 'ri{}'.format(ind)
                try:
                    df_ind = residue_positions[ri]
                    peptide_dfs[i].append(dfs[df_ind])
                except:
                    print('EXCEPTTTT {}'.format(ri))
                    peptide_dfs[i].append(blank)
                i += 1
        for key, value in peptide_dfs.items():
            combined = pd.concat(value)
            avg = combined.groupby(level=0).mean()
            peptide_dfs[key] = avg
        by_residue = peptide_dfs
        # else:
        #     pass
            #TODO: make case when there isn't multiple chains
        by_time = {}
        for res,df in by_residue.items():
            print(res)
            print(df)
            print('\n')
            for column in df.columns:
                col = str(int(column))
                if col not in by_time.keys():
                    by_time[col] = {}
                by_time[col][res] = {}
                by_time[col][res]['hbonds'] = df.loc['hydrogen bonds', column]
                by_time[col][res]['stdev'] = df.loc['stdev', column]
        for time in by_time.keys():
            dic = by_time[time]
            df = pd.DataFrame(dic)
            by_time[time] = df.T
            csv = 'hbonds_{}_{}.csv'.format(int(int(time)-block_average), int(time))
            csv_path = os.path.join(self.system.directory[rep]['hbonds'][dtype]['root'], csv)
            df.T.to_csv(csv_path)            
        return by_time

    def hbondsProteinSM(self, group=None, rep=None, period=(400,600), average_replicates=True):
        if rep is None:
            return self.iterate(self.hbondsProteinSM, group=group, period=period, average_replicates=average_replicates)
        if period is not None:
            if group is None:
                filename = os.path.join(self.system.directory[rep]['hbonds']['protein_sm']['root'], '{}_protein_{}_{}.xvg'.format(self.system.ligands, period[0], period[1]))
            else:
                filename = os.path.join(self.system.directory[rep]['hbonds']['protein_sm']['root'], '{}_protein_{}_{}_{}.xvg'.format(self.system.ligands, group, period[0], period[1]))
        else:
            if group is None:
                filename = os.path.join(self.system.directory[rep]['hbonds']['protein_sm']['root'], '{}_protein.xvg'.format(self.system.ligands))
            else:
                filename = os.path.join(self.system.directory[rep]['hbonds']['protein_sm']['root'], '{}_protein_{}.xvg'.format(self.system.ligands, group))

        convert = False
        time = []
        hbonds = []
        f = open(filename, 'r')
        for line in f:
            line_parts = line.split()
            if ('#' in line_parts[0]):
                pass
            elif ('@' in line_parts[0]):
                if line_parts[1] == 'xaxis':
                    if '(ps)"' in line_parts[-1]:
                        convert = True
            else:
                x = float(line_parts[0].strip())
                if convert is True:
                    x = x / 1000
                y = float(line_parts[1].strip())
                time.append(x)
                hbonds.append(y)
        df = pd.DataFrame()
        df['hbonds'] = hbonds
        df.index = time
        return df

    def gyration(self, group, rep=None, block_average=10, average_replicates=True):
        if rep is None:
            return self.iterate(self.gyration, group=group, block_average=block_average, average_replicates=average_replicates)
        xdata = []
        ydata = []
        convert = False
        gyr = 0
        frames = 0
        if group == 'sm':
            filepath = os.path.join(self.system.directory[rep]['gyration'][self.system.name]['root'], '{}_gyration.xvg'.format(self.system.name))
        else:
            filepath = os.path.join(self.system.directory[rep]['gyration'][group]['root'], '{}_gyration.xvg'.format(group))
        f = open(filepath, 'r')
        for line in f:
            line_parts = line.split()
            if ('#' in line_parts[0]):
                pass
            elif ('@' in line_parts[0]):
                if line_parts[1] == 'xaxis':
                    if '(ps)"' in line_parts[-1]:
                        convert = True
            else:
                x = float(line_parts[0].strip())
                if convert is True:
                    x = x / 1000
                y = float(line_parts[1].strip())
                if x == 0:
                    xdata.append(x)
                    ydata.append(y)
                    frames += 1
                elif frames != block_average:
                    gyr += y
                    frames += 1
                else:
                    gyr += y
                    gyr_avg = gyr / frames
                    xdata.append(x)
                    ydata.append(gyr_avg)
                    gyr = 0
                    frames = 1
        if gyr != 0:
            gyr_avg = gyr / frames
            xdata.append(x)
            ydata.append(gyr_avg)
        
        df = {
            'time':xdata,
            'gyr':ydata
        }
        df = pd.DataFrame(df, index=[i for i in range(0,len(xdata))])
        return df

    def eccentricity(self, group, rep=None, block_average=10):
        if group == 'sm':
            filepath = os.path.join(self.system.directory[rep]['gyration']['ecc']['root'], '{}_ecc.xvg'.format(self.system.name))
        else:
            filepath = os.path.join(self.system.directory[rep]['gyration']['ecc']['root'], '{}_ecc.xvg'.format(group))
        
        xdata = []
        ydata = []
        convert = False
        gyr = 0
        frames = 0
        f = open(filepath, 'r') 
        print(filepath) 
        for line in f:
            line_parts = line.split()
            if ('#' in line_parts[0]):
                pass
            elif ('@' in line_parts[0]):
                if line_parts[1] == 'xaxis':
                    if '(ps)"' in line_parts[-1]:
                        convert = True
            else:
                # print(line)
                x = float(line_parts[0].strip())
                if convert is True:
                    x = x / 1000
                i1, i2, i3 = map(float, line_parts[2:])
                y = self.calcEccentricity(i1,i2,i3)
                if x == 0:
                    xdata.append(x)
                    ydata.append(y)
                    frames += 1
                elif frames != block_average:
                    gyr += y
                    frames += 1
                else:
                    gyr += y
                    gyr_avg = gyr / frames
                    xdata.append(x)
                    ydata.append(gyr_avg)
                    gyr = 0
                    frames = 1
            if gyr != 0:
                gyr_avg = gyr / frames
                xdata.append(x)
                ydata.append(gyr_avg)
        f.close()
        df = {
            'time':xdata,
            'gyr':ydata
        }
        df = pd.DataFrame(df, index=[i for i in range(0,len(xdata))])
        return df

    def iterate(self, function, **kwargs):
        dfs = []
        for i in range(1, self.system.reps+1):
            df = function(**kwargs, rep=i)
            dfs.append(df)
        if 'average_replicates' in kwargs.keys():
            if kwargs['average_replicates'] == True:
                combined = pd.concat(dfs)
                average = combined.groupby(level=0).mean()
                return average
        return dfs

    def contactTypes(self, ndxt):
        '''
        SM-SM contact type sorting
        '''

        contact_types = {
            'ringA':
            {
                'ringA':'hydrophobic',
                'ringB':'hydrophobic',
                'ringC':'hydrophobic',
                'ringA_hydroxyl':None,
                'ringB_hydroxyl':None,
                'ringC_hydroxyl':None,
                'ringC_carbonyl':None,

            },
            'ringB':
            {
                'ringA':'hydrophobic',
                'ringB':'hydrophobic',
                'ringC':'hydrophobic',
                'ringA_hydroxyl':None,
                'ringB_hydroxyl':None,
                'ringC_hydroxyl':None,
                'ringC_carbonyl':None,

            },
            'ringC':
            {
                'ringA':'hydrophobic',
                'ringB':'hydrophobic',
                'ringC':'hydrophobic',
                'ringA_hydroxyl':None,
                'ringB_hydroxyl':None,
                'ringC_hydroxyl':None,
                'ringC_carbonyl':None,

            },
            'ringA_hydroxyl':
            {
                'ringA':None,
                'ringB':None,
                'ringC':None,
                'ringA_hydroxyl':'polar',
                'ringB_hydroxyl':'polar',
                'ringC_hydroxyl':'polar',
                'ringC_carbonyl':'polar',

            },  
            'ringB_hydroxyl':
            {
                'ringA':None,
                'ringB':None,
                'ringC':None,
                'ringA_hydroxyl':'polar',
                'ringB_hydroxyl':'polar',
                'ringC_hydroxyl':'polar',
                'ringC_carbonyl':'polar',

            },  
            'ringC_hydroxyl':
            {
                'ringA':None,
                'ringB':None,
                'ringC':None,
                'ringA_hydroxyl':'polar',
                'ringB_hydroxyl':'polar',
                'ringC_hydroxyl':'polar',
                'ringC_carbonyl':'polar',

            },  
            'ringC_carbonyl':
            {
                'ringA':None,
                'ringB':None,
                'ringC':None,
                'ringA_hydroxyl':'polar',
                'ringB_hydroxyl':'polar',
                'ringC_hydroxyl':'polar',
                'ringC_carbonyl':'polar',

            }
        }

        contacts = {
            'polar':
            {
                'mean':0,
                'stdev':0
            },
            'hydrophobic':{
                'mean':0,
                'stdev':0
            }
        }
        # parse ndxt for atom count data
        counts = {}
        f = open(ndxt, 'r')
        for line in f:
            parts = line.split(':')
            name = parts[0]
            if name == self.system.ligands:
                continue
            counts[name] = len(parts[1].split(',')) * len(self.system.protein.ligand_obj[0].ligands)
        f.close()

        temp = {
            'polar':[],
            'hydrophobic':[]
        }
        groups = ['ringA', 'ringB', 'ringC'] # fixing file naming fuckup for sm_sm
        for rep in range(1, self.system.reps+1):
            path = os.path.join(self.system.directory[rep]['mindist']['sm_sm']['root'])
            for filename in os.listdir(path):
                # parse groups from files
                if not filename.endswith('xvg'):
                    continue
                parts = filename.split('_')
                if len(parts) == 3:
                    group1 = parts[0]
                    group2 = parts[1]
                elif len(parts) == 5:
                    group1 = '_'.join(parts[0:2])
                    group2 = '_'.join(parts[2:4])
                else:
                    if parts[0] == self.system.name:
                        group1 = parts[0]
                        group2 = '_'.join(parts[1:3])
                    else:
                        if (parts[1] == 'hydroxyl') or (parts[1] == 'carbonyl'):
                            group1 = '_'.join(parts[0:2])
                            group2 = parts[3]
                        else:
                            group1 = parts[0]
                            group2 = '_'.join(parts[1:3])

                if (group1 == self.system.ligands) or (group2 == self.system.ligands):
                    continue
                ctype = contact_types[group1][group2]
                if ctype is None:
                    continue

                # parse xvg
                fname = os.path.join(path, filename)
                convert = False
                ydata = []
                f = open(fname, 'r')
                for line in f:
                    line_parts = line.split()
                    if ('#' in line_parts[0]):
                        pass
                    elif ('@' in line_parts[0]):
                        if line_parts[1] == 'xaxis':
                            if '(ps)"' in line_parts[-1]:
                                convert = True
                    else:
                        y = float(line_parts[1].strip())
                        ydata.append(y)
                f.close()

                # append to dict
                average_contacts = sum(ydata) / len(ydata)
                num_atoms = counts[group1] + counts[group2]
                normalized = average_contacts / num_atoms
                temp[ctype].append(normalized)
        
        # average for all reps
        mean_polar = sum(temp['polar']) / len(temp['polar'])
        stdev_polar = np.std(temp['polar'])
        contacts['polar']['mean'] = mean_polar
        contacts['polar']['stdev'] = stdev_polar

        mean_hydrophobic = sum(temp['hydrophobic']) / len(temp['hydrophobic'])
        stdev_hydrophobic = np.std(temp['hydrophobic'])
        contacts['hydrophobic']['mean'] = mean_hydrophobic
        contacts['hydrophobic']['stdev'] = stdev_hydrophobic

        df = pd.DataFrame(contacts)
        return df, temp

    def mindistInteractionRatio(self, group1, group2, period=(400,600)):
        
        df1 = pd.DataFrame()
        df2 = pd.DataFrame()
        for rep in range(1, self.system.reps+1):
            # parse data for groups
            path1 = self.system.directory[rep]['mindist'][group1]['root']
            path2 = self.system.directory[rep]['mindist'][group2]['root']
            paths = [path1, path2]
            for path in paths:
                for filename in os.listdir(path):
                    if path == path1:
                        if group1 != 'sm_sm':
                            fname = os.path.join(path, filename)
                        else:
                            if filename != 'sm_sm_mindist.xvg':
                                continue
                            else:
                                fname = os.path.join(path, filename)
                    else:
                        if group2 != 'sm_sm':
                            fname = os.path.join(path, filename)
                        else:
                            if filename != 'sm_sm_mindist.xvg':
                                continue
                            else:
                                fname = os.path.join(path, filename)
                    xdata = []
                    ydata = []
                    convert = False
                    f = open(fname, 'r')
                    for line in f:
                        line_parts = line.split()
                        if ('#' in line_parts[0]):
                            pass
                        elif ('@' in line_parts[0]):
                            if line_parts[1] == 'xaxis':
                                if '(ps)"' in line_parts[-1]:
                                    convert = True
                        else:
                            x = float(line_parts[0])
                            if convert is True:
                                x = x / 1000
                            y = float(line_parts[-1])
                            if x < period[0]:
                                continue
                            else:
                                xdata.append(x)
                                ydata.append(y)
                    f.close()
                    if path == path1:
                        df1[rep] = ydata
                        if rep == 1:
                            df1.index = xdata
                    if path == path2:
                        df2[rep] = ydata
                        if rep == 1:
                            df2.index = xdata

        df = df1 / df2
        df['mean'] = df.mean(axis=1)
        return df

    def residueContactFrequency(self, group, res1, res2, rep=None, normalize=True, average_peptides=True, average_replicates=True, period=(400,600)):
        df = self.mindist(dtype=group, normalize=normalize, rep=rep, average_peptides=True, average_replicates=True)
        df.index = [i for i in range(1, len(df.index)+1)]
        df.columns = [i for i in range(1, len(df.columns)+1)]
        index1 = self.protein.ids.index(res1) + 1
        index2 = self.protein.ids.index(res2) + 1
        contact_frequency = df.loc[index1, index2]
        return contact_frequency

    def dsspVsResidueContactFrequency(self, group, xvg, res1, res2, period=(400,600)):
        y = []
        x = []
        for rep in range(1, self.system.reps+1):
            bsheet = self.dsspOverTime(xvg=xvg, rep=rep, block_average=1)['bsheet_percent'].loc[period[0]:period[1]].mean()
            y.append(bsheet)
            res = self.residueContactFrequency(group='sidechain', res1=res1, res2=res2, rep=rep, normalize=False, average_replicates=False)
            x.append(res)
        return x, y

    @staticmethod
    def addChargeBFactors(pdb, itp):
        # parse itp
        charges = {}
        i = 0
        read = False
        f = open(itp, 'r')
        for line in f:
            if line.startswith('[ atoms ]'):
                read = True
                continue
            elif (read == False):
                continue
            elif (line.startswith(';   nr')) and (i == 0):
                i += 1
                continue
            elif (line.startswith('[ bonds ]')):
                break
            else:
                line_parts = line.strip().split()
                if line_parts == []:
                    continue
                atom_name = line_parts[4]
                charge = line_parts[6]
                charge_fixed = str(round(float(charge),2)) 
                charges[atom_name] = charge_fixed
        f.close()
        contents = []
        f = open(pdb, 'r')
        for line in f:
            if line.startswith('ATOM'):
                line_parts = line.split()
                atom_name = line_parts[2]
                charge = charges[atom_name]
                line_parts[-2] = charge
                contents.append(line_parts)
            else:
                continue
        f.close()
        writePDB(contents, pdb)

    @staticmethod
    def addAnisotropyBFactors(pdb, nmr):
        f = open(nmr, 'r')
        anis = []
        for line in f:
            if 'Anisotropy' in line:
                line_parts = line.split()
                ani = float(line_parts[-1].strip()) / 100
                ani_str = '{:1.2f}'.format(ani)
                anis.append(ani_str)

        f.close()
        contents = []
        f = open(pdb, 'r')
        i = 0
        for line in f:
            if (line.startswith('ATOM')) or (line.startswith('HETATM')):
                line_parts = line.split()
                line_parts[-2] = anis[i]
                contents.append(line_parts)
                print(line_parts)
                i += 1
            else:
                continue
        f.close()
        writePDB(contents, pdb)

    @staticmethod
    def anisotropyHistogram(nmr):
        f = open(nmr, 'r')
        anis = []
        for line in f:
            if 'Anisotropy' in line:
                line_parts = line.split()
                ani = float(line_parts[-1].strip()) / 100
                ani_str = float('{:1.2f}'.format(ani))
                anis.append(ani_str)
        f.close()
        return anis

    @staticmethod
    def calcEccentricity(i1, i2, i3):
        return np.sqrt(1 - (i1+i2-i3)/(-i1+i2+i3))

    @staticmethod
    def cluster(pdb, nterm, cterm, ligand=None):
        f = open(pdb, 'r')
        contents = f.readlines()
        f.close()
        model1 = []
        i = 0
        for line in contents:
            if 'ENDMDL' not in line:
                line_parts = line.split()
                if 'ATOM' in line_parts[0]:
                    model1.append(line_parts)
            else:
                break
            i += 1
        f.close()
        base = pdb.split('.')[0]
        newfilename = base + '_model1.pdb'
        newpdb = writePDB(model1, newfilename)
        chainid_filename = base + '_chainid.pdb'
        newpdb = editChainIDResidue(filename=newpdb, newfilename=chainid_filename, nterm=nterm, cterm=cterm, sm=ligand)

    @staticmethod
    def getQikpropData(csv):
        all_contents = []
        info = []
        f = open(csv, 'r')
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
            if drug_id == '':
                drug_id = report[2].strip()
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
                        elif ('SASA' in line_parts[2]) or ('SASA' in line_parts[3]):
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
        df = pd.DataFrame.from_dict(properties)
        return df
