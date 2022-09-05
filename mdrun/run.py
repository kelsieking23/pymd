import os
import sys
import multiprocessing as mp
import pymd
from pymd.mdrun.runfunc import pbc, rmsd
from pymd.mdrun.gmx import GMX
import time
DIRECTORY = pymd.dir_path
# if sys.platform == 'linux':
#     DIRECTORY = '/work/cascades/kelsieking23/iapp_analysis/scripts/python'
class Run:

    def __init__(self, system=None, peptides=None, cascades=True):
        if system is not None:
            self.system = system
            # self.root = system.root
            # self.reps = system.reps
            # self.name = system.name
            # self.source = system.source
            # self.peptides = system.peptides
            # self.ligands = system.ligands
            # self.cascades = system.cascades
            # self.directory = system.directory
            # self.protein = system.protein
            # self.system = system
            # self.email = system.email
            # self.peptides = system.peptides
            # self.alias = system.alias
        else:
            self.system = None
        if cascades is False:
            self.testing = True
        else:
            self.testing = False
        self.types = None
        self.gmx = GMX(self.system, self.testing)

    def cat(self, job_name, walltime, nodes=1, ntasks_per_node=24):
        pass

    def pbc(self, output=None):


        pool = mp.Pool(processes=len(self.system._reps))
        self.system.job_params = {'output':output}
        results = pool.map(pbc, self.system._reps)
        return results

    def rmsd(self, inp, top='system', job_name='rmsd', output='rmsd.xvg', selections=('backbone'), reference_index=0,
        precentered=True, parallel=True):
        '''
        Calculates RMSD. Performs RMSD for all replicats of System in-parallel.
        Parameters:
            inp : str
                Base name of trajectory (.xtc, .trr)
            top : str, optional, default System.gro
                Base name for topology file (.gro, .pdb). Number of atoms in topology must match that of trajectory.
                Note that .tpr is not a valid extension. 
            output : str, optional, default 'rmsd.xvg'
                Output filename. 
            selections : array-like, optional, default ('backbone')
                Selection(s) for RMSD calculation. Can specify multiple selections. Selection algebra is based on MDtraj 
                topology selection. See MDtraj atom selection reference (https://mdtraj.org/1.9.4/atom_selection.html)
            reference_index : int, optional, default 0
        Returns:
            df : pandas.DataFrame
                A dataframe. 
        '''
        
        self.system.job_params = {'inp':inp,'top':top,'output':output, 'selections':selections,
                                'reference_index':reference_index, 'precentered':precentered}
        if parallel:
            start = time.time()
            pool = mp.Pool(processes=len(self.system._reps))
            results = pool.map(rmsd, self.system._reps)
            i = 0
            for rep in self.system._reps:
                rep.rmsd.load(job_name)
                i += 1
            end = time.time()
            # print(end-start)
        else:
            start = time.time()
            results = []
            for rep in self.system._reps:
                rms = rmsd(rep)
                if isinstance(rms, Exception):
                  raise rms
                else:
                    rep.rmsd.load(job_name)
                    results.append(rms)
            end = time.time()
            # print(end-start)
        # self.system.update(results)
        return self.system

    def getHeader(self, job_name, walltime, nodes=1, ntasks_per_node=24):
        header = []
        vrsn = '_'.join(self.source.split('.')) + '.sh'
        base = os.path.join(DIRECTORY, 'mdrun', vrsn)
        f = open(base, 'r')
        for line in f:
            line_parts = line.split()
            try:
                if ('-t' in line_parts[-1]) and ('all' not in line_parts[-1]):
                    line = line.strip() + ' ' + str(walltime) + '\n'
                    header.append(line)
                elif '--job-name=' in line_parts[-1]:
                    line = line.strip() + job_name + '\n'
                    header.append(line)
                elif '--mail-user=' in line_parts[-1]:
                    line = line.strip() + self.email + '\n'
                    header.append(line)
                elif '--nodes=' in line_parts[-1]:
                    line = line.strip() + str(nodes) + '\n'
                    header.append(line)
                elif '--ntasks-per-node=' in line_parts[-1]:
                    if ntasks_per_node is None:
                        continue
                    else:
                        line = line.strip() + str(ntasks_per_node) + '\n'
                        header.append(line)
                else:
                    header.append(line)
            except:
                header.append(line)
        f.close()
        return header
    
    def getSystemXTC(self, run=True):
        '''
        Writes & submits .sh file to create system .xtc (cat_pbc.xtc) from cat.xtc
        '''
        if self.cascades == True:
            filename = os.path.join(self.directory['scripts'], 'cat_pbc_{}.sh'.format(self.name))
        else:
            filename = 'cat_pbc_{}.sh'.format(self.name)     
        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_catpbc'.format(self.name), walltime='15:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in range(1, self.reps+1):
            cmd = 'echo 0 0 | gmx trjconv -f {} -s {} -pbc mol -ur compact -center -o {} \n'.format(self.directory[rep]['xtc_nopbc'], self.directory[rep]['tpr'], os.path.join(self.directory[rep]['root'], 'cat_pbc.xtc'))
            f.write(cmd)
        f.close()
        if run == True:
            tmp = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(tmp)
        return filename

    def getCatPBC(self, run=True):
        '''
        Writes & submits .sh file to create system .xtc (cat_pbc.xtc) from run xtcs
        '''
        if self.cascades == True:
            filename = os.path.join(self.directory['scripts'], 'cat_{}.sh'.format(self.name))
        else:
            filename = 'cat_{}.sh'.format(self.name)
        f = open(filename, 'w')
        header = self.getHeader(job_name='cat_{}'.format(self.name), walltime='15:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in range(1, self.reps+1):
            xtcs = self.directory[rep]['run_xtcs']
            input_string = ''
            for xtc in xtcs:
                if xtc == xtcs[-1]:
                    input_string = input_string + xtc
                else:
                    input_string = input_string + xtc + ' '
            if input_string != '':
                cmd = 'echo 0 0 | gmx trjcat -f {} -o {} \n'.format(input_string, os.path.join(self.directory[rep]['root'], 'cat.xtc'))
                f.write(cmd)
                self.directory[rep]['xtc_nopbc'] = os.path.join(self.directory[rep]['root'], 'cat.xtc')
                f.write('\n')
        f.write('\n\n')
        f.write('wait')
        f.write('\n\n')
        f.close()
        cat_xtc_filename = self.getSystemXTC(run=False)
        f = open(filename, 'a')
        f.write('sbatch {}'.format(cat_xtc_filename))
        f.close()
        if run == True:
            tmp = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(tmp)
        return filename

    def getProteinXTC(self):
        '''
        Writes & submits .sh file to create protein .xtc and dump
        '''
        if self.cascades == True:
            filename = os.path.join(self.directory['scripts'], 'pro_sm_xtc_{}.sh'.format(self.name))
        else:
            filename = 'pro_sm_xtc_{}.sh'.format(self.name)
        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_prosmxtc'.format(self.name), walltime='12:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in range(1, self.reps+1):
            cmd = 'echo 1 1 | gmx trjconv -f {} -s {} -pbc mol -ur compact -center -o {} \n'.format(self.directory[rep]['xtc_system'], self.directory[rep]['tpr'], os.path.join(self.directory[rep]['root'], 'cat_pbc_pro.xtc'))
            f.write(cmd)
        f.close()
        os.system('sbatch {}'.format(filename))

    def getProteinSMXTC(self):
        if self.cascades == True:
            filename = os.path.join(self.directory['scripts'], 'pro_sm_xtc_{}.sh'.format(self.name))
        else:
            filename = 'pro_xtc_{}.sh'.format(self.name)
        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_prosm_xtc'.format(self.name), walltime='12:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in range(1, self.reps+1):
            ndx_filename = os.path.join(self.system.directory[rep]['root'], 'pro_sm.ndx')
            self.makeNDX(rep=rep, gro=self.system.gro, groups=['protein_caps or {}'.format(self.system.ligands)], filename=ndx_filename)
            xtc_system = self.directory[rep]['xtc_system']
            tpr = self.directory[rep]['tpr']
            output = os.path.join(self.directory[rep]['root'], 'cat_pbc_pro_sm.xtc')
            cmd = 'echo 0 0 | gmx trjconv -f {} -s {} -n {} -pbc mol -ur compact -center -o {} \n'.format(xtc_system, tpr, ndx_filename, output)
            f.write(cmd)
            f.write('wait\n\n')
        f.close()
        cwd = os.getcwd()
        os.chdir(self.directory['scripts'])
        os.system('sbatch {}'.format(filename))
        os.chdir(cwd)
        
    def rmsdBackbone(self, run=True):
        '''
        Writes & submits .sh file to perform rmsd calculations
        '''
        if self.cascades == True:
            filename = os.path.join(self.directory['scripts'], 'rmsd_{}.sh'.format(self.name))
        else:
            filename = 'rmsd_{}.sh'.format(self.name)
        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_rmsd'.format(self.name), walltime='10:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in range(1, self.reps+1):
            cmd = 'echo 4 4 | gmx rms -f {} -s {} -o {} -tu ns \n'.format(self.directory[rep]['xtc_system'], self.directory[rep]['tpr'], os.path.join(self.directory[rep]['root'], 'rmsd.xvg'))
            f.write(cmd)
        f.close()
        if run == True:
            tmp = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(tmp)

    def rmsdPerPeptideBackbone(self, rep=None, run=True):
        '''
        Writes & submits .sh file to perform rmsd calculations
        '''
        if self.cascades == True:
            if rep is None:
                filename = os.path.join(self.directory['scripts'], 'rmsd_per_peptide_{}.sh'.format(self.name))
            else:
                filename = os.path.join(self.directory['scripts'], 'rmsd_per_peptide_{}_rep{}.sh'.format(self.name, rep))
        else:
            filename = 'rmsd_per_peptide_{}.sh'.format(self.name)
        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_rmsd_per_peptide'.format(self.name), walltime='10:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        _rep = rep
        for rep in range(1, self.reps+1):
            if _rep is not None:
                if rep != _rep:
                    continue
            if 'rmsd_per_peptide.ndx' not in os.listdir(self.directory[rep]['root']):
                ndx_filename = self.makeNDX(rep, custom='rmsdPerPeptide')
            else:
                ndx_filename = os.path.join(self.directory[rep]['root'], 'rmsd_per_peptide.ndx')

            if self.directory[rep]['xtc_pro_sm'] is not None:
                xtc = self.directory[rep]['xtc_pro_sm']
            elif self.directory[rep]['xtc_pro'] is not None:
                xtc = self.directory[rep]['xtc_pro']
            else:
                xtc = self.directory[rep]['xtc_system']

            for pep in range(0, self.peptides):
                output_path = os.path.join(self.directory[rep]['rmsd']['root'], 'rmsd_{}.xvg'.format(str(pep+1)))
                cmd = 'echo {} {} | gmx rms -f {} -s {} -n {} -o {} -tu ns \n'.format(str(pep), str(pep), xtc, self.directory[rep]['tpr'], ndx_filename, output_path)
                f.write(cmd)
            f.write('\n\n')
        f.close()
        if run == True:
            tmp = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(tmp)
    

    def rmsf(self, group, res=False, start=None, stop=None, run=True):
        '''
        Writes & submits .sh file to perform rmsf calculations
        '''
        filename = 'rmsf_{}_{}.sh'.format(group, self.name)
        filename = os.path.join(self.directory['scripts'], filename)
        f = open(filename,'w')
        header = self.getHeader(job_name='{}_rmsf'.format(self.name), walltime='10:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in range(1, self.reps+1):
            if self.directory[rep]['xtc_pro_sm'] is not None:
                xtc = self.directory[rep]['xtc_pro_sm']
            elif self.directory[rep]['xtc_pro'] is not None:
                xtc = self.directory[rep]['xtc_pro']
            else:
                xtc = self.directory[rep]['xtc_system']
            tpr = self.directory[rep]['tpr']

            ndx_filename = os.path.join(self.directory[rep]['root'], 'rmsf_{}.ndx'.format(group))
            self.makeNDX(rep, groups=(group,), filename=ndx_filename)
            if (start is None) and (stop is None):
                output_path = os.path.join(self.directory[rep]['rmsf']['root'], 'rmsf_{}.xvg'.format(group))
                od_path = os.path.join(self.directory[rep]['rmsf']['root'], 'rmsdev_{}.xvg'.format(group))
                oc_path = os.path.join(self.directory[rep]['rmsf']['root'], 'correl_{}.xvg'.format(group))
            else:
                strt = int(start / 1000)
                stp = int(stop / 1000)
                output_path = os.path.join(self.directory[rep]['rmsf']['root'], 'rmsf_{}_{}_{}.xvg'.format(group, strt, stp))
                od_path = os.path.join(self.directory[rep]['rmsf']['root'], 'rmsdev_{}_{}_{}.xvg'.format(group, strt, stp))
                oc_path = os.path.join(self.directory[rep]['rmsf']['root'], 'correl_{}_{}_{}.xvg'.format(group, strt, stp))
            if res is False:
                if (start is None) and (stop is None):
                    cmd = 'echo 0 | gmx rmsf -f {} -s {} -n {} -o {} -od {} -oc {}\n'.format(xtc, tpr, ndx_filename, output_path, od_path, oc_path)
                else:
                   cmd = 'echo 0 | gmx rmsf -f {} -s {} -n {} -o {} -od {} -oc {} -b {} -e {}\n'.format(xtc, tpr, ndx_filename, output_path, od_path, oc_path, start, stop) 
            else:
                if (start is None) and (stop is None):
                    cmd = 'echo 0 | gmx rmsf -res -f {} -s {} -n {} -o {} -od {} -oc {}\n'.format(xtc, tpr, ndx_filename, output_path, od_path, oc_path)
                else:
                    cmd = 'echo 0 | gmx rmsf -res -f {} -s {} -n {} -o {} -od {} -oc {} -b {} -e {}\n'.format(xtc, tpr, ndx_filename, output_path, od_path, oc_path, start, stop)
            f.write(cmd)
            f.write('\nwait\n\n')
        f.close()
        if run is True:
            home = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(home)
   
    def rmsfSidechain(self, run=True):
        '''
        Writes & submits .sh file to perform rmsf calculations
        '''
        filename = 'rmsf_{}.sh'.format(self.name)
        filename = os.path.join(self.directory['scripts'], filename)
        f = open(filename,'w')
        header = self.getHeader(job_name='{}_rmsf'.format(self.name), walltime='10:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in range(1, self.reps+1):
            if 'rmsf.ndx' not in os.listdir(self.directory[rep]['root']):
                ndx_filename = self.makeNDX(rep, custom='rmsfSidechain')
            else:
                ndx_filename = os.path.join(self.directory[rep]['root'], 'rmsf.ndx')
            ndx_filename = self.makeNDX(rep, custom='rmsfSidechain')
            if self.peptides is not None:
                for pep in range(0, self.peptides):
                    output_path = os.path.join(self.directory[rep]['rmsf']['root'], 'rmsf_{}.xvg'.format(str(pep+1)))
                    cmd = 'echo {} | gmx rmsf -res -f {} -s {} -n {} -o {} \n'.format(str(pep), self.directory[rep]['xtc_system'], self.directory[rep]['tpr'], ndx_filename, output_path)
                    f.write(cmd)
        f.close()
        if run is True:
            os.system('sbatch {}'.format(filename))
        
    def cluster(self, stop, start=None, interval=None, cutoff=0.2, rep=None, sm=True, run=True):
        '''
        Writes & submits .sh file to perform clustering
        '''
        if rep is None:
            filename = 'cluster_{}.sh'.format(self.name)
        else:
            filename = 'cluster_{}_rep{}.sh'.format(self.name, rep)
        filename = os.path.join(self.directory['scripts'], filename)
        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_cluster'.format(self.name), walltime='20:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        _rep = rep
        for rep in range(1, self.reps+1):
            if _rep is not None:
                if rep != _rep:
                    continue
            # if 'cluster.ndx' not in os.listdir(self.directory[rep]['root']):
            ndx_filename = self.makeNDX(rep, custom='cluster', sm=sm)
            # else:
            #     ndx_filename = os.path.join(self.directory[rep]['root'], 'cluster.ndx')

            if self.directory[rep]['xtc_pro_sm'] is not None:
                xtc = self.directory[rep]['xtc_pro_sm']
            elif self.directory[rep]['xtc_pro'] is not None:
                xtc = self.directory[rep]['xtc_pro']
            else:
                xtc = self.directory[rep]['xtc_system']
            if interval is None:
                for i in range(0, stop, 100000):
                    root = self.directory[rep]['clusters']['root']
                    xpm_filename = 'rmsd_clust_{}_{}.xpm'.format(str(int(i/1000)), str(int((i + 100000)/1000)))
                    xpm_filepath = os.path.join(root, xpm_filename)
                    log_filename = 'cluster_{}_{}.log'.format(str(int(i/1000)), str(int((i + 100000)/1000)))
                    log_filepath = os.path.join(root, log_filename)
                    pdb_filename = 'clusters_{}_{}.pdb'.format(str(int(i/1000)), str(int((i + 100000)/1000)))
                    pdb_filepath = os.path.join(root, pdb_filename)
                    size_filename = 'cluster_size_{}_{}.xvg'.format(str(int(i/1000)), str(int((i + 100000)/1000)))
                    size_filepath = os.path.join(root, size_filename)
                    dist_filename = 'rmsd_dist_{}_{}.xvg'.format(str(int(i/1000)), str(int((i + 100000)/1000)))
                    dist_filepath = os.path.join(root, dist_filename)
                    cmd = 'echo {} {} | gmx cluster -n {} -f {} -s {} -method gromos -o {} -minstruct 200 -g {} -cl {} -wcl 10 -cutoff 0.2 -sz {} -dist {} -b {} -e {}\n'.format(str(0), str(1), ndx_filename, xtc, self.directory[rep]['tpr'], xpm_filepath, log_filepath, pdb_filepath, size_filepath, dist_filepath, str(i), str(i+100000))
                    f.write(cmd)
                    f.write('\nwait\n\n')
            else:
                if start is not None:
                    for i in range(start, stop, interval):
                        root = self.directory[rep]['clusters']['root']
                        xpm_filename = 'rmsd_clust_{}_{}.xpm'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        xpm_filepath = os.path.join(root, xpm_filename)
                        log_filename = 'cluster_{}_{}.log'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        log_filepath = os.path.join(root, log_filename)
                        pdb_filename = 'clusters_{}_{}.pdb'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        pdb_filepath = os.path.join(root, pdb_filename)
                        size_filename = 'cluster_size_{}_{}.xvg'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        size_filepath = os.path.join(root, size_filename)
                        dist_filename = 'rmsd_dist_{}_{}.xvg'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        dist_filepath = os.path.join(root, dist_filename)
                        cmd = 'echo {} {} | gmx cluster -n {} -f {} -s {} -method gromos -o {} -minstruct 200 -g {} -cl {} -wcl 10 -cutoff {} -sz {} -dist {} -b {} -e {}\n'.format(str(0), str(1), ndx_filename, xtc, self.directory[rep]['tpr'], xpm_filepath, log_filepath, pdb_filepath, cutoff, size_filepath, dist_filepath, str(i), str(i+interval))
                        f.write(cmd)
                        f.write('\nwait\n\n')
                else:
                    for i in range(0, stop, interval):
                        root = self.directory[rep]['clusters']['root']
                        xpm_filename = 'rmsd_clust_{}_{}.xpm'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        xpm_filepath = os.path.join(root, xpm_filename)
                        log_filename = 'cluster_{}_{}.log'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        log_filepath = os.path.join(root, log_filename)
                        pdb_filename = 'clusters_{}_{}.pdb'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        pdb_filepath = os.path.join(root, pdb_filename)
                        size_filename = 'cluster_size_{}_{}.xvg'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        size_filepath = os.path.join(root, size_filename)
                        dist_filename = 'rmsd_dist_{}_{}.xvg'.format(str(int(i/1000)), str(int((i + interval)/1000)))
                        dist_filepath = os.path.join(root, dist_filename)
                        cmd = 'echo {} {} | gmx cluster -n {} -f {} -s {} -method gromos -o {} -minstruct 200 -g {} -cl {} -wcl 10 -cutoff {} -sz {} -dist {} -b {} -e {}\n'.format(str(0), str(1), ndx_filename, xtc, self.directory[rep]['tpr'], xpm_filepath, log_filepath, pdb_filepath, cutoff, size_filepath, dist_filepath, str(i), str(i+interval))
                        f.write(cmd)
                        f.write('\nwait\n\n')
            f.write('\n\n')
        f.close()
        if run is True:
            home = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            print('sbatch {}'.format(filename))
            os.chdir(home)

    def mindist(self, groups, start=None, stop=None, reps=None, run=True, parallel=True, ndxt=None, runsingle=True):
        '''
        Writes and submits .sh file for mindist.
        Groups: Tuple (group1, group2). Examples: ('sidechain','sidechain'), ('mainchain', 'mainchain'), ('residue', 'residue'), ('sidechain', 'sm'), etc
        '''
        scripts = []
        for rep in range(1, self.reps+1):
            if reps is not None:
                if isinstance(reps, list):
                    if rep in reps:
                        continue
                if isinstance(reps, int):
                    if rep == reps:
                        continue
            mindist_path = os.path.join(self.directory[rep]['root'], 'mindist')
            if not os.path.isdir(mindist_path):
                os.mkdir(mindist_path)
            paths = {
                'sidechain':os.path.join(mindist_path, 'sidechain'),
                'mainchain':os.path.join(mindist_path, 'mainchain'),
                'residue':os.path.join(mindist_path, 'residue'),
                'peptide':os.path.join(mindist_path, 'peptide'),
                'protein':os.path.join(mindist_path, 'protein')
                #'sidechain_sm':os.path.join('')
            }
            if self.ligands is not None:
                paths['sidechain_sm'] = os.path.join(mindist_path, 'sidechain_sm')
                paths['mainchain_sm'] = os.path.join(mindist_path, 'mainchain_sm')
                paths['residue_sm'] = os.path.join(mindist_path, 'residue_sm')
                paths['sm_sm'] = os.path.join(mindist_path, 'sm_sm')
                paths['protein_sm'] = os.path.join(mindist_path, 'protein_sm')
            for path in paths.values():
                if not os.path.isdir(path):
                    os.mkdir(path)
            if (groups[0] == 'sidechain') or (groups[1] == 'sidechain'):
                if (groups[0] == 'sm') or (groups[1] == 'sm'):
                    filename = '{}_rep{}_sm_sidechain_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_sm_sc_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                    out_path = paths['sidechain_sm']
                else:
                    filename = '{}_rep{}_sidechain_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_sc_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                    out_path = paths['sidechain']
            if (groups[0] == 'mainchain') or (groups[1] == 'mainchain'):
                if (groups[0] == 'sm') or (groups[1] == 'sm'):
                    filename = '{}_rep{}_sm_mainchain_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_sm_mc_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                    out_path = paths['mainchain_sm']
                else:
                    filename = '{}_rep{}_mainchain_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_mc_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                    out_path = paths['mainchain']
            if (groups[0] == 'residue') or (groups[1] == 'residue'):
                if (groups[0] == 'sm') or (groups[1] == 'sm'):
                    filename = '{}_rep{}_sm_res_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_sm_res_mindist'.format(self.name. str(rep)), walltime='60:00:00')
                    out_path = paths['residue_sm']
                else:
                    filename = '{}_rep{}_res_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_res_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                    out_path = paths['residue']
            if (groups[0] == 'protein') or (groups[1] == 'protein'):
                if (groups[0] == 'sm') or (groups[1] == 'sm'):
                    filename = '{}_rep{}_sm_protein_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_sm_protein_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                    out_path = paths['protein_sm']
                else:
                    filename = '{}_rep{}_protein_mindist.sh'.format(self.name, str(rep))
                    header = self.getHeader(job_name='{}_r{}_protein_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                    out_path = paths['protein']
            if (groups[0] == 'peptide') and (groups[1] == 'peptide'):
                filename = '{}_rep{}_peptide_mindist.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_pep_mindist'.format(self.name), walltime='60:00:00')
                out_path = paths['peptide']
            if (groups[0] == 'sm') and (groups[1] == 'sm'):
                filename = '{}_rep{}_sm_sm_mindist.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_sm_mindist'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['sm_sm']

            
            _filename = filename
            filename = os.path.join(self.directory['scripts'], filename)

            if parallel is True:
                header = ['#!/bin/bash', '\n\n']
            f = open(filename, 'w')
            for line in header:
                f.write(line)
            f.write('\n')
            f.write('export GMX_MAXBACKUP=-1\n')
            
            if parallel is True:
                string = '{}{}()'.format(self.name, rep) + r' {' + '\n'
                f.write(string)
            
            if ndxt is None:
                ndx_filename = self.makeNDX(rep, custom='mindist', groups=groups)
            else:
                ndx_filename = self.makeNDX(rep, gro=self.system.directory[rep]['gro'], custom='mindist', ndxt=os.path.abspath(ndxt), groups=groups)

            if self.directory[rep]['xtc_pro_sm'] is not None:
                xtc = self.directory[rep]['xtc_pro_sm']
            elif self.directory[rep]['xtc_pro'] is not None:
                xtc = self.directory[rep]['xtc_pro']
            else:
                xtc = self.directory[rep]['xtc_system']

            ntasks = 0
            if ((groups[0] == 'sm') or (groups[1] == 'sm')) and (groups != ('sm', 'sm')):
                g = open(ndx_filename, 'r')
                residue_labels = []
                sm_labels = []
                sm_positions = []
                residue_positions = []
                current_label = None
                i = -1
                cont = False
                for line in g:
                    if '[' in line:
                        current_label = line.strip()[1:-2].strip()
                        i += 1
                        cont = False
                    else:
                        if cont is False:
                            test_num = line.strip().split()[0]
                            if test_num in self.types[self.ligands]:
                                sm_labels.append(current_label)
                                sm_positions.append(i)
                                cont = True
                            else:
                                res_label = current_label.split('_')[0]
                                residue_labels.append(res_label)
                                residue_positions.append(i)
                                cont = True
                        else:
                            continue
                g.close()
                sm_pos = 0
                res_pos = 0
                for i in sm_positions:
                    for j in residue_positions:
                        out_base = '{}_{}_mindist.xvg'.format(sm_labels[sm_pos], residue_labels[res_pos])
                        out_filepath = os.path.join(out_path, out_base)
                        if (start is not None) and (stop is not None):
                            cmd = 'echo {} {} | gmx mindist -f {} -s {} -n {} -on {} -b {} -e {}\n\n'.format(str(i), str(j), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath, str(start), str(stop))
                            f.write(cmd)
                        else:
                            cmd = 'echo {} {} | gmx mindist -f {} -s {} -n {} -on {}\n\n'.format(str(i), str(j), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath)
                            f.write(cmd)
                        res_pos += 1
                        ntasks += 1
                    sm_pos += 1
                    res_pos = 0
            elif (groups[0] == 'peptide') and (groups[1] == 'peptide'):
                ndx_groups = []
                peptide_labels = []
                g = open(ndx_filename, 'r')
                i = 0
                for line in g:
                    if '[' in line:
                        ndx_groups.append(str(i))
                        p_label_str = 'p_{}'.format(str(i+1))
                        peptide_labels.append(p_label_str)
                        i += 1
                g.close()
                for i in range(len(ndx_groups)):
                    for j in range(len(ndx_groups)):
                        if i != j:
                            out_base = '{}_{}_mindist.xvg'.format(peptide_labels[i], peptide_labels[j])
                            out_filepath = os.path.join(out_path, out_base)
                            if (start is not None) and (stop is not None):
                                cmd = 'echo {} {} | gmx mindist -f {} -s {} -n {} -on {} -b {} -e {}\n\n'.format(str(i), str(j), self.directory[rep]['xtc_system'], self.directory[rep]['tpr'], ndx_filename, out_filepath, str(start), str(stop))
                                f.write(cmd)
                            else:
                                cmd = 'echo {} {} | gmx mindist -f {} -s {} -n {} -on {}\n\n'.format(str(i), str(j), self.directory[rep]['xtc_system'], self.directory[rep]['tpr'], ndx_filename, out_filepath)
                                f.write(cmd)  
            elif (groups == ('sm', 'sm')):
                ndx_groups = []
                labels = []
                g = open(ndx_filename, 'r')
                i = 0
                for line in g:
                    if line.startswith('['):
                        ndx_groups.append(i)
                        label = line.strip()[1:-2].strip()
                        labels.append(label)
                        i += 1
                    else:
                        continue
                g.close()
                echo_pairs = {}
                for i in range(0, len(ndx_groups)):
                    if i not in echo_pairs.keys():
                        echo_pairs[i] = {}
                    for k in range(len(ndx_groups)):
                        if k not in echo_pairs.keys():
                            echo_pairs[k] = {}
                        if k not in echo_pairs[i].keys():
                            echo_pairs[i][k] = 1
                        if i not in echo_pairs[k].keys():
                            echo_pairs[k][i] = 1
                        if (echo_pairs[k][i] == 1):
                            out_base = '{}_{}_mindist.xvg'.format(labels[i], labels[k])
                            out_filepath = os.path.join(out_path, out_base)
                            if (start is not None) and (stop is not None):
                                cmd = 'echo {} {} | gmx mindist -f {} -s {} -n {} -on {} -b {} -e {}\n\n'.format(str(i), str(k), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath, str(start), str(stop))
                                f.write(cmd)
                            else:
                                cmd = 'echo {} {} | gmx mindist -f {} -s {} -n {} -on {}\n\n'.format(str(i), str(k), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath)
                                f.write(cmd)     
                            echo_pairs[i][k] += 1
                            echo_pairs[k][i] += 1  
                            ntasks += 1
            else:
                ndx_groups = []
                res_indeces = []
                g = open(ndx_filename, 'r')
                i = 0
                for line in g:
                    if '[' in line:
                        ndx_groups.append(str(i))
                        i += 1
                        line_parts = line.split()
                        for part in line_parts:
                            parts = part.split('_')
                            if 'ri' in parts[0]:
                                res_indeces.append(parts[0])
                    else:
                        continue
                g.close()
                echo_pairs = {}
                for i in range(len(ndx_groups)):
                    if i not in echo_pairs.keys():
                        echo_pairs[i] = {}
                    for k in range(len(ndx_groups)):
                        if k not in echo_pairs.keys():
                            echo_pairs[k] = {}
                        if k not in echo_pairs[i].keys():
                            echo_pairs[i][k] = 1
                        if i not in echo_pairs[k].keys():
                            echo_pairs[k][i] = 1
                        if (echo_pairs[k][i] == 1) and (echo_pairs[i][k] == 1):
                            out_base = '{}_{}_mindist.xvg'.format(res_indeces[i], res_indeces[k])
                            out_filepath = os.path.join(out_path, out_base)
                            if (start is not None) and (stop is not None):
                                cmd = '\techo {} {} | gmx mindist -f {} -s {} -n {} -on {} -b {} -e {}\n\n'.format(str(i), str(k), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath, str(start), str(stop))
                                f.write(cmd)
                            else:
                                cmd = '\techo {} {} | gmx mindist -f {} -s {} -n {} -on {}\n\n'.format(str(i), str(k), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath)
                                f.write(cmd)
                            echo_pairs[i][k] += 1
                            echo_pairs[k][i] += 1
                            ntasks += 1
            if parallel is True:
                f.write('}')
            f.close()
            scripts.append(_filename)
            if (run is True) and (parallel is False):
                os.system('sbatch {}'.format(filename))

        if parallel is True:
            if (groups[0] == 'mainchain') or (groups[1] == 'mainchain'):
                if 'sm' not in groups:
                    filepath = os.path.join(self.directory['scripts'], '{}_mainchain_mindist.sh'.format(self.name))
                    filename = '{}_mainchain_mindist.sh'.format(self.name)
                    jobname = '{}_mainchain_mindist'.format(self.name)
                else:
                    filepath = os.path.join(self.directory['scripts'], '{}_mainchain_sm_mindist.sh'.format(self.name))
                    filename = '{}_mainchain_sm_mindist.sh'.format(self.name)
                    jobname = '{}_sm_mc_mindist'.format(self.name)
            if (groups[0] == 'residue') or (groups[1] == 'residue'):
                if 'sm' not in groups:
                    filepath = os.path.join(self.directory['scripts'], '{}_residue_mindist.sh'.format(self.name))
                    filename = '{}_residue_mindist.sh'.format(self.name)
                    jobname = '{}_residue_mindist'.format(self.name)
                else:
                    filepath = os.path.join(self.directory['scripts'], '{}_residue_sm_mindist.sh'.format(self.name))
                    filename = '{}_residue_sm_mindist.sh'.format(self.name)
                    jobname = '{}_sm_res_mindist'.format(self.name)
            if (groups[0] == 'sidechain') or (groups[1] == 'sidechain'):
                if 'sm' not in groups:
                    filepath = os.path.join(self.directory['scripts'], '{}_sidechain_mindist.sh'.format(self.name))
                    filename = '{}_sidechain_mindist.sh'.format(self.name)
                    jobname = '{}_sc_mindist'.format(self.name)
                else:
                    filepath = os.path.join(self.directory['scripts'], '{}_sidechain_sm_mindist.sh'.format(self.name))
                    filename = '{}_sidechain_sm_mindist.sh'.format(self.name)
                    jobname = '{}_sm_sc_mindist'.format(self.name)
            if (groups == ('sm', 'sm')):
                filepath = os.path.join(self.directory['scripts'], '{}_sm_sm_mindist.sh'.format(self.alias))
                filename = '{}_sm_sm_mindist.sh'.format(self.alias)
                jobname = '{}_sm_sm_mindist'.format(self.alias)
            if (groups == ('peptide', 'peptide')):
                filepath = os.path.join(self.directory['scripts'], '{}_peptide_mindist.sh'.format(self.name))
                filename = '{}_peptide_mindist.sh'.format(self.name)
                jobname = '{}_peptide_mindist'.format(self.name)   
            if (groups == ('protein', 'sm')) or (groups == ('sm', 'protein')):
                filepath = os.path.join(self.directory['scripts'], '{}_protein_sm.sh'.format(self.name))
                filename =   '{}_protein_sm.sh'.format(self.name)
                jobname = '{}_pro_sm_mindist'.format(self.name)
            ntasks = ntasks * self.reps
            if run is True:
                self.runMulti(scripts=scripts, filename=filename, filepath=filepath, job_name=jobname, nodes=2, ntasks=ntasks, run=True, runsingle=runsingle)
            else:
                self.runMulti(scripts=scripts, filename=filename, filepath=filepath, job_name=jobname, nodes=2, ntasks=ntasks, run=False, runsingle=runsingle)

    def hbonds(self, group, start=None, stop=None, reps=None, ndxt=None, run=True, parallel=True, runsingle=True):
        scripts = []
        for rep in range(1, self.reps+1):
            if reps is not None:
                if isinstance(reps, list):
                    if rep in reps:
                        continue
                if isinstance(reps, int):
                    if rep == reps:
                        continue
            hbonds_path = os.path.join(self.directory[rep]['root'], 'hbonds')
            if not os.path.isdir(hbonds_path):
                os.mkdir(hbonds_path)
            paths = {
                'sidechain_pro':os.path.join(hbonds_path, 'sidechain_pro'),
                'mainchain_pro':os.path.join(hbonds_path, 'mainchain_pro'),
                'residue_pro':os.path.join(hbonds_path, 'residue'),
                'backbone_pro':os.path.join(hbonds_path, 'backbone_pro')
            }
            if self.ligands is not None:
                paths['sidechain_sm'] = os.path.join(hbonds_path, 'sidechain_sm')
                paths['mainchain_sm'] = os.path.join(hbonds_path, 'mainchain_sm')
                paths['residue_sm'] = os.path.join(hbonds_path, 'residue_sm')
                paths['sm'] = os.path.join(hbonds_path, 'sm')
                paths['protein_sm'] = os.path.join(hbonds_path, 'protein_sm')
                paths['backbone_sm'] = os.path.join(hbonds_path, 'backbone_sm')
            for path in paths.values():
                if not os.path.isdir(path):
                    os.mkdir(path)
            if (group == 'sidechain_sm'):
                filename = '{}_rep{}_sm_sidechain_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_sm_sc_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['sidechain_sm']
                outfile_suffix = '_sm_sc'
            if (group == 'sidechain_pro'):
                filename = '{}_rep{}_pro_sidechain_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_pro_sc_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['sidechain_pro']
                outfile_suffix = '_pro_sc'
            if (group == 'mainchain_sm'):
                filename = '{}_rep{}_sm_mainchain_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_sm_sc_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['mainchain_sm']
                outfile_suffix = '_sm_mc'
            if (group == 'mainchain_pro'):
                filename = '{}_rep{}_pro_mainchain_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_pro_mc_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['mainchain_pro']
                outfile_suffix = '_pro_mc'
            if (group == 'residue_sm'):
                filename = '{}_rep{}_sm_res_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_sm_res_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['residue_sm']
                outfile_suffix = '_sm_res'
            if (group == 'protein_sm'):
                filename = '{}_rep{}_pro_sm_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_pro_sm_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['protein_sm']
                outfile_suffix = 'pro_sm'
            if (group == 'residue_pro'):
                filename = '{}_rep{}_pro_res_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_pro_res_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['residue_pro']
                outfile_suffix = '_pro_res'
            if (group == 'sm'):
                filename = '{}_rep{}_sm_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_sm_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['sm']
                outfile_suffix = None
            if (group == 'backbone_pro'):
                filename = '{}_rep{}_pro_backbone_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_pro_bb_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['backbone_pro']
            if (group == 'backbone_sm'):
                filename = '{}_rep{}_backbone_sm_hbonds.sh'.format(self.name, str(rep))
                header = self.getHeader(job_name='{}_r{}_sm_bb_hbonds'.format(self.name, str(rep)), walltime='60:00:00')
                out_path = paths['backbone_sm']



            _filename = filename
            filename = os.path.join(self.directory['scripts'], filename)
            # ndx_filename = self.makeNDX(rep, custom='hbonds', groups=group)
            if parallel is True:
                header = ['#!/bin/bash', '\n\n']
            f = open(filename, 'w')
            for line in header:
                f.write(line)
            f.write('\n')
            f.write('export GMX_MAXBACKUP=-1\n')
            string = '{}{}()'.format(self.name, rep) + r' {' + '\n' # function to source
            f.write(string)

            if ndxt is None:
                ndx_filename = self.makeNDX(rep, custom='hbonds', groups=group)
            else:
                ndx_filename = self.makeNDX(rep, gro=self.system.directory[rep]['gro'], custom='hbonds', ndxt=os.path.abspath(ndxt), groups=groups)
            echo_pairs = {}
            ntasks = 0

            ndx_groups = []
            pro_groups = []
            g = open(ndx_filename, 'r')
            first = False

            if self.directory[rep]['xtc_pro_sm'] is not None:
                xtc = self.directory[rep]['xtc_pro_sm']
            elif self.directory[rep]['xtc_pro'] is not None:
                xtc = self.directory[rep]['xtc_pro']
            else:
                xtc = self.directory[rep]['xtc_system']

            if ('pro' not in group) and (group != 'sm'):
                for line in g:
                    if '[' in line:
                        if first is False:
                            first = True
                        elif line.startswith('[ protein'):
                            continue
                        else:
                            ri = line[1:-2].strip().split('_')[0] # get residue index 
                            ndx_groups.append(ri)
            elif (group == 'sm') or (group == 'protein_sm'):
                for line in g:
                    if '[' in line:
                        ri = line[1:-2].strip()
                        ndx_groups.append(ri)
            else:
                for line in g:
                    if line.startswith('[ ri'):
                        ri = line[1:-2].strip().split('_')[0]
                        ndx_groups.append(ri)  
                for ri in ndx_groups:
                    index = int(ri[2:])
                    pro_group_index = int(ri[2:]) + len(ndx_groups)-1
                    pro_groups.append(pro_group_index)
            g.close()        

            if ('pro' not in group) and (group != 'sm'):
                for i in range(1, len(ndx_groups)+1):
                    if (start is not None) and (stop is not None):
                        out_base = '{}{}_{}_{}.xvg'.format(ndx_groups[i-1], outfile_suffix, start, stop)
                        out_filepath = os.path.join(out_path, out_base)
                        cmd = '\techo 0 {} | gmx hbond -f {} -s {} -n {} -num {} -b {} -e {}\n\n'.format(str(i), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath, str(start), str(stop))
                        f.write(cmd)
                    else:
                        out_base = '{}{}.xvg'.format(ndx_groups[i-1], outfile_suffix)
                        out_filepath = os.path.join(out_path, out_base)
                        cmd = '\techo 0 {} | gmx hbond -f {} -s {} -n {} -num {}\n\n'.format(str(i), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath)
                        f.write(cmd)
                    ntasks += 1
            elif (group == 'sm') or (group == 'protein_sm'):
                pairs = []
                for i in range(0, len(ndx_groups)):
                    for k in range(0, len(ndx_groups)):
                        if i == k:
                            continue
                        pair = (i,k)
                        if pair in pairs:
                            continue
                        else:
                            pairs.append(pair)
                        pair = (k,i)
                        if pair in pairs:
                            continue
                        else:
                            pairs.append(pair)
                        if (start is not None) and (stop is not None):
                            out_base = '{}_{}_{}_{}.xvg'.format(ndx_groups[i], ndx_groups[k], int(start/1000), int(stop/1000))
                            out_filepath = os.path.join(out_path, out_base)
                            cmd = '\techo {} {} | gmx hbond -f {} -s {} -n {} -num {} -b {} -e {}\n\n'.format(str(i), str(k), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath, str(start), str(stop))
                            f.write(cmd)
                        else:
                            out_base = '{}_{}.xvg'.format(ndx_groups[i], ndx_groups[k])
                            out_filepath = os.path.join(out_path, out_base)
                            cmd = '\techo {} {} | gmx hbond -f {} -s {} -n {} -num {}\n\n'.format(str(i), str(k), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath)
                            f.write(cmd)
            else:
                for i in range(0, len(pro_groups)):
                    if (start is not None) and (stop is not None):
                        out_base = '{}{}_{}_{}.xvg'.format(ndx_groups[i-1], outfile_suffix, start, stop)
                        out_filepath = os.path.join(out_path, out_base)
                        cmd = '\techo {} {} | gmx hbond -f {} -s {} -n {} -num {} -b {} -e {}\n\n'.format(pro_groups[i-1], str(i), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath, str(start), str(stop))
                        f.write(cmd)
                    else:
                        out_base = '{}{}.xvg'.format(ndx_groups[i], outfile_suffix)
                        out_filepath = os.path.join(out_path, out_base)
                        # print(pro_groups[i], str(i+(len(ndx_groups)-1)))
                        cmd = '\techo {} {} | gmx hbond -f {} -s {} -n {} -num {}\n\n'.format(i, str(i+(len(ndx_groups))), xtc, self.directory[rep]['tpr'], ndx_filename, out_filepath)
                        f.write(cmd)
                    ntasks += 1
            f.write('}')
            f.close()
            scripts.append(_filename)
            if (run is True) and (parallel is False):
                os.system('sbatch {}'.format(filename))

        if parallel is True:
            if group != 'sm':
                filename = _filename.split('_')[0] + '_' + _filename.split('_')[2] + '_' + _filename.split('_')[3] + '_hbonds.sh'
            else:
                filename = '{}_sm_hbonds.sh'.format(self.name)
            filepath = os.path.join(self.directory['scripts'], filename)
            jobname = '{}_{}'.format(self.name, group)
            if run is True:
                self.runMulti(scripts=scripts, filename=filename, filepath=filepath, job_name=jobname, nodes=2, ntasks=ntasks, run=True, runsingle=runsingle)
            else:
                self.runMulti(scripts=scripts, filename=filename, filepath=filepath, job_name=jobname, nodes=2, ntasks=ntasks, run=False, runsingle=runsingle)


    def hbondsSMProtein(self, group, start=None, stop=None, run=True, parallel=True, runsingle=True):
        filename = os.path.join(self.directory['scripts'], 'hbonds_fprotein_{}_sm.sh'.format(group))

        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_{}_sm_hbonds'.format(self.name, group), walltime='10:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')

        for rep in range(1, self.reps+1):
            if self.directory[rep]['xtc_pro_sm'] is not None:
                xtc = self.directory[rep]['xtc_pro_sm']
            elif self.directory[rep]['xtc_pro'] is not None:
                xtc = self.directory[rep]['xtc_pro']
            else:
                xtc = self.directory[rep]['xtc_system']

            # get ndx
            ndx_filename = os.path.join(self.system.directory[rep]['root'], 'hbonds_fprotein_{}_sm.ndx'.format(group))
            self.makeNDX(rep=rep, custom='hbonds', groups=('protein', group), filename=ndx_filename)
            output_path = os.path.join(self.system.directory[rep]['hbonds']['protein_sm']['root'], '{}_protein_{}.xvg'.format(self.system.ligands, group))
            cmd = '\techo 0 1 | gmx hbond -f {} -s {} -n {} -num {} -b 400000 -e 600000\n\n'.format(xtc, self.directory[rep]['tpr'], ndx_filename, output_path)
            f.write(cmd)
        f.close()
        if run == True:
            tmp = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(tmp)     
    def dssp(self, dssp_path, rep=None, start=None, stop=None, run=True):
        if self.cascades == True:
            if (start is None) and (stop is None):
                if rep is None:
                    filename = os.path.join(self.directory['scripts'], 'dssp_{}.sh'.format(self.name))
                else:
                    filename = os.path.join(self.directory['scripts'], 'dssp_{}_rep{}.sh'.format(self.name, rep))
            else:
                if rep is None:
                    filename = os.path.join(self.directory['scripts'], 'dssp_{}_{}_{}.sh'.format(self.name, start, stop))
                else:
                   filename = os.path.join(self.directory['scripts'], 'dssp_{}_{}_{}_rep{}.sh'.format(self.name, start, stop, rep)) 
            if filename in os.listdir(self.directory['scripts']):
                prev = 0
                for item in os.listdir(self.directory['scripts']):
                    if 'dssp' in item:
                        if 'prev' in item:
                            prev += 1
                backup = os.path.join(self.directory['scripts'], 'dssp_{}_prev{}.sh'.format(self.name, str(int(prev+1))))
                os.system('mv {} {}'.format(filename, backup))
        else:
            filename = 'dssp_{}.sh'.format(self.name)
        f = open(filename, 'w')
        if (start is None) and (stop is None):
            header = self.getHeader(job_name='{}_dssp'.format(self.name), walltime='8:00:00')
        else:
            header = self.getHeader(job_name='{}_dssp_{}_{}'.format(self.name, start, stop), walltime='8:00:00')
        for line in header:
            f.write(line)
        f.write('\n\n')
        cmd = 'export DSSP={} \n'.format(dssp_path)
        f.write(cmd)
        f.write('export GMX_MAXBACKUP=-1\n\n')
        _rep = rep
        for rep in range(1, self.reps+1):
            if _rep is not None:
                if rep != _rep:
                    continue
            dssp_path = self.directory[rep]['dssp']['root']
            if self.directory[rep]['xtc_pro_sm'] is not None:
                xtc = self.directory[rep]['xtc_pro_sm']
            elif self.directory[rep]['xtc_pro'] is not None:
                xtc = self.directory[rep]['xtc_pro']
            else:
                xtc = self.directory[rep]['xtc_system']
            tpr = self.directory[rep]['tpr']
            if (start is None) and (stop is None):
                xpm = os.path.join(dssp_path, 'dssp.xpm')
                xvg = os.path.join(dssp_path, 'dssp.xvg')
                cmd = 'echo 1 | gmx do_dssp -f {} -s {} -o {} -sc {} \n'.format(xtc, tpr, xpm, xvg)
                f.write(cmd)
                f.write('wait\n\n') 
            else:
                xpm = os.path.join(dssp_path, 'dssp_{}_{}.xpm').format(str(int(start/1000)), str(int(stop/1000)))
                xvg = os.path.join(dssp_path, 'dssp_{}_{}.xvg').format(str(int(start/1000)), str(int(stop/1000)))
                cmd = 'echo 1 | gmx do_dssp -f {} -s {} -b {} -e {} -o {} -sc {} \n'.format(xtc, tpr, start, stop, xpm, xvg)
                f.write(cmd)          
                f.write('wait\n\n')     
        f.close()
        if run is True:
            home = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(home)
    
    def gyration(self, groups, ecc=False, run=True):
        filenames = []
        for group in groups:
            if self.cascades == True:
                # filename 
                if ecc is False:
                    filename = os.path.join(self.directory['scripts'], '{}_gyration_{}.sh'.format(self.name, group))
                else:
                    filename = os.path.join(self.directory['scripts'], '{}_ecc_{}.sh'.format(self.name, group))
                if filename in os.listdir(self.directory['scripts']):
                    prev = 0
                    for item in os.listdir(self.directory['scripts']):
                        if 'gyration' in item:
                            if 'prev' in item:
                                prev += 1
                    backup = os.path.join(self.directory['scripts'], '{}_gyration_{}_prev{}.sh'.format(self.name, group, str(int(prev+1))))
                    os.system('mv {} {}'.format(filename, backup))

            #header
            if ecc is False:
                header = self.getHeader(job_name='{}_gyration_{}'.format(self.name, group), walltime='8:00:00')
            else:
                header = self.getHeader(job_name='{}_ecc_{}'.format(self.name, group), walltime='8:00:00')

            # init file
            f = open(filename, 'w')
            for line in header:
                f.write(line)
            f.write('\n\n')

            for rep in range(1, self.reps+1):
                # ndx
                if group == 'sm':
                    group = self.system.ligands
                if ecc is False:
                    ndx_filename = os.path.join(self.system.directory[rep]['root'], '{}_gyration_{}.ndx'.format(self.name, group))
                else:
                    ndx_filename = os.path.join(self.system.directory[rep]['root'], '{}_ecc_{}.ndx'.format(self.name, group))
                self.makeNDX(rep=rep, groups=[group], filename=ndx_filename)

                # files
                if self.directory[rep]['xtc_pro_sm'] is not None:
                    xtc = self.directory[rep]['xtc_pro_sm']
                elif self.directory[rep]['xtc_pro'] is not None:
                    xtc = self.directory[rep]['xtc_pro']
                else:
                    xtc = self.directory[rep]['xtc_system']
                tpr = self.directory[rep]['tpr']

                # make paths
                gyr_path = os.path.join(self.system.directory[rep]['root'], 'gyration')
                if not os.path.isdir(gyr_path):
                    os.mkdir(gyr_path)
                if ecc is False:
                    group_path = os.path.join(gyr_path, group)
                else:
                    group_path = os.path.join(gyr_path, 'ecc')
                if not os.path.isdir(group_path):
                    os.mkdir(group_path)
                
                # write
                if ecc is False:
                    output = os.path.join(group_path, '{}_gyration.xvg'.format(group))
                else:
                    output = os.path.join(group_path, '{}_ecc.xvg'.format(group))
                if ecc is False:
                    cmd = 'echo 0 | gmx gyrate -f {} -s {} -n {} -o {}\n\n'.format(xtc, tpr, ndx_filename, output)
                else:
                    cmd = 'echo 0 | gmx gyrate -f {} -s {} -n {} -o {} -moi \n\n'.format(xtc, tpr, ndx_filename, output)
                f.write(cmd)
            f.close()
            filenames.append(filename)
        if run is True:
            for filename in filenames:
                os.system('sbatch {}'.format(filename))

 
    def makeNDX(self, rep, gro=None, indeces=None, custom=None, sm=True, groups=None, filename=None, testing=False, ndxt=None):
        '''
        Makes .ndx file
        indeces: the standard index names as seen in getTypes, as a list
        custom: syntax options below
        ** can be name of function that will use the ndx (ex: rmsfSidechain)
        Testing can be passed to write a sample file to the current working directory. 
        '''    
        def chunkGenerator(lis, n):
            '''
            Internal function to split lists into groups of 15
            '''
            for i in range(0, len(lis), n):
                try:
                    yield lis[i:i+n]
                except:
                    yield lis[i:]
        
        def writeLines(indeces):
            '''
            Internal function to generate lines for NDX file.
            '''
            strings = []
            for key, value in indeces.items():
                if len(value) > 0:
                    string = '[ {} ]\n'.format(key)
                    strings.append(string)
                    string = ''
                    for line in chunkGenerator(value, 15):
                        length = 0
                        for item in line:
                            if len(str(item)) > length:
                                length = len(str(item))
                        if length < 4:
                            string = ''.join(['%4s' % i for i in line])
                        else:
                            string = ''.join(['%4s ' % i for i in line])
                        string = string + '\n'
                        strings.append(string)
            return strings

        if (custom is not None) and (ndxt is None):
            if custom == 'rmsfSidechain':
                indeces = self.rmsfSidechainNDX(self.directory[rep]['gro'])
                filename = os.path.join(self.directory[rep]['root'],'rmsf.ndx')
            if custom == 'rmsdPerPeptide':
                indeces = self.rmsdPerPeptideNDX(self.directory[rep]['gro'])
                filename = os.path.join(self.directory[rep]['root'], 'rmsd_per_peptide.ndx')
            if custom == 'cluster':
                indeces = self.clustersNDX(self.directory[rep]['gro'], sm=sm)
                filename = os.path.join(self.directory[rep]['root'],'cluster.ndx')
            if custom == 'mindist':
                if testing is True:
                    indeces = self.mindistNDX(gro, groups)
                    filename = 'test.ndx'
                else:
                    try:
                        indeces = self.mindistNDX(self.directory[rep]['gro'], groups)
                    except:
                        indeces = self.mindistNDX(gro, groups)
                    if (groups[0] == 'sidechain') or (groups[1] == 'sidechain'):
                        if (groups[0] == 'sm') or (groups[1] == 'sm'):
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_sidechain_sm.ndx') 
                        else:
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_sidechain.ndx') 
                    if (groups[0] == 'mainchain') or (groups[1] == 'mainchain'):
                        if (groups[0] == 'sm') or (groups[1] == 'sm'):
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_mainchain_sm.ndx') 
                        else:
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_mainchain.ndx') 
                    if (groups[0] == 'residue') or (groups[1] == 'residue'):
                        if (groups[0] == 'sm') or (groups[1] == 'sm'):
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_residue_sm.ndx') 
                        else:
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_residue.ndx') 
                    if (groups[0] == 'peptide') or (groups[1] == 'peptide'):
                        if self.testing is False:
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_peptide_peptide.ndx')
                        else:
                            filename = 'testing_ndx.ndx'
                    if (groups[0] == 'protein') or (groups[1] == 'protein'):
                        if (groups[0] == 'sm') or (groups[1] == 'sm'):
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_protein_sm.ndx') 
                        else:
                            filename = os.path.join(self.directory[rep]['root'], 'mindist_protein.ndx') 
                    if groups == ('sm', 'sm'):
                        filename = os.path.join(self.directory[rep]['root'], 'mindist_sm_sm.ndx') 
            if custom == 'residuePeptide':
                indeces = self.residueNDX(gro=gro)
                if filename is None:
                    filename = 'res_peptide.ndx'
            if custom == 'hbonds':
                change = False
                if (groups != 'sm') and (isinstance(groups, str)):
                    groups = (groups.split('_')[0], groups.split('_')[1])
                elif isinstance(groups, tuple):
                    change = True
                    groups = groups[1]
                else:
                    groups = ('sm', 'sm')
                if change is False:
                    indeces = self.hbondsNDX(gro=self.directory[rep]['gro'], groups=groups)
                    filename = os.path.join(self.directory[rep]['root'], 'hbonds_{}_{}.ndx'.format(groups[0], groups[1]))
                else:
                    indeces = self.hbondsSMProteinNDX(gro=self.directory[rep]['gro'], groups=groups)

        elif (groups is not None) and (ndxt is None):
            indeces = {}
            if gro is not None:
                types = self.getTypes(gro)
            else:
                types = self.getTypes(self.system.gro)
            for group in groups: # need to write a proper interpereter one of these days!!!!!
                finished = False
                operators = ['and', 'or', '&', '|']
                for op in operators:
                    logic = group.split(op)
                    if len(logic) > 1:
                        g1 = types[logic[0].strip()]
                        g2 = types[logic[1].strip()]
                        if (op == 'and') or (op == '&'):
                            label = '{}_and_{}'.format(logic[0].strip(), logic[1].strip())
                            indeces[label] = [i for i in g1 if i in g2]
                            finished = True
                            break
                        if (op == 'or') or (op == '|'):
                            label = '{}_or_{}'.format(logic[0], logic[1])
                            for item in g2:
                                g1.append(item)
                            indeces[label] = g1
                            finished=True
                            break
                if finished is False:
                    indeces[group] = types[group]
        elif ndxt is not None:
            if filename is None:
                fname = os.path.basename(ndxt).split('.')[0] + '.ndx'
                filename = os.path.join(self.directory[rep]['root'], fname)
            if (groups is not None) and ('sidechain' in groups):
                indeces = self.ndxFromNDXT(gro=self.directory[rep]['gro'], ndxt=ndxt, groups='sidechain', itype=custom)
                # for key,value in indeces.items():
                #     print(key,value)
            else:
                indeces = self.ndxFromNDXT(gro=self.directory[rep]['gro'], ndxt=ndxt, groups=None, itype=None)
        else:
            print("you're fucked LOL")
        

        f = open(filename, 'w')
        strings = writeLines(indeces)
        for string in strings:
            f.write(string)
        f.close()
        return filename

    def ndxFromNDXT(self, gro, ndxt, groups=None, itype=None):
        ndxt_groups = {}
        f = open(ndxt, 'r')
        for line in f:
            line_parts = line.split(':')
            key = line_parts[0]
            vals = line_parts[1].strip().split(',')
            ndxt_groups[key] = vals
        f.close()

        types = self.getTypes(gro=gro, ndxt_groups=ndxt_groups, set_types=True)
        indeces = {}
        indeces[self.system.ligands] = types[self.system.ligands]
        for key in ndxt_groups.keys():
            indeces[key] = types[key]
        if (groups is not None) and (itype is not None):
            if itype == 'mindist':
                group_indeces = self.mindistNDX(gro=gro, groups=(groups, groups))
            # TODO: write case for hbonds
                for key in group_indeces.keys():
                    indeces[key] = group_indeces[key]
        return indeces

    def getIndecesByGroup(self, gro, group):
        pass
        # indeces = {}
        # types = self.getTypes(gro)
        # if len(group) == 2:
        #     if 'peptide' in group:
        #     group0 = types[group[0]]
        #     group1 = types[group[1]]

    def mindistNDX(self, gro, groups):
        '''
        Write NDX file for specific mindist groups.
        '''
        print(groups)
        types = self.getTypes(gro)
        indeces = {}
        if (groups[0] == 'sidechain') or (groups[1] == 'sidechain'):
            resgroup = types['sidechain_h']
            label_add = '_and_Sidechain'
        if (groups[0] == 'mainchain') or (groups[1] == 'mainchain'):
            resgroup = types['mainchain_h_nocaps']
            label_add = '_and_Mainchain'
        if (groups[0] == 'residue') or (groups[1] == 'residue'):
            resgroup = types['protein']
            label_add = None
        if (groups[0] == 'protein') or (groups[1] == 'protein'):
            indeces['protein'] = types['protein']
        if (groups[0] == 'sm') or (groups[1] == 'sm'):
            indeces['sm'] = types[self.system.ligands]
        if (groups[0] == 'peptide') and (groups[1] == 'peptide'):
            for key, atoms in types.items():
                if ('p' in key) and (len(key) == 2):
                    indeces[key] = atoms
            return indeces
        elif 'protein' in groups:
            return indeces
        elif groups == ('sm', 'sm'):
            return indeces
        else:
            for _type in types.keys():
                if 'ri' in _type:
                    atoms = types[_type]
                    if label_add is not None:
                        label = _type + label_add
                    else:
                        label = _type
                    indeces[label] = [value for value in atoms if value in resgroup]
            return indeces

    def interactionsNDX(self, gro, groups, itype):
        '''
        Write NDX file for specific mindist groups.
        '''
        types = self.getTypes(gro)
        indeces = {}
        glycine_indeces = []
        for chain in self.system.protein.chains.keys(): # get glycine location info
            i = 0
            for ndx in self.protein.chains[chain]['indeces']:
                if 'GLY' in self.protein.chains[chain]['ids'][i]:
                    glycine_indeces.append(ndx)
                i += 1
        if itype == 'hbonds':
            print(groups)
            if groups == ('sm', 'sm'):
                for key in types.keys():
                    if (self.ligands in key) and (key != self.ligands):
                        indeces[key] = types[key]
                return indeces
            if 'pro' in groups:
                indeces['protein'] = types['protein']
            if 'sm' in groups:
                indeces['sm'] = types[self.ligands]
        if (groups[0] == 'sidechain') or (groups[1] == 'sidechain'):
            resgroup = types['sidechain_h']
            label_add = '_and_Sidechain'
        if (groups[0] == 'mainchain') or (groups[1] == 'mainchain'):
            resgroup = types['mainchain_h_nocaps']
            label_add = '_and_Mainchain'
        if (groups[0] == 'residue') or (groups[1] == 'residue'):
            resgroup = types['protein']
            label_add = None
        if (groups[0] == 'sm') or (groups[1] == 'sm'):
            indeces['sm'] = types['nonprotein']
        if (groups[0] == 'peptide') and (groups[1] == 'peptide'):
            for key, atoms in types.items():
                if ('p' in key) and (len(key) == 2):
                    indeces[key] = atoms
            return indeces
        else:
            for _type in types.keys():
                if 'ri' in _type:
                    atoms = types[_type]
                    if label_add is not None:
                        label = _type + label_add
                    else:
                        label = _type
                    indeces[label] = [value for value in atoms if value in resgroup]
            if (itype == 'hbonds') and ('pro' in groups):
                for _type in types.keys():
                    if 'ri' in _type:
                        if 'sidechain' in groups:
                            current_ndx = int(_type[2:])
                            if current_ndx in glycine_indeces:
                                continue
                        atoms = types[_type]
                        label = 'protein_not_{}'.format(_type)
                        indeces[label] = [value for value in indeces['protein'] if value not in atoms]
                del indeces['protein']
            return indeces

    def hbondsNDX(self, gro, groups):
        types = self.getTypes(gro)
        indeces = {}
        if 'sm' in groups:
            if isinstance(self.ligands, str):
                indeces[self.ligands] = types[self.ligands]
            else:
                for ligand in self.ligands:
                    indeces[ligand] = types[ligand]
        if 'sidechain' in groups:
            resgroup = types['sidechain_h']
            label_add = '_and_Sidechain'
        if 'mainchain' in groups:
            resgroup = types['mainchain_h']
            label_add = '_and_Mainchain'
        if 'residue' in groups:
            resgroup = types['protein']
            label_add = None
        if 'backbone' in groups:
            resgroup = types['backbone']
        indeces['protein'] = [value for value in types['protein'] if value in resgroup]
        if 'protein' not in groups:
            for _type in types.keys():
                if 'ri' in _type:
                    atoms = types[_type]
                    if label_add is not None:
                        label = _type + label_add
                    else:
                        label = _type
                    indeces[label] = [value for value in atoms if value in resgroup]
        else:
            indeces['protein'] = types['protein']
        return indeces

    def hbondsSMProteinNDX(self, gro, groups):
        indeces = {}
        types = self.getTypes(gro)
        if isinstance(self.ligands, str):
            indeces[self.ligands] = types[self.ligands]
        else:
            for ligand in self.ligands:
                indeces[ligand] = types[ligand]
        if 'sidechain' in groups:
            resgroup = types['sidechain_h']
            label_add = '_and_Sidechain'
        if 'mainchain' in groups:
            resgroup = types['mainchain_h']
            label_add = '_and_Mainchain'
        if 'residue' in groups:
            resgroup = types['protein']
            label_add = None
        if 'backbone' in groups:
            resgroup = types['backbone']
            label_add = '_and_Backbone'
        if label_add is not None:
            label = 'protein' + label_add
        else:
            label = _type
        indeces[label] = [value for value in types['protein'] if value in resgroup]
        return indeces
    def rmsfSidechainNDX(self, gro):
        types = self.getTypes(gro)
        sidechain = types['sidechain_h']
        indeces = {}
        if self.peptides is None:
            protein = types['protein']
            indeces['protein_sidechain_H'] = [value for value in protein if value in sidechain]
        else:
            for key,atoms in types.items():
                if ('p' in key) and len(key) == 2:
                    new_key = '{}_sidechain_H'.format(key)
                    indeces[new_key] = [value for value in atoms if value in sidechain]
        return indeces

    def clustersNDX(self, gro, sm):
        print(sm)
        types = self.getTypes(gro)
        indeces = {}
        indeces['backbone'] = types['backbone']
        protein_sm = []
        for item in types['protein']:
            protein_sm.append(item)
        for item in types['caps']:
            protein_sm.append(item)
        if sm == True:
            for item in types[self.system.ligands]:
                protein_sm.append(item)
                print(item)
        protein_sm = [int(i) for i in protein_sm]
        protein_sm = sorted(protein_sm)
        protein_sm = [str(i) for i in protein_sm]
        indeces['protein_sm'] = protein_sm
        return indeces

    def rmsdPerPeptideNDX(self, gro):
        types = self.getTypes(gro)
        indeces = {}
        backbone = types['backbone']
        for key, atoms in types.items():
            if ('p' in key) and len(key) == 2:
                new_key = '{}_backbone'.format(key)
                indeces[new_key] = [value for value in atoms if value in backbone]
        return indeces

    def residueNDX(self, gro):
        residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD']
        lipids = ['POPC', 'CHL1', 'SDPE', 'POPE', 'PSM', 'SOPS']
        types = self.getTypes(gro)
        indeces = {}
        residue_indeces = {}
        peptide_indeces = {}
        for key, value in types.items():
            if (len(key)) > 3:
                if key[-3:] in residues:
                    residue_indeces[key] = value
            if ('p' in key) and (len(key) == 2):
                peptide_indeces[key] = value
            if key in lipids:
                indeces[key] = value
        for peptide_id, p_indeces in peptide_indeces.items():
            for residue_id, r_indeces in residue_indeces.items():
                new_key = '{}_{}'.format(peptide_id, residue_id)
                indeces[new_key] = [atom for atom in r_indeces if atom in p_indeces]
        return indeces

    def smPeptideNDX(self, gro):
        pass
        # types = self.getTypes(gro)
        # c = 0
        # o = 0
        # h = 0
        # indeces = {}
        # for key, atoms in types.items():
        #     if 'ri' in key:
        #         indeces[key] = atoms
        #     if 'key' == 'sm':
        #         for atom in atoms:
                    
    def lineGenerator(self, filename):
        '''
        For memory saving, generates lines for any given file
        '''
        f = open(filename, 'r')
        contents = f.readlines()
        f.close()
        for line in contents:
            yield line

    def getTypes(self, gro, ndxt_groups=None, set_types=False):
        '''
        Assigns atoms as/by:
        ** backbone
        ** mainchain,  mainchain_h, mainchain_nocaps
        ** sidechain, sidechain_h
        ** caps
        ** nonprotein (nonprotein really specifies small molecules, lumped together into one index)
        ** ions
        ** solvent
        ** residue index
        ** residue id
        ** residue name
        ** small molecule
        Returns:
        ** dict {str:[str]}
        '''
        residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD']
        caps = ['ACE', 'NH2']
        lipids = ['POPC', 'CHL1', 'SDPE', 'POPE', 'PSM', 'SOPS']
        backbone = ['CA', 'C', 'N']
        mainchain = ['CA', 'C', 'O', 'N']
        ions = ['K', 'CL', 'NA']
        solvent = ['SOL', 'TIP3']
        types = {
            'protein':[],
            'protein_caps':[],
            'backbone':[],
            'backbone_nocaps':[],
            'mainchain':[],
            'mainchain_nocaps':[],
            'mainchain_h':[],
            'mainchain_h_nocaps':[],
            'sidechain':[],
            'sidechain_h':[],
            'caps':[],
            'nonprotein':[],
            'ions':[],
            'solvent':[],
            'lipids':[]
        }
        i = 0
        k = 0
        ligands = []
        if self.system.ligands is not None:
            if isinstance(self.system.ligands, str):
                types[self.system.ligands] = []
                ligands.append(self.system.ligands)
            if isinstance(self.system.ligands, list):
                for l in self.system.ligands:
                    types[l] = []
                    ligands.append(l)
        if ndxt_groups is not None:
            to_del = []
            for key in ndxt_groups.keys():
                if key not in types.keys():
                    types[key] = []
                else:
                    to_del.append(key)
            if to_del != []:
                for item in to_del:
                    del ndxt_groups[item]
        restarted = False
        atom_index = 100000
        last_res_id = None
        visited_residues = []
        peptides = 0
        for line in self.lineGenerator(gro):
            if i < 2:
                i += 1
                continue
            else:
                line_parts = line.split()
                if len(line_parts) <= 3:
                    break
                
                res_id = line[0:10].strip()
                atom_name = line[8:15].strip()
                res_name = line[5:10].strip()
                res_num = line[:5].strip()
                atom_num = line[15:20].strip()
                if atom_num == '0':
                    atom_num = str(atom_index)
                    restarted = True
                    atom_index += 1
                elif restarted is True:
                    atom_num = str(atom_index)
                    atom_index += 1
                else:
                    pass
                # print(res_id, atom_name, atom_num, res_name)
                # if not res_name.isalpha():
                #     if (res_name != 'NH2') or (res_name !=):
                #         res_name = res_id[-2:]
                if res_name in residues:
                    # protein atoms
                    types['protein'].append(atom_num)
                    types['protein_caps'].append(atom_num)

                    # backbone
                    if atom_name in backbone:
                        types['backbone'].append(atom_num)
                        types['backbone_nocaps'].append(atom_num)

                    # mainchain
                    if atom_name in mainchain:
                        types['mainchain'].append(atom_num)
                        types['mainchain_nocaps'].append(atom_num)
                        types['mainchain_h'].append(atom_num)
                        types['mainchain_h_nocaps'].append(atom_num)

                    # hydrogens
                    if (atom_name not in backbone) and (atom_name not in mainchain):
                        if atom_name == 'H':
                            types['mainchain_h'].append(atom_num)
                            types['mainchain_h_nocaps'].append(atom_num)
                        elif 'H' in atom_name:
                            types['sidechain_h'].append(atom_num)
                        else:
                            types['sidechain'].append(atom_num)
                            types['sidechain_h'].append(atom_num)

                    # res id 
                    if res_id not in list(types.keys()):
                        types[res_id] = []
                        types[res_id].append(atom_num)
                    else:
                        types[res_id].append(atom_num)

                    # residue index
                    if last_res_id != res_id:
                        k += 1
                        ri = 'ri{}'.format(str(k))
                        types[ri] = []
                        types[ri].append(atom_num)
                    else:
                        types[ri].append(atom_num)

                    # res name
                    if res_name not in list(types.keys()):
                        types[res_name] = []
                        types[res_name].append(atom_num)
                    else:
                        types[res_name].append(atom_num)

                    # peptide
                    if self.peptides is not None:
                        if last_res_id is None:
                                peptides += 1
                                pi = 'p{}'.format(str(peptides))
                                types[pi] = []
                                types[pi].append(atom_num)
                        else:
                            if res_id == visited_residues[0]:
                                if res_id != last_res_id:
                                    peptides += 1
                                    pi = 'p{}'.format(str(peptides))
                                    types[pi] = []
                                    types[pi].append(atom_num)
                                else:
                                    types[pi].append(atom_num)
                            else:
                                types[pi].append(atom_num)
                    if res_id not in visited_residues:
                        visited_residues.append(res_id)
                    last_res_id = res_id
                # caps
                elif res_name in caps:
                    types['protein_caps'].append(atom_num)
                    types['caps'].append(atom_num)
                    if atom_name in backbone:
                        types['backbone'].append(atom_num)
                    if atom_name in mainchain:
                        types['mainchain'].append(atom_num)
                    if 'H' in atom_name:
                        types['mainchain_h'].append(atom_num)
                # lipids
                elif res_name in lipids:
                    types['lipids'].append(atom_num)
                    if res_name not in types.keys():
                        types[res_name] = [atom_num]
                    else:
                        types[res_name].append(atom_num)
                # ions
                elif res_name in ions:
                    types['ions'].append(atom_num)
                # solvent
                elif res_name in solvent:
                    # types['solvent'].append(atom_num)
                    break
                # ligands
                elif res_name in ligands:
                    types[res_name].append(atom_num)
                    if res_id not in list(types.keys()):
                        types[res_id] = []
                        types[res_id].append(atom_num)
                    else:
                        types[res_id].append(atom_num)
                # other non protein
                else:
                    types['nonprotein'].append(atom_num)
                    if res_name not in list(types.keys()):
                        types[res_name] = []
                        types[res_name].append(atom_num)
                    else:
                        types[res_name].append(atom_num)
                # ndxt groups
                if ndxt_groups is not None:
                    for key in ndxt_groups.keys():
                        if atom_name in ndxt_groups[key]:
                            if atom_num not in types[key]:
                                types[key].append(atom_num)
                        if atom_num in ndxt_groups[key]:
                            if atom_num not in types[key]:
                                types[key].append(atom_num)
        # if set_types is True:
        #     self.types = types 
        self.types = types      
        return types

    def runMulti(self, scripts, filename, filepath, job_name, nodes, ntasks, run=True, runsingle=False):
        t = ((ntasks*5)/60)*self.reps
        if t > 144:
            walltime='144:00:00'
        else:
            wt = round(t)
            walltime = '80:00:00'
        header = self.getHeader(job_name=job_name, walltime=walltime, nodes=2, ntasks_per_node=None)
        f = open(filepath, 'w')
        for line in header:
            f.write(line)

        if runsingle is False:
            i = 0
            for script in scripts:
                if i % 2 == 0:
                    cmd = '. {}\n'.format(script)
                    f.write(cmd)
                    cmd = '. {}\n'.format(scripts[i+1])
                    f.write(cmd)
                    cmd = '{}{} & {}{}'.format(self.name, i+1, self.name, i+2)
                    f.write(cmd)
                else:
                    f.write('\nwait\n')
                i += 1
        else:
            i = 1
            for script in scripts:
                cmd = '. {}\n'.format(script)
                f.write(cmd)
                cmd = '{}{}'.format(self.name, i)
                f.write(cmd)
                f.write('\n\nwait\n\n')
                i += 1
        f.close()

        if run is True:
            home = os.getcwd()
            os.chdir(self.directory['scripts'])
            os.system('sbatch {}'.format(filename))
            os.chdir(home)


# r = Run(system=None, peptides=3)
# r.makeNDX(rep=None, gro='900_1000.part0011.gro', custom='residuePeptide', filename='residue_peptide_lipids.ndx')
# types = r.getTypes(gro='900_1000.part0011.gro')
# for key in types.keys():
#     if not key.startswith('ri'):
#         print(key)
# r.residueNDX(gro='900_1000.part0011.gro')
# r.makeNDX(rep=None, custom='residuePeptide', gro='900_1000.part0011.gro')
