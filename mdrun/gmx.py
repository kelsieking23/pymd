import os
import pymd
DIRECTORY = pymd.dir_path

class GMX:

    def __init__(self, system, testing):
        self.system = system
        self.testing = testing
        self.__dict__.update(self.system.__dict__)

    def getHeader(self, job_name, walltime, nodes=1, ntasks_per_node=24):
        header = []
        vrsn = self.system.source + '.sh'
        base = os.path.join(DIRECTORY, 'mdrun', 'sh',vrsn)
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
                    line = line.strip() + self.system.email + '\n'
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
    
    def pbc(self, inp='system', job_name='pbc', walltime='30:00:00', nodes=1, ntasks_per_node=24, output='pbc.xtc', run=False):
        '''
        Writes & submits .sh file to create system .xtc (cat_pbc.xtc) from cat.xtc
        '''
        if output is None:
            output = 'pbc.xtc'
        if self.testing == True:
            sh = os.path.join(self.system.scripts, '{}.pbc.sh'.format(self.system.name))
        else:
            sh = os.path.join(self.system.scripts, '{}.pbc.sh'.format(self.system.name))
        f = open(sh, 'w')
        header = self.getHeader(job_name='{}.{}'.format(self.system.name, job_name), walltime=walltime, nodes=nodes, ntasks_per_node=ntasks_per_node)
        for line in header:
            f.write(line)
        f.write('\n\n')
        print('+++++++++++++++++++++++')
        for rep in self.system._reps:
            if inp == 'system':
                inp = self.system.xtc
            print('***********')
            print(os.path.join(rep.root, inp))
            print(rep.tpr)
            print('***********')
            _output = os.path.join(rep.root, output)
            cmd = 'echo 0 0 | gmx trjconv -f {} -s {} -pbc mol -ur compact -center -o {} \n'.format(os.path.join(rep.root, inp), rep.tpr, _output)
            f.write(cmd)
        f.close()
        print('+++++++++++++++++++++++')
        print('Wrote {}'.format(sh))
        if run == True:
            tmp = os.getcwd()
            os.chdir(self.system.scripts)
            os.system('sbatch {}'.format(sh))
            os.chdir(tmp)
        return sh

    def stripwat(self, inp='system', job_name='stripwat', walltime='30:00:00', nodes=1, ntasks_per_node=24, output='nowat.xtc', run=False):
        if output is None:
            output = 'nowat.xtc'
        if self.testing == True:
            sh = os.path.join(self.system.scripts, '{}.stripwat.sh'.format(self.system.name))
        else:
            sh = os.path.join(self.system.scripts, '{}.stripwat.sh'.format(self.system.name))
        f = open(sh, 'w')
        header = self.getHeader(job_name='{}.{}'.format(self.system.name, job_name), walltime=walltime, nodes=nodes, ntasks_per_node=ntasks_per_node)
        for line in header:
            f.write(line)
        f.write('\n\n')
        for rep in self.system._reps:
            if inp == 'system':
                inp = self.system.xtc
            _output = os.path.join(rep.root, output)
            cmd = 'echo 1 | gmx trjconv -f {} -s {} -o {} \n'.format(os.path.join(rep.root, inp), rep.tpr, _output)
            f.write(cmd)
        f.close()
        print('Wrote {}'.format(sh))
        if run == True:
            tmp = os.getcwd()
            os.chdir(self.system.scripts)
            os.system('sbatch {}'.format(sh))
            os.chdir(tmp)
        return sh

    #TODO::edit this so that reflfects new structure
    def rmsd(self, group, ndx, run=True, start=None, stop=None):
        if self.cascades == True:
            filename = os.path.join(self.directory['scripts'], 'rmsd_{}_{}.sh'.format(group, self.name))
        else:
            filename = 'rmsd_{}_{}.sh'.format(group, self.name)


        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_{}_rmsd'.format(self.name, group), walltime='10:00:00')
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
            ndx_filename = os.path.join(self.system.directory[rep]['root'], ndx)
            ndx_num = 0
            g = open(ndx_filename, 'r')
            for line in g:
                if line.startswith('['):
                    if group not in line:
                        ndx_num += 1
                    else:
                        break
            g.close()
            output_path = os.path.join(self.system.directory[rep]['rmsd']['root'], '{}.xvg'.format(group))
            cmd = 'echo {} {} | gmx rms -f {} -s {} -n {} -o {} -tu ns \n'.format(ndx_num, ndx_num, xtc, self.directory[rep]['tpr'], ndx_filename, output_path)
            f.write(cmd)
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

    def dssp(self, dssp_path, rep=None, start=None, stop=None, run=True):
        if self.testing == False:
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