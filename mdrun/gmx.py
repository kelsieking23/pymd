import os
import pymd
DIRECTORY = pymd.dir_path

class GMX:

    def __init__(self, system, testing):
        self.system = system
        self.testing = testing

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