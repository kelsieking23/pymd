import os
import sys

if os.name == 'posix':
    sys.path.append('/work/cascades/kelsieking23/iapp_analysis/scripts/python/')
    from mdanalysis.system import System
    from structure.protein import Protein
    from utilities.rewritepdb import writePDB
    from utilities.rewritepdb import editChainIDResidue
    from mdanalysis.run import Run

class MD:

    def __init__(self, root, reps, version, ff, email, name=None, peptides=None):
        self.root = root
        self.reps = reps
        self.version = version
        self.ff = ff
        self.email = email
        self.name = name
        self.peptides = peptides

        self.submission = os.path.join(self.root, 'submission')
        if not os.path.isdir(self.submission):
            os.mkdir(self.submission)


        if os.name == 'posix':
            self.BASE_PATH = '/work/cascades/kelsieking23/iapp_analysis/scripts/python/'
        self.headers = {
            '2020.3':os.path.join(self.BASE_PATH, 'mdsetup', 'base_2020_3.sh')
        }

        self.NVT = os.path.join(self.root, 'NVT')
        self.NPT = os.path.join(self.root, 'NPT')
        self.MDrun = os.path.join(self.root, 'MDrun')
        self.mdps =[]
        self.mdp_dir = os.path.join(self.BASE_PATH, 'mdsetup', 'mdp', self.version, self.ff)

    def getHeader(self, job_name, walltime, nodes=2, ntasks_per_node=24):
        header = []
        base = self.headers[self.version]
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

    def nvt(self):
        if self.version == '4.6.5':
            f = open('465base.sh', 'r')
        else:
            f = open('base.sh', 'r')
        header = f.readlines()
        f.close()
        for rep in self.reps:
            filename = os.path.join(self.submission, str(rep), 'nvt_{}.sh'.format(str(rep)))
            f = open(filename, 'w')
            for line in header:
                f.write(line)
                f.write('\n')
            f.write('\n\n')
            f.write('cd ../../../NVT/{}\n'.format(str(rep)))
            # grompp 
            if self.version == '4.6.5':
                cmd = 'grompp -f ../nvt.mdp -c ../em.gro -p ../topol.top -o nvt.tpr -po out.mdp\n'
                f.write(cmd)
                cmd = 'mdrun -c nvt.gro -s nvt.tpr -mp ../../EM/topol.top -nt 12 -deffmn nvt\n'
                f.write(cmd)
            f.close()
    
    def mdrun(self, stop):
        # get topology, itps from nvt
        os.system('cp {}/*itp {}'.format(self.NVT, self.MDrun))
        os.system('cp {} {}'.format(os.path.join(self.NVT, 'topol.top'), self.MDrun))

        # get data from NPT for each rep 
        npt = NPT(root=self.NPT, reps=self.reps)
        
        # get mdp file 
        mdp = os.path.join(self.mdp_dir, 'md.mdp')
        os.system('cp {} {}'.format(mdp, self.MDrun))

        # write files
        for rep in range(1, self.reps+1, 2):
            start = 0
            end = 25
            last_start = None
            last_end = None
            while not (end > stop):
                if start == 0:
                    walltime='30:00:00'
                else:
                    walltime='80:00:00'
                
                rep2 = rep + 1
                pair = '{}{}'.format(rep, rep2)
                filename = os.path.join(self.submission, 'md_{}_{}_{}.sh'.format(pair, start, end))
                f = open(filename, 'w')

                header = self.getHeader(job_name='{}{}_{}_{}_md'.format(self.name, pair, start, end), walltime=walltime)
                for line in header:
                    f.write(line)

                # rep 1
                f.write('\n\n')
                f.write('cd ../MDrun/{}/\n'.format(rep))
                f.write('\n')

                # grompp rep 1
                if start == 0:
                    npt_index = rep - 1
                    npt_gro = npt.gro[npt_index]
                    npt_cpt = npt.cpt[npt_index]
                    tpr = 'md_{}_{}.tpr'.format(start, end)
                    cmd = 'gmx grompp -f ../md.mdp -c {} -r {} -t {} -p ../topol.top -o {} -maxwarn 1\n'.format(npt_gro, npt_gro, npt_cpt, tpr)
                    f.write(cmd)
                else:
                    # move files
                    newdir = '{}_{}ns'.format(last_start, last_end)
                    cmd = 'mkdir {}\n'.format(newdir)
                    f.write(cmd)
                    last_md = 'md_{}_{}'.format(last_start, last_end)
                    cmd = 'mv {}* {}\n\n'.format(last_md, newdir)
                    f.write(cmd)
                    # convert tpr
                    last_tpr = 'md_{}_{}.tpr'.format(last_start, last_end)
                    last_tpr = os.path.join(newdir, last_tpr)
                    tpr = 'md_{}_{}.tpr'.format(start, end)
                    nsteps = self.nsteps(end)
                    cmd = 'gmx convert-tpr -s {} -o {} -nsteps {}\n'.format(last_tpr, tpr, nsteps)
                    f.write(cmd)


                # mdrun rep 1
                deffnm = 'md_{}_{}'.format(start, end)
                if start == 0:
                    cmd = 'mdrun_gpu -gpu_id 0 -nt 24 -mp ../topol.top -s {} -deffnm {} & \n'.format(tpr, deffnm)
                    f.write(cmd)
                else:
                    cpi = os.path.join(newdir, '{}.cpt'.format(last_md))
                    cmd = 'mdrun_gpu -gpu_id 0 -nt 24 -cpi {} -mp ../topol.top -deffnm {} -noappend &\n'.format(cpi, deffnm)
                    f.write(cmd)
                
                # rep 2
                f.write('\n\n')
                f.write('cd ../{}/\n'.format(rep2))
                f.write('\n')

                # grompp rep 2
                if start == 0:
                    npt_index = rep2 - 1
                    npt_gro = npt.gro[npt_index]
                    npt_cpt = npt.cpt[npt_index]
                    tpr = 'md_{}_{}.tpr'.format(start, end)
                    cmd = 'gmx grompp -f ../md.mdp -c {} -r {} -t {} -p ../topol.top -o {} -maxwarn 1\n'.format(npt_gro, npt_gro, npt_cpt, tpr)
                    f.write(cmd)
                else:
                    # move files
                    newdir = '{}_{}ns'.format(last_start, last_end)
                    cmd = 'mkdir {}\n'.format(newdir)
                    f.write(cmd)
                    last_md = 'md_{}_{}'.format(last_start, last_end)
                    cmd = 'mv {}* {}\n\n'.format(last_md, newdir)
                    f.write(cmd)
                    # convert tpr
                    last_tpr = 'md_{}_{}.tpr'.format(last_start, last_end)
                    last_tpr = os.path.join(newdir, last_tpr)
                    tpr = 'md_{}_{}.tpr'.format(start, end)
                    nsteps = self.nsteps(end)
                    cmd = 'gmx convert-tpr -s {} -o {} -nsteps {}\n'.format(last_tpr, tpr, nsteps)
                    f.write(cmd)

                # mdrun rep 2
                deffnm = 'md_{}_{}'.format(start, end)
                if start == 0:
                    cmd = 'mdrun_gpu -gpu_id 1 -pinoffset 24 -pinstride 1 -mp ../topol.top -s {} -deffnm {}\n'.format(tpr, deffnm)
                    f.write(cmd)
                else:
                    cpi = os.path.join(newdir, '{}.cpt'.format(last_md))
                    cmd = 'mdrun_gpu -gpu_id 1 -pinoffset 24  -pinstride 1 -cpi {} -mp ../topol.top -deffnm {} -noappend\n'.format(cpi, deffnm)
                    f.write(cmd)


                f.write('\n\nwait\n')
                
                if end != stop:
                    f.write('cd $SLURM_SUBMIT_DIR\n')
                    if start == 0:
                        next_script = 'md_{}_25_100.sh'.format(pair)
                    elif start == 25:
                        next_script = 'md_{}_100_200.sh'.format(pair)
                    else:
                        next_start = start + 100
                        next_end = end + 100
                        next_script = 'md_{}_{}_{}.sh'.format(pair, next_start, next_end)
                    cmd = 'sbatch {}'.format(next_script)
                    f.write(cmd)
                f.close()

                last_start = start
                last_end = end

                if start == 0:
                    start = 25
                elif start == 25:
                    start = 100
                else:
                    start += 100
                
                if end == 25:
                    end = 100
                else:
                    end += 100
    
    @staticmethod
    def nsteps(ns):
        return 500000*ns

class NPT:

    def __init__(self, root, reps):
        self.root = root
        self.reps = reps

    @property
    def gro(self):
        gro = []
        for rep in range(1, self.reps+1):
            gro.append(os.path.join(self.root, str(rep), 'npt.gro'))
        return gro

    @property
    def cpt(self):
        cpt = []
        for rep in range(1, self.reps+1):
            cpt.append(os.path.join(self.root, str(rep), 'npt.cpt'))
        return cpt



