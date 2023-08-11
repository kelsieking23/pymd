
import os
import sys
import shutil
import subprocess
import pandas as pd
from pymd.mdsetup.mdsetup_dir import scriptdir
from pymd.structure.protein import Protein
# if sys.platform == 'win32':
#     sys.path.append('D:/Work')
#     import pymd
#     from pymd.mdsetup.mdsetup_dir import scriptdir
# if sys.platform == 'linux':
#     sys.path.append('/work/kelsieking23/software/')
#     from pymd.mdsetup.mdsetup_dir import scriptdir
    # sys.path.append('/mnt/d/Work')
    # from pymd.mdsetup.mdsetup_dir import scriptdir

# export MODULEPATH=$MODULEPATH:/projects/bevanlab/software/tinkercliffs/modules/modules/tinkercliffs-rome/all 
class Setup:
    def __init__(self, cluster, root, reps, home, email, partition, modulepath, version='2020.3', ff='charmm36m', top=None, name=None, interactive=True):
        '''
        cluster (str): computing cluster ('cascades', 'infer', etc)
        root (str): path to directory to setup MD in
        top: path to directory containing topology files
        '''
        self.cluster = cluster
        self.version = version
        self._ff = ff
        self.reps = reps
        self.email = email
        self.partition = partition
        self.modulepath = modulepath
        try:
            self.root = os.path.abspath(root)
        except:
            pass
        self.name = name
        self.home = home
        self.interactive = interactive
        self.top = top
        self.directory = self.populateDirectory()
        self.buildpath = self.directory['build']['root']

    @property
    def ff(self):
        if self._ff == 'charmm36m':
            return 'charmm36-mar2019'
    @property
    def modules(self):
        if self.cluster == 'cascades':
            return ['gcc/7.3.0', 'cuda/9.2.148', 'Anaconda', f'source /groups/bevanlab/software/cascades/gromacs/{self.version}/bin/GMXRC']
        if self.cluster == 'infer':
            return ['site/infer/easybuild/setup', 'CMake/3.15.3-GCCcore-8.3.0', 'CUDA/10.1.243-GCC-8.3.0', 'GROMACS']
        if self.cluster == 'tinkercliffs':
            return ['apps site/tinkercliffs/easybuild/setup', 'apps site/tinkercliffs-rome_a100/easybuild/setup', 'CUDA/10.1.243-iccifort-2019.5.281', 'CMake/3.15.3-GCCcore-8.3.0', f'source /projects/bevanlab/software/tinkercliffs/gromacs/{self.version}/bin/GMXRC']
    
    @property
    def mdrunPrefix(self):
        if self.cluster == 'cascades':
            return 'mdrun_gpu'
        elif self.cluster == 'tinkercliffs':
            return 'mdrun_gpu'
        else:
            return 'gmx mdrun'

    @property
    def ffSelectOption(self):
        if (self.ff == 'charmm36-mar2019'):
            if (self.cluster == 'cascades') or (self.cluster == 'tinkercliffs') or (self.cluster == 'infer'):
                return 9
            else:
                return 1  
        else:
            return 14
    def populateDirectory(self):
        '''
        Make directories: em, nvt, npt, mdrun for n reps
        '''

        directory = {}
        dirs = ['EM', 'build', 'NVT', 'NPT', 'MDrun']
        print('building directories...')
        for _dir in dirs:
            # make EM, NVT, NPT, MDrun folders
            directory[_dir] = {}
            run_dir = os.path.join(self.root, _dir)
            if os.path.isdir(run_dir):
                directory[_dir]['root'] = run_dir
            elif os.path.isdir(run_dir.lower()):
                directory[_dir]['root'] = run_dir.lower()
            else:
                os.mkdir(run_dir)
                directory[_dir]['root'] = run_dir
            directory[_dir]['files'] = []
            # make run subfolders
            if (_dir != 'EM') and (_dir != 'build'):
                for rep in range(1, self.reps+1):
                    rep_path = os.path.join(directory[_dir]['root'], 'rep{}'.format(str(rep)))
                    if not os.path.isdir(rep_path):
                        os.mkdir(rep_path)
                    directory[_dir][rep] = rep_path

        # gather topology & move files to build
        directory['topol'] = []
        for file in os.listdir(directory['build']['root']):
            filepath = os.path.join(directory['build']['root'], file)
            if (not file.endswith('py')) or (not file.endswith('sh')):
                directory['build']['files'].append(filepath)
            if (file.endswith('top')) or (file.endswith('itp')):
                directory['topol'].append(filepath)
        bases = [os.path.basename(fpath) for fpath in directory['build']['files']]
        for file in os.listdir(self.root):
            filepath = os.path.join(self.root, file)
            if (file.endswith('top')) or (file.endswith('itp')): 
                directory['topol'].append(filepath)
            if (not os.path.isdir(filepath)) and (file not in bases):
                shutil.copy(filepath, directory['build']['root'])
                directory['build']['files'].append(filepath)
    

        # make scripts folder
        directory['scripts'] = os.path.join(self.root, 'scripts')
        if not os.path.isdir(directory['scripts']):
            os.mkdir(directory['scripts'])

        # make mdp folder & gather mdp files
        directory['mdps'] = {}
        mdp_dir = os.path.join(self.root, 'mdp')
        if not os.path.isdir(mdp_dir):
            os.mkdir(mdp_dir)
        software_mdps_dir = os.path.join(scriptdir, 'mdp', self.version, self.ff)
        for filename in os.listdir(software_mdps_dir):
            software_mdp_path = os.path.join(software_mdps_dir, filename)
            shutil.copy(software_mdp_path, mdp_dir)
            if filename == 'em.mdp':
                directory['mdps']['em'] = os.path.join(mdp_dir, 'em.mdp')
            if filename == 'nvt.mdp':
                directory['mdps']['nvt'] = os.path.join(mdp_dir, 'nvt.mdp')
            if filename == 'npt.mdp':
                directory['mdps']['npt'] = os.path.join(mdp_dir, 'npt.mdp')
            if filename == 'md.mdp':
                directory['mdps']['md'] = os.path.join(mdp_dir, 'md.mdp')
            if filename == 'ions.mdp':
                directory['mdps']['ions'] = os.path.join(mdp_dir, 'ions.mdp')

        # # fix force field for infer
        # if (self.cluster == 'infer') and (self.ff == 'charmm36-mar2019'):
        #     charmmpath = os.path.join(scriptdir, 'ff', 'charmm36-mar2019.ff')
        #     build_charmm = os.path.join(directory['build']['root'], 'charmm36-mar2019.ff')
        #     if (not os.path.isdir(build_charmm)):
        #         shutil.copytree(charmmpath, os.path.join(directory['build']['root'], 'charmm36-mar2019.ff'))

        return directory

    def getHeader(self, job_name, walltime, nodes=1, ntasks_per_node=24, gpus=1):
        header = []
        header.append('#!/bin/bash\n')
        header.append(f'#SBATCH --nodes={nodes}\n')
        header.append(f'#SBATCH --ntasks-per-node={ntasks_per_node}\n')
        if (gpus > 0):
            header.append('#SBATCH --gres=gpu:{}\n'.format(gpus))
        if self.cluster == 'infer':
            header.append('#SBATCH -p {}\n'.format(self.partition))
        header.append(f'#SBATCH -t {walltime}\n')
        header.append('#SBATCH -A bevanlab\n')
        header.append('#SBATCH --mail-type=all\n')
        header.append(f'#SBATCH --mail-user={self.email}\n')
        header.append(f'#SBATCH --job-name={job_name}\n\n')
        # header.append('# load modules, export library path, load gromacs\n')
        # header.append('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/bevanlab/software/cascades/fftw/3.3.8/lib:/groups/bevanlab/software/cascades/gromacs/2019.3/lib64 \n')
        header.append('export MODULEPATH=$MODULEPATH:{}\n'.format(self.modulepath))
        if self.cluster == 'infer':
            header.append('module load gromacs-{}/{}\n\n'.format(self.partition.split('_')[0], self.version))
        else:
            header.append('module load gromacs/{}\n\n'.format(self.version))
        # for module in self.modules:
        #     if not module.startswith('source'):
        #         line = f'module load {module}\n'
        #     else:
        #         line = f'{module}\n'
        #     header.append(line)
        return header

    def getNumberChains(self, structure):
        f = open(structure, 'r')
        chain_ids = []
        for line in f:
            if not line.startswith('ATOM'):
                continue
            line_parts = line.split()
            chain_id = line_parts[4]
            if not (chain_id.isalpha()):
                raise ValueError('No Chain ID')
            if not (chain_id in chain_ids):
                chain_ids.append(chain_id)
        f.close()
        return len(chain_ids)

    def run_cmd(self, cmd):
        if not cmd.startswith('echo'):
            process = subprocess.Popen(cmd.strip().split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif cmd.startswith('grep'):
            process = subprocess.Popen(cmd.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            process = subprocess.Popen(cmd.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('********************')
        print(cmd.strip())
        if cmd.startswith(self.mdrunPrefix):
            print('Waiting...')
            process.wait()
        output, error = process.communicate()
        print(output.decode())
        print(error.decode())
        return output, error
    
    def build(self, structure, pname='NA', nname='CL', nterm=0, cterm = 0, run=True):
        cwd = os.getcwd()
        os.chdir(self.directory['build']['root'])
        filename = os.path.abspath(os.path.join(self.directory['scripts'], 'build.sh'))
        f = open(filename, 'w')
        header = self.getHeader(job_name='{}_build'.format(self.name), walltime='1:00:00', gpus=0)
        for line in header:
            f.write(line)
        f.write('\n\n')
        _structure = os.path.join(self.root, structure)
        if _structure not in os.listdir(self.directory['build']['root']):
            shutil.copy(os.path.abspath(_structure),self.directory['build']['root'])
        base = self.directory['build']['root']
        cmd = f'cd {self.buildpath}\n'
        f.write(cmd)
        # clean_pdb = os.path.join(base, 'clean.pdb')
        # cmd = 'grep -v HOH {} > {}\n'.format(_structure, clean_pdb)
        # f.write(cmd)
        # if (self.interactive) and (run):
        #     print('Cleaning pdb...')
        #     self.run_cmd(cmd)
        #     print('Cleaning done')
        processed = os.path.join(base, 'processed.gro')
        if self.ff == 'charmm36-mar2019':
            water = 'tip3p'
        else:
            water = 'spc'
        option = self.ffSelectOption
        topol = os.path.join(base, 'topol.top')
        if self.ff == 'charmm36-mar2019':
            ff = 'charmm36-mar2019'
        chains = self.getNumberChains(_structure)
        if (cterm > 0) or (nterm > 0):
            cmd = f'echo {nterm} {cterm} | gmx pdb2gmx -ff {ff} -f {_structure} -o {processed} -ignh -water {water} -ter\n'
            if chains > 1:
                frag = cmd.split('|')[1]
                frag = ' {} {} |'.format(nterm, cterm) + frag
                for k in range(1, chains):
                    frag = ' {} {}'.format(nterm, cterm) + frag
                cmd = 'echo' + frag
        else:
            cmd = f'gmx pdb2gmx -ff {ff} -f {_structure} -o {processed} -ignh -water {water}\n'
        f.write(cmd)
        if (self.interactive) and (run):
            print('doing pdb2gmx...')
            self.run_cmd(cmd)
            print('pdb2gmx done')
        newbox = os.path.join(base, 'box.gro')
        cmd = 'gmx editconf -f {} -o {} -c -d 1.0 -bt cubic\n'.format(processed, newbox)
        f.write(cmd)
        if (self.interactive) and (run):
            print('doing editconf...')
            self.run_cmd(cmd)
            print('editconf done')
        solv = os.path.join(base, 'solv.gro')
        cmd = 'gmx solvate -cp {} -cs spc216.gro -o {} -p {}\n'.format(newbox, solv, topol)
        f.write(cmd)
        if (self.interactive) and (run):
            print('solvating...')
            self.run_cmd(cmd)
            print('solvation done')
        ions_tpr = os.path.join(base, 'ions.tpr')
        cmd = 'gmx grompp -f {} -c {} -p {} -o {}\n'.format(self.directory['mdps']['ions'], solv, topol, ions_tpr)
        f.write(cmd)
        if (self.interactive) and (run):
            print('running grompp for ions...')
            self.run_cmd(cmd)
            print('grompp done')
        solv_ions = os.path.join(base,'ions.gro')
        cmd = 'echo 13 | gmx genion -s {} -o {} -p {} -pname {} -nname {} -conc 0.150 -neutral\n'.format(ions_tpr, solv_ions, topol, pname, nname)
        f.write(cmd)
        if (self.interactive) and (run):
            print('doing genion...')
            self.run_cmd(cmd)
            print('running genion done')
        f.close()
        if run is True:
            if not self.interactive:
                os.chdir(self.directory['scripts'])
                process = subprocess.Popen('sbatch {}'.format(filename).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                output, error = process.communicate()
                print(output.decode())
                print(error.decode())
        os.chdir(cwd)

    def em(self, gpu=1, maxwarn=0, run=True):
        cwd = os.getcwd()
        # os.chdir(self.directory['scripts'])
        # get header
        header = self.getHeader(job_name='{}_em'.format(self.name), walltime='10:00:00', ntasks_per_node=12, gpus=gpu)
        filename = os.path.abspath(os.path.join(self.directory['scripts'], 'em.sh'))
        f = open(filename, 'w')
        for line in header:
            f.write(line)
        f.write('\n\n')

        cmd = 'cd {}\n\n'.format(self.directory['EM']['root'])
        f.write(cmd)
        os.chdir(self.directory['EM']['root'])
        
        solv_ions = os.path.join(self.buildpath,'ions.gro')
        topol = os.path.join(self.buildpath, 'topol.top')
        self.directory['EM']['tpr'] = os.path.join(self.directory['EM']['root'], 'em.tpr')
        cmd = 'gmx grompp -f {} -c {} -p {} -o {} -maxwarn {}'.format(self.directory['mdps']['em'], solv_ions, topol, 'em.tpr', maxwarn)
        f.write(cmd)
        if (self.interactive) and (run):
            print('Running grompp...')
            self.run_cmd(cmd)
            print('grompp done')
        f.write('\n\n')
        if gpu > 0:
            prefix = self.mdrunPrefix
            cmd = '{} -s {} -mp {} -v -deffnm em -gpu_id 0'.format(prefix, self.directory['EM']['tpr'], topol)
        else:
            prefix = 'gmx mdrun'
            cmd = '{} -s {} -mp {} -v -deffnm em'.format(prefix, self.directory['EM']['tpr'], topol)
        f.write(cmd)
        if (self.interactive) and (run):
            print('Running EM...')
            self.run_cmd(cmd)
            print('EM complete')
        f.write('\n\n')
        f.write('wait\n\n')
        f.write('\n\n')
        cmd = 'cd {}\n'.format(self.directory['scripts'])
        f.close()
        if run is True:
            if not self.interactive:
                process = subprocess.Popen('sbatch {}'.format(filename).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                output, error = process.communicate()
                print(output.decode())
                print(error.decode())
        os.chdir(cwd)
        
    def nvt(self, run=True, nodes=1, ntasks_per_node=12, gpu=1, maxwarn=0):
        cwd = os.getcwd()
        os.chdir(self.directory['scripts'])
        scripts = []
        for rep in range(1, self.reps+1):
            filename = os.path.join(self.directory['scripts'], 'nvt{}.sh'.format(rep))
            gpu_id = 0
            f = open(filename, 'w')
            header = self.getHeader(job_name='{}{}_NVT'.format(self.name,rep), walltime='20:00:00', nodes=nodes, ntasks_per_node=ntasks_per_node, gpus=gpu)
            for line in header:
                f.write(line)
            f.write('\n\n')
            cmd = 'cd ../NVT/{}\n'.format(rep)
            os.chdir(os.path.join(self.directory['NVT']['root'], 'rep{}'.format(rep)))
            f.write(cmd)
            em_gro = os.path.join(self.directory['EM']['root'], 'em.gro')
            # nvt_tpr = os.path.join(self.directory['NVT']['root'], str(rep), 'nvt.tpr')
            topol = os.path.join(self.buildpath, 'topol.top')
            nvt_tpr = 'nvt.tpr'
            cmd = 'gmx grompp -f {} -c {} -r {} -p {} -o {} -maxwarn {}\n'.format(self.directory['mdps']['nvt'], em_gro, em_gro, topol, nvt_tpr, maxwarn)
            f.write(cmd)
            if (self.interactive) and (run):
                print('Running grompp for NVT, rep{}...'.format(rep))
                self.run_cmd(cmd)
                print('Grompp done')
            cmd = '{} -gpu_id {} -nt {} -deffnm nvt\n\n'.format(self.mdrunPrefix, str(gpu_id), ntasks_per_node)
            f.write(cmd)
            if (self.interactive) and (run):
                print('Running NVT, rep{}...'.format(rep))
                self.run_cmd(cmd)
                print('NVT for rep{} done'.format(rep))
            f.write('wait\n\n')
            cmd = 'cd {}\n'.format(self.directory['scripts'])
            f.write(cmd)
            f.close()
            scripts.append(os.path.basename(filename))
        if run is True:
            if not self.interactive:
                for script in scripts:
                    process = subprocess.Popen('sbatch {}'.format(script).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    process.wait()
                    output, error = process.communicate()
                    print(output.decode())
                    print(error.decode())
        return filename

    def npt(self, nodes=1, ntasks_per_node=12, gpu=1, maxwarn=0, run=True):
        cwd = os.getcwd()
        os.chdir(self.directory['scripts'])
        scripts = []
        for rep in range(1, self.reps+1):
            header = self.getHeader(job_name='{}{}_NPT'.format(self.name, rep), walltime='20:00:00', nodes=nodes, ntasks_per_node=ntasks_per_node,gpus=gpu)
            filename = os.path.join(self.directory['scripts'], 'npt{}.sh'.format(rep))
            f = open(filename, 'w')
            for line in header:
                f.write(line)
            f.write('\n\n')
            cmd = 'cd ../NPT/{}/\n'.format(str(rep))
            f.write(cmd)
            os.chdir(os.path.join(self.directory['NPT']['root'], f'rep{str(rep)}'))
            nvt_gro = os.path.join(self.directory['NVT']['root'], f'rep{str(rep)}', 'nvt.gro')
            npt_tpr = 'npt.tpr'
            topol = os.path.join(self.buildpath, 'topol.top')
            cmd = 'gmx grompp -f {} -c {} -r {} -p {} -o {} -maxwarn {}\n'.format(self.directory['mdps']['npt'], nvt_gro, nvt_gro, topol, npt_tpr, maxwarn)
            f.write(cmd)
            if (self.interactive) and (run):
                print('Running grompp for NPT, rep{}...'.format(rep))
                self.run_cmd(cmd)
                print('Grompp done')
            cmd = '{} -gpu_id {} -nt {} -deffnm npt\n\n'.format(self.mdrunPrefix, str(0), ntasks_per_node)
            f.write(cmd)
            if (self.interactive) and (run):
                print('Running NPT, rep{}...'.format(rep))
                self.run_cmd(cmd)
                print('NPT for rep{} done'.format(rep))
            cmd = 'cd {}\n'.format(self.directory['scripts'])
            f.write(cmd)
            f.close()
            scripts.append(os.path.basename(filename))
        if run is True:
            for script in scripts:
                if not self.interactive:
                    process = subprocess.Popen('chmod u+x {}'.format(script).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    process = subprocess.Popen('./{}'.format(script).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    process.wait()
                    output, error = process.communicate()
                    print(output.decode())
                    print(error.decode())
        return filename
    
    def production(self, nodes=1, ntasks_per_node=12, gpu=1, run=False, start=0, length=1000, interval=100, walltime='45:00:00', auto_submit=False):
        cwd = os.getcwd()
        scripts = []
        mdrun_directories = []
        os.chdir(self.directory['scripts'])
        # edit MDP to run for the specified interval
        steps = self.nsToSteps(ns=interval)
        self.editMDP(steps)
        sim_start = start
        for rep in range(1, self.reps+1):
            runs = []
            part = 1
            for i in range(sim_start,length, interval):
                start = i
                stop = i + interval
                steps = self.nsToSteps(ns=stop)
                filename = os.path.join(self.directory['MDrun']['root'], 'rep{}'.format(rep), 'production_{}_{}.sh'.format(start, stop))
                f = open(filename, 'w')
                header = self.getHeader(job_name='{}{}_{}_{}'.format(self.name, rep, start, stop),nodes=nodes, ntasks_per_node=ntasks_per_node, walltime=walltime, gpus=gpu)
                for line in header:
                    f.write(line)
                f.write('\n\n')
                mdp = self.directory['mdps']['md']
                topol = os.path.join(self.buildpath, 'topol.top')
                md_tpr = 'md.{}.{}.tpr'.format(start, stop)
                deffnm = 'md.{}.{}'.format(start, stop)
                current_interval = '{}_{}ns'.format(start, stop)
                cmd = f'mkdir {current_interval}\n'
                f.write(cmd)
                cmd = f'cd {current_interval}\n'
                f.write(cmd)
                if i == sim_start:
                    cr = os.path.join(self.directory['NPT']['root'], f'rep{str(rep)}', 'npt.gro')
                    cpt = os.path.join(self.directory['NPT']['root'], f'rep{str(rep)}', 'npt.cpt')
                    cmd = f'gmx grompp -f {mdp} -c {cr} -r {cr} -t {cpt} -p {topol} -o {md_tpr} \n'
                    f.write(cmd)
                    f.write('\n')
                    cmd = f'{self.mdrunPrefix} -deffnm {deffnm} -nt {ntasks_per_node} -gpu_id 0 \n'
                    f.write(cmd)
                else:
                    cmd = f'gmx convert-tpr -s {last_tpr} -o {md_tpr} -nsteps {steps}\n'
                    f.write(cmd)
                    cmd = f'{self.mdrunPrefix} -deffnm {deffnm} -nt {ntasks_per_node} -gpu_id 0 -cpi {last_cpi} -noappend\n'
                    f.write(cmd)
                f.write('\nwait\n\n')
                runs.append(current_interval)
                if auto_submit:
                    if (stop + interval) <= length:
                        cmd = 'cd {}\n'.format(self.directory['scripts'])
                        f.write(cmd)
                        cmd = 'sbatch {}\n'.format(filename)
                        f.write(cmd)
                        f.write('\nexit;\n')
                f.close()
                last_tpr = os.path.join(current_interval, md_tpr)
                if part > 1:
                    cpi_basename = '{}.cpt'.format(deffnm)
                    last_cpi = os.path.join(current_interval, cpi_basename)
                else:
                    last_cpi = os.path.join(current_interval, '{}.cpt'.format(deffnm))
                part += 1
                if i == sim_start:
                    scripts.append(filename)
                    mdrun_directories.append(os.path.join(self.directory['MDrun']['root'], 'rep{}'.format(rep)))

        if run is True:
            for script, mdrun_directory in zip(scripts, mdrun_directories):
                os.chdir(mdrun_directory)
                process = subprocess.Popen('sbatch {}'.format(script).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                output, error = process.communicate()
                print(output.decode())
                print(error.decode())
        os.chdir(cwd)


    def editMDP(self, steps):
        f = open(self.directory['mdps']['md'], 'r')
        contents = f.readlines()
        f.close()
        newcontents = []
        for line in contents:
            if not line.startswith('nsteps'):
                newcontents.append(line)
                continue
            line_parts = line.split('=')
            ns = self.nsToSteps(steps=steps)
            newline = ''.join([line_parts[0], '= ', str(steps), '    ; ({} ns)\n'.format(ns)])
            newcontents.append(newline)
        f = open(self.directory['mdps']['md'], 'w')
        for line in newcontents:
            f.write(line)
        f.close()

    def nsToSteps(self, ns=None, steps=None):
        if ns is not None:
            return ns * 500000
        if steps is not None:
            return steps / 500000
        if (ns is None) and (steps is None):
            return None
    
    @staticmethod 
    def minmax(_min, _max, val):
        _realmin = _min
        _realmax = _max
        if val < _min:
            _realmin = _min
        if val > _max:
            _realmax = _max
        return _realmin, _realmax
    
    @staticmethod
    def calcBoxSize(structure, unit='angstrom', decimals=5, edge_distance=10):
        protein = Protein(structure)
        coords = []
        for atom in protein.atoms:
            coords.append(atom.coordinates)
        df = pd.DataFrame(coords, columns=['x', 'y', 'z'])
        print(df)
        dx = (df['x'].max() + edge_distance) - (df['x'].min() - edge_distance) 
        dy = (df['y'].max() + edge_distance) - (df['y'].min() - edge_distance)
        dz = (df['z'].max() + edge_distance) - (df['z'].min() - edge_distance) 
        if unit == 'angstrom':
            dx = dx / 10
            dy = dy / 10
            dz = dz / 10
        return round(dx, decimals), round(dy, decimals), round(dz, decimals)
        



            

pdb = 'D:/Work/Projects/amyloidNSF/membrane/ab_control/decamer/deca_v1.pdb'
print(Setup.calcBoxSize(pdb, unit='nm'))
# ast = Setup('/work/cascades/kelsieking23/fetub/asticin', reps=4, name='ast')
# ast.em('model1.pdb')
# wt = Setup('/work/cascades/kelsieking23/fetub/wt', reps=4, name='wt', home='/home/kelsieking23/software')
# wt.production(run=False)
# mut = Setup('/work/cascades/kelsieking23/fetub/mut', reps=4, name='mut', home='/home/kelsieking23/software')
# mut.nvt(run=False)

# setup = Setup(cluster='cascades', root='/mnt/d/Work/amyloidNSF/md/bend', reps=4, home='/home/kelsieking23', email='kelsieking23@vt.edu', name='be')
# setup.editMDP(steps=50000000)
# setup.build(structure='bend_monomer.pdb')
# setup.em(gpu=0)
# setup = Setup(cluster='cascades', root='/mnt/d/Work/per2/md/itasser2021_per2_wt', reps=4, name='it2021', home='/home/kelsieking23/software')
# setup.em(structure='prep2_model3.pdb', run=False)