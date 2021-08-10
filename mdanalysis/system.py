import os
import pandas as pd
import numpy as np

if os.name == 'nt':
    import sys
    sys.path.append('D:/Work/iapp')
    from pymd.plot.plotter import Plotter, SystemPlot
    from pymd.mdanalysis.postprocess import PostProcess
    from pymd.structure.protein import Protein
    from pymd.mdanalysis.run import Run
    from pymd.mdanalysis.analyze import Analyzer
else:
    import sys
    sys.path.append('/work/cascades/kelsieking23/iapp_analysis/scripts/python')
    from plot.plotter import Plotter, SystemPlot
    from mdanalysis.postprocess import PostProcess
    from structure.protein import Protein
    from mdanalysis.run import Run
    from mdanalysis.analyze import Analyzer

class System:
    '''
    System performs data collection operations on completed MD systems.
    Primarily writes & executes bash scripts for obtaining data such as RMSD, RMSF, etc.
    '''

    def __init__(self, root=None, reps=None, gro=None, name=None, alias=None, source=None, peptides=None, ligands=None, cascades=True, email=None):
        try:
            self.root = os.path.abspath(root)
        except:
            self.root = None
        self.reps = reps
        self.name = name
        self.source = source
        self.peptides = peptides
        self.ligands = ligands
        self.cascades = cascades
        self.email = email
        self.alias = alias
        if cascades == True:
            self.directory = self.populateDirectory()
            self.__dict__.update(self.directory)
        else:
            self.directory = None
        if (self.root is not None) and (gro is None):
            # try:
            from random import randrange
            if self.reps > 1:
                randrep = randrange(1, self.reps)
                gro = self.directory[randrep]['gro']
                self.gro = gro
                self.protein = Protein(gro, ignore=['ACE', 'NH2'], ligands=self.ligands)
            else:
                if gro is not None:
                    self.gro = gro 
                    self.protein = Protein(structure=gro, ligands=ligands, ignore=['ACE', 'NH2'])
                else:
                    self.gro = None
                    self.protein = None
            # except:
            #     print('except')
            #     if gro is not None:
            #         self.gro = None
            #         self.protein = None
            #     else:
            #         self.gro = gro
            #         self.protein = Protein(structure=gro, ligands=ligands, ignore=['ACE', 'NH2'])
        else:
            if gro is not None:
                self.gro = gro
                self.protein = Protein(structure=gro, ligands=ligands, ignore=['ACE', 'NH2'])
            else:
                self.gro = None
                self.protein = None
        self.post = PostProcess(self)
        self.run = Run(self, cascades)
        self.plot = SystemPlot(self)
        self.analyze = Analyzer(self)

    def populateDirectory(self):
        '''
        Populates self.directory
        '''
        directory = {}

        # find .tpr file & xtc file; add to directory

        if not os.path.isdir(os.path.join(self.root, 'scripts')):
            os.mkdir(os.path.join(self.root, 'scripts'))
        directory['scripts'] = os.path.join(self.root, 'scripts')
        

        # images 

        if not os.path.isdir(os.path.join(self.root, 'images')):
            os.mkdir(os.path.join(self.root, 'images'))
        directory['images'] = {}
        directory['images']['root'] = os.path.join(self.root, 'images')
        
        analysis_types = ['dssp', 'rmsd', 'rmsf', 'mindist', 'hbonds', 'gyration']
        for atype in analysis_types:
            if not os.path.isdir(os.path.join(directory['images']['root'], atype)):
                os.mkdir(os.path.join(directory['images']['root'], atype))
            directory['images'][atype] = os.path.join(directory['images']['root'], atype)

        for rep in range(1, self.reps+1):
            _dir = os.path.join(self.root, str(rep))
            if not os.path.isdir(_dir):
                _dir = os.path.join(self.root, 'rep{}'.format(rep))
            directory[rep] = {}
            directory[rep]['root'] = _dir

            # tpr & gro
            tprs = []
            gros = []
            xtcs = []
            for roots, dirs, files in os.walk(_dir):
                for _file in files:
                    if 'tpr' in _file.split('.')[-1]:
                        tpr_path = os.path.join(roots, _file)
                        tprs.append(tpr_path)
                    if 'gro' in _file.split('.')[-1]:
                        gro_path = os.path.join(roots, _file)
                        gros.append(gro_path)
                    if 'xtc' in _file.split('.')[-1]:
                        xtc_path = os.path.join(roots, _file)
                        xtcs.append(xtc_path)
                    
            _tpr = None
            highest_ns = 0
            for tpr in tprs:
                path_parts = tpr.split(os.sep)
                base = path_parts[-1].split('.')[0]
                try:
                    last_ns = int(base.split('_')[-1])
                except:
                    continue
                if last_ns > highest_ns:
                    highest_ns = last_ns
                    _tpr = tpr
                else:
                    continue
            directory[rep]['tpr'] = _tpr


            _gro = None
            highest_ns = 0
            for gro in gros:
                path_parts = gro.split(os.sep)
                base = path_parts[-1].split('.')[0]
                try:
                    last_ns = int(base.split('_')[-1])
                except:
                    continue
                if last_ns > highest_ns:
                    highest_ns = last_ns
                    _gro = gro
                else:
                    continue
            directory[rep]['gro'] = _gro
            # xtc
            if 'cat_pbc.xtc' in os.listdir(_dir):
                directory[rep]['xtc_system'] = os.path.join(_dir, 'cat_pbc.xtc')
            else:
                directory[rep]['xtc_system'] = None  
            if 'cat.xtc' in os.listdir(_dir):
                directory[rep]['xtc_nopbc'] = os.path.join(_dir, 'cat.xtc')
            else:
                directory[rep]['xtc_nopbc'] = None
            if 'cat_pbc_pro.xtc' in os.listdir(_dir):
                directory[rep]['xtc_pro'] = os.path.join(_dir, 'cat_pbc_pro.xtc') 
            else:
                directory[rep]['xtc_pro'] = None
            if 'cat_pbc_pro_sm.xtc' in os.listdir(_dir):
                directory[rep]['xtc_pro_sm'] = os.path.join(_dir, 'cat_pbc_pro_sm.xtc')
            else:
                directory[rep]['xtc_pro_sm'] = None
            
            # individual run .xtc files
            directory[rep]['run_xtcs'] = [] 
            for xtc in xtcs:
                path_parts = xtc.split(os.sep)
                base = path_parts[-1].split('.')[0]
                try:
                    base_split = base.split('_')
                    if (base_split[1].isnumeric()) and (base_split[2].isnumeric()):
                        if xtc not in directory[rep]['run_xtcs']:
                            directory[rep]['run_xtcs'].append(xtc)
                except:
                    continue

            temp = []       
            nums = []
            prev = None
            for xtc in directory[rep]['run_xtcs']:
                path_parts = xtc.split(os.sep)
                base = path_parts[-1].split('.')[0]
                base_split = base.split('_')
                end = int(base_split[-1])
                nums.append(end)
            zipped = list(zip(nums, directory[rep]['run_xtcs']))
            for item in zipped:
                end = item[0]
                if prev is None:
                    temp.append(end)
                    prev = True
                else:
                    start = len(temp)
                    for i in range(len(temp)):
                        if end > temp[i]:
                            continue
                        if end < temp[i]:
                            if end not in temp:
                                temp.insert(i, end)
                    if len(temp) == start:
                        if end not in temp:
                            temp.append(end)
            sort = []
            for num in temp:
                for item in zipped:
                    if num == item[0]:
                        sort.append(item[1])

            # deal with duplicates, if any
            for xtc in sort:
                _break = False
                path_parts = xtc.split(os.sep)
                base = path_parts[-1].split('.')[0]
                base_split = base.split('_')
                end = int(base_split[-1])
                if (nums.count(end) > 1):
                    for item in zipped:
                        if end == item[0]:
                            for _xtc in sort:
                                if item[1] == _xtc:
                                    sort.remove(_xtc)
                                    nums.remove(end)
                                    _break = True
                                    break
                        if _break == True:
                            break
            directory[rep]['run_xtcs'] = sort

            # get mindist data, if any
            mindist_path = os.path.join(directory[rep]['root'], 'mindist')
            if os.path.isdir(mindist_path):
                if len(os.listdir(mindist_path)) == 0:
                    directory[rep]['mindist'] = None
                else:
                    directory[rep]['mindist'] = {}
                    directory[rep]['mindist']['root'] = mindist_path
                    directory[rep]['mindist']['data'] = {}
                    for folder in os.listdir(mindist_path):
                        directory[rep]['mindist'][folder] = os.path.join(mindist_path, folder)
                        if 'distances.csv' in os.listdir(directory[rep]['mindist'][folder]):
                            directory[rep]['mindist']['data'][folder] = os.path.join(mindist_path, folder, 'distances.csv')
                        else:
                            directory[rep]['mindist']['data'][folder] = None
            else:
                directory[rep]['mindist'] = None 


            directory = self.setupDirBasic(directory, rep, 'rmsd')
            directory = self.setupDirBasic(directory, rep, 'rmsf')
            directory = self.setupDirBasic(directory, rep, 'dssp')
            directory = self.setupDirBasic(directory, rep, 'clusters')
            directory = self.setupDirInteractions(directory, rep, 'mindist')
            directory = self.setupDirInteractions(directory, rep, 'hbonds')
            directory = self.setupDirOther(directory, rep, 'gyration', mkdir=False)

        return directory

    def setupDirBasic(self, directory, rep, location, mkdir=True):
        path = os.path.join(directory[rep]['root'], location)
        if not os.path.isdir(path):
            if mkdir is True:
                os.mkdir(path)
        directory[rep][location] = {}
        directory[rep][location]['root'] = path
        directory[rep][location]['data'] = []
        for filename in os.listdir(path):
            p = os.path.join(directory[rep][location]['root'], filename)
            directory[rep][location]['data'].append(p)
        return directory

    def setupDirInteractions(self, directory, rep, location):
        path = os.path.join(directory[rep]['root'], location)
        if location == 'mindist':
            dtypes = ['sidechain', 'mainchain', 'residue', 'peptide', 'sidechain_sm', 'mainchain_sm', 'residue_sm', 'sm_sm', 'protein_sm']
        if location == 'hbonds':
            dtypes = ['mainchain_sm', 'sidechain_sm', 'residue_sm', 'mainchain_pro', 'sidechain_pro', 'residue', 'sm', 'protein_sm', 'backbone_pro', 'backbone_sm']
        if not os.path.isdir(path):
            os.mkdir(path)
        directory[rep][location] = {}
        for dtype in dtypes:
            subpath = os.path.join(path, dtype)
            if not os.path.isdir(subpath):
                os.mkdir(subpath)
            directory[rep][location][dtype] = {}
            directory[rep][location][dtype]['root'] = subpath
            data = []
            for filename in os.listdir(subpath):
                p = os.path.join(directory[rep][location][dtype]['root'], filename)
                data.append(p)
            directory[rep][location][dtype]['data'] = data
        return directory
        

    def setupDirOther(self,directory,rep,location,mkdir=True):
        path = os.path.join(directory[rep]['root'], location)
        if not os.path.isdir(path):
            if mkdir is True:
                os.mkdir(path)
            else:
                return directory
        else:
            directory[rep][location] = {}
            directory[rep][location]['root'] = path
            for folder in os.listdir(path):
                folder_path = os.path.join(path, folder)
                directory[rep][location][folder] = {}
                directory[rep][location][folder]['root'] = folder_path
                directory[rep][location][folder]['data'] = []
                for filename in os.listdir(folder_path):
                    p = os.path.join(folder_path, filename)
                    directory[rep][location][folder]['data'].append(p)
            return directory
    
          

        
        

                



# mor = System('/work/cascades/kelsieking23/iapp_analysis/HIAPP_Trimer_MOR/MDrun/', 6, name='MOR', peptides=3)
# mor.rmsfSidechain()
# mor.rmsdPerPeptideBackbone()
# qur = System('/work/cascades/kelsieking23/iapp_analysis/HIAPP_Trimer_QUR/MDrun/', 6, name='QUR', peptides=3)
# qur.rmsfSidechain()
# qur.rmsdPerPeptideBackbone()
# myr = System('/work/cascades/kelsieking23/iapp_analysis/HIAPP_Trimer_MYR/MDrun/', 6, name='MYR', peptides=3, ligands='MYR')
# myr.plotMindist('sidechains')
# # myr.mindist(groups=('residue', 'residue'), start=400000, stop=600000)
# myr.mindist(groups=('sidechain', 'sidechain'), start=400000, stop=600000)
# myr.plotRMSF()
# myr.rmsdPerPeptideBackbone()
# epi = System('/work/cascades/kelsieking23/iapp_analysis/HIAPP_Trimer_EPI/', 6, name='EPI', peptides=3)
# epi.mindist(groups=('sidechain', 'sidechain'), start=400000, stop=600000)
# epi.rmsdPerPeptideBackbone()
# ctrl = System('/work/cascades/kelsieking23/iapp_analysis/HIAPP_Trimer/MDrun/', 6, name='CTRL', peptides=3)
# ctrl.rmsdPerPeptideBackbone()
# ctrl.rmsfSidechain()
# ctrl.cluster(stop=600000, sm=False)
# mut = System('/work/cascades/kelsieking23/iapp_analysis/HIAPP_Trimer_Leu23/MDrun/', 6, name='MUT', peptides=3)
# mut.getCatPBC(run=True)
# from pathlib import Path
# folder = Path('../tests/')
# gro = folder / 'md_550_600.gro'
# test = System(None, None, 'test', None, 3, cascades=False, ligands='MYR')
# test.plotMindist('md_500_600.gro')