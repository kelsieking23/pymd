import os
import sys

if sys.platform == 'win32':
    sys.path.append('D:/Work')
    import pymd
    from pymd.mdsetup.mdsetup_dir import scriptdir
if sys.platform == 'linux':
    # sys.path.append('/work/cascades/kelsieking23/iapp_analysis/scripts/python')
    # from mdsetup.mdsetup_dir import scriptdir
    sys.path.append('/mnt/d/Work')
    from pymd.mdsetup.mdsetup_dir import scriptdir

class DirectoryBuilder:

    def __init__(self, setup):
        self.wizard = setup
        directory = {}

    def setupRunFolders(self):
        '''
        setup run folders. internal function used by Setup.populateDirectory()
        '''
        dirs = ['EM', 'build', 'NVT', 'NPT', 'MDrun']
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
            # make run subfolders
            if (_dir != 'EM') or (_dir != 'build'):
                for rep in range(1, self.reps+1):
                    rep_path = os.path.join(directory[_dir]['root'], str(rep))
                    if not os.path.isdir(rep_path):
                        os.mkdir(rep_path)
                    directory[_dir][rep] = rep_path
        # make scripts folder
        directory['scripts'] = os.path.join(self.root, 'scripts')
        if not os.path.isdir(directory['scripts']):
            os.mkdir(directory['scripts'])

        # make mdp folder
        directory['mdps'] = {}
        mdp_dir = os.path.join(self.root, 'mdp')
        if not os.path.isdir(mdp_dir):
            os.mkdir(mdp_dir)
    
    def setupTopology(self):
        '''
        gather topology files from build folder
        '''
        if self.wizard.top is not None:
            directory['topol'] = [os.path.join(self.wizard.top, file) for file in os.listdir(self.wizard.top)]
        else:
            directory['topol'] = []
