import os
import json
import pandas as pd
from datetime import datetime
import mdtraj
from pymd.mdanalysis.postprocess import PostProcess

class Analysis:

    def __init__(self, inp, top, parent=None, **kwargs):
        '''
        For now, takes a dict of parameters from run.py call
        '''
        self.parent = parent
        self._inp = inp
        self._topfile = top
        self.traj = None
        self._traj = None
        self._output = None
        self.job_name = 'analysis'
        self.job_params = {}
        self.__dict__.update(kwargs)

    
    @property
    def inp(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self._inp)  # type: ignore
        return self._inp

    @property
    def topfile(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self._topfile)  # type: ignore
        return self._topfile

    def save(self):
        params = {}
        manual_keys = ['parent', 'df', 'matrix', 'traj', '_traj']
        for key, value in self.__dict__.items():
            if key not in manual_keys:
                params[key] = value
        filename = os.path.join(self.root, 'job_params.json')
        with open(filename, 'w') as f:
            params_dict = json.dumps(params)
            f.write(params_dict)
        self.job_params = params_dict

    @classmethod
    def from_json(cls, path, inp=None, top=None, parent=None):
        with open(os.path.join(path, 'job_params.json'), 'r') as f:
            params = json.load(f)
        dic = {}
        manual_keys = ['inp', 'top', 'parent', 'df', 'matrix']
        for key, value in params.items():
            if key not in manual_keys:
                dic[key] = value
        for filename in os.listdir(path):
            if filename.endswith('csv'):
                dic['df'] = pd.read_csv(os.path.join(path, filename), index_col=0)
            if filename.endswith('xvg'):
                post = PostProcess()
                dic['df'] = post.getDataFrame(os.path.join(path, filename))
        if parent is not None:
            inp = parent.inp
            top = parent.topfile
        return cls(inp, top, parent, **dic)

    
    @property
    def output(self):
        if self.parent is not None:
            if not os.path.isdir(os.path.join(self.parent.root, self.job_name)): # type: ignore
                os.mkdir(os.path.join(self.parent.root, self.job_name)) # type: ignore
            return os.path.join(self.parent.root, self.job_name, self._output) # type: ignore
        return self._output
    
    @property
    def root(self):
        if self.parent is not None:
            return os.path.join(self.parent.root, self.job_name)
        else:
            return os.path.join(os.getcwd(), self.job_name)

    def loadTrajectory(self, stride=100, selection='backbone', b=0, e=-1):
        if stride != 0:
            traj = mdtraj.load(self.inp, top=self.topfile, stride=stride)
        else:
            traj = mdtraj.load(self.inp, top=self.topfile)
        traj = traj.superpose(traj)

        self._traj = traj
        sele = traj.top.select(selection)
        traj = traj.atom_slice(sele)
        if (e == -1):
            self.traj = traj.center_coordinates()[b:]
        else:
            self.traj = traj.center_coordinates()[b:e]
        self.stride = stride
        self.selection=selection
        self.b = b
        self.e = e
        return self

    @staticmethod
    def now():
        return datetime.now().strftime("%d/%m/%Y %H:%M:%S")