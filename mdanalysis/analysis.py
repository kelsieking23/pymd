import os
import json
import pandas as pd
import numpy as np
from datetime import datetime
import mdtraj
from pymd.mdanalysis.postprocess import PostProcess
import multiprocessing as mp
from pymd.utilities.rewritepdb import writePDB
from collections.abc import Iterable
from pymd.utilities.library import lipids, residues

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
        self.traj_iter = None
        self._output = None
        self.job_name = 'analysis'
        self.job_params = {}
        self.load_state = False
        self.verbose = True
        self._top = mdtraj.load(self.topfile).topology
        self._root = os.path.dirname(self._inp)
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

    @property
    def top(self):
        if self.traj is None:
            return self._top
        else:
            return self.traj.topology

    @property
    def chain_idx(self):
        chain_idx = []
        for chain in self.top.chains:
            idx = chain.residue(0).chain.index
            if idx not in chain_idx:
                chain_idx.append(idx)
        return np.array(chain_idx)
    
    @property
    def protein_idx(self):
        return self.top.select('protein')
    
    @property
    def protein_atoms(self):
        return np.array([self.top.atom(idx) for idx in self.protein_idx])
    
    def save(self, **kwargs):
        params = {}
        manual_keys = ['parent', 'df', 'matrix', 'traj', '_traj', 'top', 'frames']
        to_dump = {}
        for key, value in self.__dict__.items():
            try:
                json.dumps(value)
                to_dump[key] = value
            except:
                continue
        for k, v in kwargs.items():
            to_dump[k] = v
        filename = os.path.join(self.root, 'job_params.json')
        with open(filename, 'w') as f:
            params_dict = json.dumps(to_dump)
            f.write(params_dict)
        self.job_params = params_dict

    @classmethod
    def from_json(cls, path, inp=None, top=None, parent=None, load_traj=False):
        with open(os.path.join(path, 'job_params.json'), 'r') as f:
            params = json.load(f)
        dic = {}
        manual_keys = ['inp', 'top', 'parent', 'df', 'matrix']
        for key, value in params.items():
            if key not in manual_keys:
                dic[key] = value
        for filename in os.listdir(path):
            if filename == params['_output']:
                if filename.endswith('csv'):
                    dic['df'] = pd.read_csv(os.path.join(path, filename), index_col=0)
                if filename.endswith('xvg'):
                    post = PostProcess()
                    dic['df'] = post.getDataFrame(os.path.join(path, filename))
        if parent is not None:
            inp = parent.inp
            top = parent.topfile
        else:
            inp = params['_inp']
            top = params['_topfile']
        if load_traj:
            stride = dic['stride']
            selection = dic['selection']
            b = dic['b']
            e = dic['e']
            if stride != 0:
                traj = mdtraj.load(inp, top=top, stride=stride)
            else:
                traj = mdtraj.load(inp, top=top)
            _traj = traj
            traj = traj.superpose(traj)

            if selection != 'all':
                sele = traj.top.select(selection)
                traj = traj.atom_slice(sele)
            if (e == -1):
                traj = traj.center_coordinates()[b:]
                # self.traj = traj[b:]
            else:
                traj = traj.center_coordinates()[b:e]
            dic['traj'] = traj
            dic['_traj'] = _traj
            dic['frames'] = traj._xyz
        return cls(inp, top, parent, **dic)

    
    @property
    def output(self):
        if self._output is not None:
            return os.path.join(self.root, self._output) # type: ignore
        return os.path.join(self.root, 'analysis.csv')
        # if self.parent is not None:
        #     if not os.path.isdir(os.path.join(self.parent.root, self.job_name)): # type: ignore
        #         os.mkdir(os.path.join(self.parent.root, self.job_name)) # type: ignore
        #     return os.path.join(self.parent.root, self.job_name, self._output) # type: ignore
        # return self._output
    
    @property
    def root(self):
        if self.parent is not None:
            root = os.path.join(self.parent.root, self.job_name)
        else:
            if os.path.dirname(self.inp) == '':
                root =  os.path.join(os.getcwd(), self.job_name)
            else:
                root = os.path.join(os.path.dirname(self.inp), self.job_name)
        if not os.path.isdir(root):
            os.mkdir(root)
        return root

    def loadTrajectory(self, stride=100, selection='backbone', b=0, e=-1):
        if self.verbose:
            print('Loading {}...'.format(self.inp))
        if stride != 0:
            self._traj = mdtraj.load(self.inp, top=self.topfile, stride=stride)
        else:
            self._traj = mdtraj.load(self.inp, top=self.topfile)

        if self.verbose:
            print(f'Slicing selection: "{selection}" ...')
        if selection != 'all':
            if isinstance(selection, str):
                if (not selection == ''):
                    sele = self._traj.top.select(selection)
                    self.traj = self._traj.atom_slice(sele)
                else:
                    self.traj = self._traj
            elif isinstance(selection, (list, tuple, pd.DataFrame)):
                self.traj = self._traj.atom_slice(np.array(selection))
            elif isinstance(selection, np.ndarray):
                self.traj = self._traj.atom_slice(selection)
            elif (selection is None):
                self.traj = self._traj
            else:
                raise ValueError('Selection must be string, list, tuple, np.ndarray, or None')
        else:
            self.traj = self._traj

        if (e == -1):
            self.traj = self.traj[b:]
            # self.traj = traj[b:]
        else:
            if self.verbose:
                print(f'Slicing interval: ({b}, {e})')
            self.traj = self._traj[b:e]
            # self.traj = traj[b:e]
        self.stride = stride
        self.selection=selection
        self.b = b
        self.e = e
        self.frames = self.traj._xyz
        if self.verbose:
            print('Trajectory loaded.')
            print(self.traj)
            print(f'Trajectory shape {self.traj._xyz.shape}')
        # self.top = self.traj.topology
        return self

    def iterloadTrajectory(self, stride=1, selection='all', chunk=0):
        if self.verbose:
            print('Loading trajectory iterator...')
        if isinstance(selection, str):
            if (selection == '') or (selection == 'all'):
                sele = None
            else:
                sele = self.top.select(selection)
        elif isinstance(selection, (list, tuple)):
            sele = np.array(selection)
        elif isinstance(selection, np.ndarray):
            sele = selection
        elif selection is None:
            sele = None
        else:
            raise ValueError('Selection must be string, list, tuple, np.ndarray, or None')
        if sele is None:
            self.traj_iter = mdtraj.iterload(self.inp, top=self.top, chunk=chunk, stride=stride)
        else:
            self.traj_iter = mdtraj.iterload(self.inp, top=self.top, chunk=chunk, stride=stride, atom_indices=sele)
        if self.verbose:
            print('Trajectory loaded.')
            print(self.traj_iter)
        return self.traj_iter
    def superpose(self):
        if self.traj is not None:
            self.traj = self.traj.superpose(self.traj)
            return self.traj
        else:
            print('Trajectory not loaded')
            return None
    
    def center(self):
        if self.traj is not None:
            self.traj = self.traj.center_coordinates()
            return self.traj
        else:
            print('Trajectory not loaded')
            return None

    def select(self, selection):
        return self.top.select(selection)
    
    @staticmethod
    def now():
        return datetime.now().strftime("%d/%m/%Y %H:%M:%S")


    def getPartitions(self, nprocs='auto'):
        if nprocs == 'auto':
            nprocs = int(mp.cpu_count() // 2)
        else:
            nprocs = nprocs
        if self.traj is None:
            self.loadTrajectory()
        nframes, _, _ = self.traj._xyz.shape # type: ignore
        interval = int(nframes // nprocs)
        partitions = []
        procid=1
        for i in range(0, nframes, interval):
            data = {
                'b':i,
                'e':i+interval,
                'procid':procid,
            }
            partitions.append(data)
            procid+=1
            if ((i + interval + interval) > nframes) and (i+interval != nframes):
                data = {
                    'b':i+interval,
                    'e':nframes,
                    'procid':procid
                }
                partitions.append(data)
                break
        return partitions, nprocs

    def toPDB(self, index, output, full_traj=False, renumber=False, remark=None):
        if full_traj:
            frame = self._traj._xyz[index]
        else:
            frame = self.traj._xyz[index]
        chain_index = 0
        chain_id = 'A'
        contents = []
        if remark is not None:
            contents.append('{}\n'.format(remark))
        for z in range(0, len(frame)):
            if full_traj:
                atom = self._traj.topology._atoms[z]
            else:
                atom = self.traj.topology._atoms[z]
            if atom.residue.name in residues():
                if atom.residue.chain.index > chain_index:
                    chain_index = atom.residue.chain.index
                    chain_id = chr(ord(chain_id) + 1)
            else:
                chain_id = ' '
            if renumber:
                res_num = str(atom.residue.index + 1)
            else:
                res_num = str(atom.residue.resSeq)
            x, y, z = map(self.fixCoordinates, frame[z])
            line = ['ATOM', str(atom.index), atom.name, atom.residue.name, chain_id, res_num, x, y, z, '1.00', '0.00', atom.element.symbol, atom.residue.segment_id]
            contents.append(line)
        writePDB(contents, output)
        print('Wrote {}'.format(output))

    def fixCoordinates(self, xyz):
        return xyz*10
    
    @staticmethod
    def chain_conversions():
        keys = [i for i in range(0,26)]
        values = list(map(chr, range(ord('A'), ord('Z')+1)))
        return {k:v for (k,v) in zip(keys, values)}