import os
import types
import pandas as pd
import mdtraj

from pymd.mdanalysis.postprocess import PostProcess
from pymd.mdanalysis.cluster import Cluster
from pymd.mdanalysis.dist import Distance
from pymd.mdanalysis.dssp import DSSP
from pymd.plot import PlotData, Plotter

class Subsystem(PostProcess):

    def __init__(self, root, inp, top, id=None, parent=None, tu='ns'):
        self.root = root
        self._inp = inp
        self._topfile = top
        self.job = None
        self.df = pd.DataFrame()
        self.data = pd.DataFrame()
        self.id = id
        self.files = {}
        # if dict is not None:
        #     self.__dict__.update(dict)
        if parent is not None:
            self.parent = parent
        self.traj = None
        self._traj = None
        self._tu = tu
        self.plotter = Plotter()
    
    @classmethod
    def from_dict(cls, dic, parent):
        return cls(dic['root'], dic['xtc'], dic['gro'], dic['id'], parent)

    @property
    def inp(self):
        return self._inp
    
    @property
    def topfile(self):
        return self._topfile

    @property
    def top(self):
        return mdtraj.load(self.topfile).topology

    @property
    def tu(self):
        conversions = {
            1:'ps',
            1000:'ns'
        }
        if (self.traj is not None) and (not isinstance(self.traj, types.GeneratorType)):
            return conversions[self.traj.timestep]
        else:
            return self._tu

    def load(self, job_name, select=-1):
        self.job = job_name
        for file in os.listdir(self.root):
            ext = os.path.splitext(file)[-1][1:]
            base = os.path.splitext(file)[0]
            if ext not in self.files.keys():
                self.files[ext] = []
            self.files[ext].append(os.path.join(self.root, file))
            if (ext == 'csv') or (ext == 'xvg'):
                if base == job_name:
                    self.df = self.getDataFrame(os.path.join(self.root, file))
        if 'cluster' in job_name:
            pass
            # self.cluster = Cluster(self.inp, self.topfile)
            # self.cluster.df = self.df
            # self.cluster.root = os.path.join(self.root, job_name)
            # self.cluster
        if 'dist' in job_name:
            self.dist = Distance.from_json(os.path.join(self.root, job_name), parent=self)
            self.df = self.dist.df
            self.data = self.dist.df
        if 'dssp' in job_name:
            self.dssp = DSSP.from_json(os.path.join(self.root, job_name), parent=self)
            self.df = self.dssp.df
            self.data = self.dssp.df
        return self

    def load_trajectory(self, stride=1, selection='all', b=0, e=-1):
        if stride != 0:
            traj = mdtraj.load(self.inp, top=self.topfile, stride=stride)
        else:
            traj = mdtraj.load(self.inp, top=self.topfile)
            
        traj = traj.superpose(traj)
        self._traj = traj
        if (e == -1):
            self.traj = traj.center_coordinates()[b:]
        else:
            self.traj = traj.center_coordinates()[b:e]
        if selection != 'all':
            sele = self.traj.top.select(selection)
            self.traj = self.traj.atom_slice(sele)
        self.top = self.traj.top
        return self


    def timeAverage(self, interval=200):
        '''
        parameters:
        df (dataframe)
        interval (int, default=200): interval for averages
        '''
        try:
            assert isinstance(interval, int)
        except:
            raise ValueError('Interval must be integer. Input {} cannot be interpreted as integer'.format(interval))
        end = int(self.df.index[-1])
        start = int(self.df.index[0])
        mean = pd.DataFrame()
        sd = pd.DataFrame()
        for i in range(start, end, interval):
            mean['{}-{} {}'.format(i,i+interval, self.tu)] = self.df.loc[i:i+interval, :].mean()
            sd['{}-{} {}-std'.format(i,i+interval, self.tu)] =self.df.loc[i:i+interval, :].std()
        data = pd.concat([mean, sd], axis=1)
        return data

    def plot(self, ptype='infer', df = None, output=None, show=True, **kwargs):
        # this needs to create a plotdata object or something
        if output is not None:
            output = os.path.join(self.root, self.job, output)
        if df is None:
            df = self.df
        if ptype == 'infer':
            if 'ptype' in df.attrs.keys():
                ptype = df.attrs['ptype']
                if ptype is None:
                    if 'type' in df.attrs.keys():
                        if df.attrs['type'] == 'xy':
                            ptype = 'timeseries'
            else:
                if 'type' in df.attrs.keys():
                    if df.attrs['type'] == 'xy':
                        ptype = 'timeseries'
        if ptype == 'timeseries':
            pdata = PlotData.timeseries(df, output=output, **kwargs)
            self.plotter.timeseries(pdata, show=show)
        elif ptype == 'heatmap':
            pdata = PlotData.heatmap(df, output=output, **kwargs)
            self.plotter.heatmap(pdata, show=show)
        elif ptype == 'dssp':
            pdata = PlotData.dsspOverTime(df, output=output, **kwargs)
            self.plotter.timeseries(pdata, show=show)
        else:
            raise ValueError('no valid type')

    def plot_with(self, systems, ptype='infer', output=None, show=True, **kwargs):
        pdatas = []
        if ptype == 'infer':
            if 'ptype' in self.df.attrs.keys():
                ptype = self.df.attrs['ptype']
        if ptype == 'timeseries':
            pdatas.append(PlotData.timeseries(self.df, **kwargs))
            for system in systems:
                pdatas.append(PlotData.timeseries(system.df, **kwargs))
            self.plotter.timeseries(pdatas, show=show)
        elif ptype == 'heatmap':
            pdatas.append(PlotData.heatmap(self.df, **kwargs))
            for system in systems:
                pdatas.append(PlotData.heatmap(system.df, **kwargs))
            self.plotter.heatmap(pdatas, show=show, output=output)
        elif ptype == 'dssp':
            pdata = PlotData.dsspOverTime(df, output=output, **kwargs)
            self.plotter.timeseries(pdata, show=show)
        else:
            raise ValueError('no valid type')