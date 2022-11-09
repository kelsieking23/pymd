import os

from pymd.mdanalysis.postprocess import PostProcess
from pymd.plot.data import PlotData
from pymd.plot.plotter import Plotter
from pymd.mdanalysis.cluster import Cluster
from pymd.mdanalysis.dist import Distance


# my idea is that you can do an analysis like this kinda:
#   system = System(args)
#   system = system.run.rmsd()
# then you have a couple of options from there:
#   system.rmsd.plot() -- this would maybe plot the system average & all replicates in the system on a single figure.
#   system.rep1.rmsd.plot() --- this would plot an individual replicate and output a single plot
# i want this thing to hold the 'rmsd' in those function calls 
class Analysis:

    def __init__(self, name, parent=None, **kwargs):
        self.parent = parent
        self._df = None
        self.name = name
        self.root = None
        self.files = {}
        if self.parent is not None:
            path = os.path.join(self.parent.root, self.name)
            if not os.path.isdir(path):
                    os.mkdir(path)
            self.root = path
            # self.files = [os.path.join(path, f) for f in os.listdir(path)]
        self.process = PostProcess(self.parent)
        self.plotter = Plotter()
        self.load()
        self.__dict__.update(**kwargs) # this is very general right now, idk what it will look like
    
    @property
    def df(self):
        return self._df
    
    @df.setter
    def df(self, d):
        if self._df is not None:
            if (hasattr(self._df, 'attrs')):
                _attrs = self._df.attrs
                self._df = d
                self._df.attrs = _attrs
                return self._df
            else:
                print('pymd.mdanalysis.analysis.Analysis @df.setter')
        else:
            if not (hasattr(d, 'attrs')):
                s = list(zip(range(len(d.columns)), d.columns))
                s = {key:value for (key,value) in s}
                attrs = {
                            'title':None,
                            'x_label':None,
                            'y_label':None, 
                            's':s,
                            'ptype':None,
                            'type':None,
                            's_list':list(s.values())
                        }
                self._df = d
                self._df.attrs = attrs
                return self._df
            else:
                self._df = d
                return self._df



    def plot(self, ptype='infer', output=None, show=True, **kwargs):
        # this needs to create a plotdata object or something
        if ptype == 'infer':
            if 'ptype' in self.df.attrs.keys():
                ptype = self.df.attrs['ptype']
                if ptype is None:
                    if 'type' in self.df.attrs.keys():
                        if self.df.attrs['type'] == 'xy':
                            ptype = 'timeseries'
            else:
                if 'type' in self.df.attrs.keys():
                    if self.df.attrs['type'] == 'xy':
                        ptype = 'timeseries'
        if ptype == 'timeseries':
            pdata = PlotData.timeseries(self.df, output=output, **kwargs)
            self.plotter.timeseries(pdata, show=show)
        elif ptype == 'heatmap':
            pdata = PlotData.heatmap(self.df, output=output, **kwargs)
            self.plotter.heatmap(pdata, show=show)
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
        else:
            raise ValueError('no valid type')

    def load(self, select=-1):
        for file in os.listdir(self.root):
            ext = os.path.splitext(file)[-1][1:]
            base = os.path.splitext(file)[0]
            if ext not in self.files.keys():
                self.files[ext] = []
            self.files[ext].append(os.path.join(self.root, file))
            if (ext == 'csv') or (ext == 'xvg'):
                if base == self.name:
                    self.df = self.process.getDataFrame(os.path.join(self.root, file))
        if 'cluster' in self.name:
            self.cluster = Cluster(self.parent.xtc, self.parent.gro)
            self.cluster.df = self.df
            self.cluster.root = os.path.join(self.parent.root, self.name)
            self.parent.cluster = self.cluster
        if 'dist' in self.name:
            self.dist = Distance.from_json(os.path.join(self.parent.root, self.name), parent=self.parent)
            self.parent.dist = self.dist
        return self
        # extensions = ['csv', 'xvg']
        # print(self.name)
        
        # if self.root is not None:
        #     for ext in extensions:
        #         file = '{}.{}'.format(job_name, ext)
        #         print(file)
        #         if file in os.listdir(self.root):
        #             print('yes')
        #             self.df = self.process.getDataFrame(os.path.join(self.root, file))
        #         if self.name == 'cluster':
        #             self.analyze = ClusterData(parent=self, n_clusters=None)
        #         return self
        #     return self
        # else:
        #     if (isinstance(files, tuple)) or (isinstance(files, list)):
        #         self.df = self.process.getDataFrame(*files, select=select)
        #         self.root = os.path.dirname(os.path.abspath(files[0]))
        #         for file in files:
        #             if file not in self.files:
        #                 self.files.append(file)
        #     else:
        #         self.root = os.path.dirname(os.path.abspath(files))
        #         self.df = self.process.getDataFrame(files, select=select)
        #         if files not in self.files:
        #             self.files.append(files)
        # return self
                
    @classmethod
    def from_xvg(cls, xvg, select=-1, **kwargs):
        proc = PostProcess()
        files = []
        if (isinstance(xvg, tuple)) or (isinstance(xvg, list)):
            df = proc.getDataFrame(*xvg, select=select)
            root = os.path.dirname(os.path.abspath(xvg[0]))
            for file in xvg:
                if file not in files:
                    files.append(file)
        else:
            root = os.path.dirname(os.path.abspath(files))
            df = proc.process.getDataFrame(files, select=select)
            if xvg not in files:
                xvg.append(file)
        return cls(parent=None, df=df, root=root, files=files)
            
        

        


        
        



    
    


