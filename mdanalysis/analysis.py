import os

from pymd.mdanalysis.postprocess import PostProcess
from pymd.plot.data import PlotData
from pymd.plot.plotter import Plotter


# my idea is that you can do an analysis like this kinda:
#   system = System(args)
#   system = system.run.rmsd()
# then you have a couple of options from there:
#   system.rmsd.plot() -- this would maybe plot the system average & all replicates in the system on a single figure.
#   system.rep1.rmsd.plot() --- this would plot an individual replicate and output a single plot
# i want this thing to hold the 'rmsd' in those function calls 
class Analysis:

    def __init__(self, parent, **kwargs):
        self.parent = parent
        self.df = None
        self.name = None
        self.__dict__.update(**kwargs) # this is very general right now, idk what it will look like
        path = os.path.join(self.parent.root, self.name)
        if not os.path.isdir(path):
                os.mkdir(path)
        self.root = path
        self.files = [os.path.join(path, f) for f in os.listdir(path)]
        self.process = PostProcess(self.parent)
        self.plotter = Plotter()

    def plot(self, ptype='infer', output=None, show=True, **kwargs):

        # this needs to create a plotdata object or something
        if ptype == 'infer':
            if 'ptype' in self.df.attrs.keys():
                ptype = self.df.attrs['ptype']
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

    def load(self, job_name):
        extensions = ['csv', 'xvg']
        for ext in extensions:
            file = '{}.{}'.format(job_name, ext)
            if file in os.listdir(self.root):
                self.df = self.process.getDataFrame(os.path.join(self.root, file))
                return self
        return self
    
    def systemLoad(self):
        pass


        
        



    
    


