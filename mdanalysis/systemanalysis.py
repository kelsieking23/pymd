import os
import pandas as pd

from pymd.mdanalysis.postprocess import PostProcess
from pymd.plot.data import PlotData
from pymd.plot.plotter import Plotter


class SystemAnalysis:

    def __init__(self, parent, **kwargs):
        self.parent = parent
        self.name = None
        self.__dict__.update(**kwargs) # this is very general right now, idk what it will look like
        self.process = PostProcess(self.parent)
        self.plotter = Plotter()

    def plot(self, average_by=None, ptype='infer', output=None, show=True, **kwargs):
        # average_by should be 
        #   None (meaning dont average, make a panel plot
        #   All (meaning average each dataframe in its entirety)
        #   column/1 (meaning average by column) (are these the same thing)
        #   index/0 (average by index)
        if ptype == 'infer':
            tmp = self.parent.rep1.getChildByJobName(self.name)
            if hasattr(tmp, 'df'):
                ptype = tmp.df.attrs['ptype']
        else:
            raise ValueError('could not detect ptype from data. please speticy ptype in function call')
        if ptype == 'timeseries':
            self.plotTimeseries(average_by=average_by, output=output, show=show, **kwargs)
        elif ptype == 'heatmap':
            pdata = PlotData.heatmap(self.df, output=output, **kwargs)
            self.plotter.heatmap(pdata, show=show)
        else:
            raise ValueError('no valid type')

    def plotTimeseries(self, average_by=None, output=None, show=True, **kwargs):
        if average_by is None:
            # plot panel
            pdatas = []
            title = kwargs['title']
            for rep in self.parent._reps:
                analysis = rep.getChildByJobName(self.name)
                pdatas.append(PlotData.timeseries(analysis.df, **kwargs))
            ax = self.plotter.timeseriesPanel(pdatas, title=title, output=output, show=show)
        elif average_by == 'all':
            # plot average overall
            dfs = []
            for rep in self.parent._reps:
                analysis = rep.getChildByJobName(self.name)
                dfs.append(analysis.df)
            df = pd.concat(dfs).groupby(level=0).mean()
            pdata = PlotData.timeseries(df, output=output, **kwargs)
            ax = self.plotter.timeseries(pdata, show=show)
        elif (average_by == 0) or (average_by == 'index'):
            dfs = []
            for rep in self.parent._reps:
                analysis = rep.getChildByJobName(self.name)
                dfs.append(analysis.df)
            df = pd.concat(dfs, axis=1).mean(axis=1)
            pdata = PlotData.timeseries(df, output=output, **kwargs)
            ax = self.plotter.timeseries(pdata, show=show)
        elif (average_by == 1) or (average_by == 'column'):
            # (i think this is the same thing as all but i literally do not know idk idk help)
            dfs = []
            for rep in self.parent._reps:
                analysis = rep.getChildByJobName(self.name)
                dfs.append(analysis.df)
            df = pd.concat(dfs).groupby(level=0).mean()
            pdata = PlotData.timeseries(df, output=output, **kwargs)
            ax = self.plotter.timeseries(pdata, show=show)
        else:
            raise ValueError('average_by must be none, all, index, columns')
        return ax

    def averageAll(self):
        dfs = []
        for rep in self.parent._reps:
            analysis = rep.getChildByJobName(self.name)
            dfs.append(analysis.df)
        attrs = analysis.df.attrs
        df = pd.concat(dfs).groupby(level=0).mean()
        df.attrs = attrs
        return df
    def averageIndex(self):
        df = pd.DataFrame()
        i = 0
        for rep in self.parent._reps:
            analysis = rep.getChildByJobName(self.name)
            for col in analysis.df.columns:
                df[i] = analysis.df[col]
                i += 1
        df['mean'] = df.mean(axis=1)
        df = df.drop([col for col in df.columns if col != 'mean'], axis=1)
        return df

    def plot_with(self, systems, nrows, ncols, average_by='all', ptype='infer', sharex=True, sharey=True, suptitle=None, titles=[], output=None, show=True, **kwargs):
        if ptype == 'infer':
            tmp = self.parent.rep1.getChildByJobName(self.name)
            if hasattr(tmp, 'df'):
                if not hasattr(tmp.df, 'attrs'):
                    ptype='infer'
                else:
                    ptype = tmp.df.attrs['ptype']
        if ptype == 'timeseries':
            # axes.append(self.plotTimeseries(average_by=average_by, output=None, show=False))
            pdatas = []
            if average_by == 'all':
                df = self.averageAll()
            if average_by == 'index':
                df = self.averageIndex()
            if titles != []:
                title = titles[0]
            else:
                title = None
            print(df)
            pdatas.append(PlotData.timeseries(df, output=None, title=title, **kwargs))
            i = 1
            for system in systems:
                if titles != []:
                    title = titles[i]
                else:
                    title = None
                analysis = system.getChildByJobName(self.name)
                df = analysis.averageAll()
                pdatas.append(PlotData.timeseries(df, output=None, title=title, **kwargs))
                i += 1
            self.plotter.timeseriesPanel(pdatas, ncols, nrows, sharex=sharex, sharey=sharey, output=output, title=suptitle, show=show)
        elif ptype == 'heatmap':
            pass
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


        
        