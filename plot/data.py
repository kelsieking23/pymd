import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

if os.name == 'nt':
    sys.path.append('D:/Work/iapp/')
    from pymd.mdanalysis.postprocess import PostProcess
else:
    sys.path.append('/work/cascades/kelsieking23/iapp_analysis/scripts/python')
    from mdanalysis.postprocess import PostProcess


class PlotData:

    def __init__(self, data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto):
        self.data = data
        self.fig = fig
        self.xticks = xticks
        self.yticks = yticks
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        self.axes = axes
        self.legend = legend
        self.annotations = annotations
        self.saveto = saveto

    '''
    Single-System Class Methods
    '''

    @classmethod
    def dsspOverTime(cls, df, title=None, annotations=None, output=None):
        colors = {
            'coil_percent':'#000080',
            'bsheet_percent':'#a00000',
            'helix_percent':'#000000'
        }
        if 'coil_std' in df.columns:
            coil = Data(df=df, x=df.index, y=df['coil_percent'], color=colors['coil_percent'], label='Coil', fill_between=df['coil_std'])
            bsheet = Data(df=df, x=df.index, y=df['bsheet_percent'], color=colors['bsheet_percent'], label=r'$\beta$-Strand', fill_between=df['bsheet_std'])
            helix = Data(df=df, x=df.index, y=df['helix_percent'], color=colors['helix_percent'], label='Helix', fill_between=df['helix_std'])
        else:
            coil = Data(df=df, x=df.index, y=df['coil_percent'], color=colors['coil_percent'], label='Coil')
            bsheet = Data(df=df, x=df.index, y=df['bsheet_percent'], color=colors['bsheet_percent'], label=r'$\beta$-Strand')
            helix = Data(df=df, x=df.index, y=df['helix_percent'], color=colors['helix_percent'], label='Helix')
        
        fig = None

        data = [coil, bsheet, helix]

        xticks = ElementParam(xmin=0, xmax=df.index[-1], locs=500, fontsize=20, minor_locs=100)
        yticks = ElementParam(ymin=0, ymax=100, locs=20, fontsize=20, minor_locs=10)

        xlabel = ElementParam(label='Time (ns)', fontsize=22)
        ylabel = ElementParam(label='Percentage (%)', fontsize=22)

        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=True)

        legend = ElementParam(loc='upper right', fontsize=14)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def dsspPerResidue(cls, dfs, colors, title=None, xlabels=None, output=None):
        data = []
        s = 500
        i = 0
        for df in dfs:
            d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label=df.name, s=s)
            data.append(d)
            s -= 100
            i += 1
        fig = None
        
        if xlabels is not None:
            if xlabels[0] != 0:
                xlabels.insert(0,0)
        xmin = (3*df.index[0] - df.index[1])/2
        xmax = (3*df.index[-1] - df.index[-2])/2
        xticks = ElementParam(locs=1, fontsize=17, xlabels=xlabels, xmin=xmin, xmax=xmax)
        
        yticks = ElementParam(fontsize=17)

        xlabel = ElementParam(label='Residue', fontsize=20)
        ylabel = ElementParam(label='Probability', fontsize=20)
        
        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=True)

        legend = ElementParam(loc='upper right', markerscale=0.60, ncol=3, individual=False)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def gyration(cls, df, title=None, color=None, label=None, output=None):
        data = Data(df=df, x=df['time'], y=df['gyr'], color=color, label=label)
        data = [data]
        fig = None

        xticks = ElementParam(xmin=0, xmax=list(df['time'])[-1], locs=100, fontsize=18, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=10, locs=2, minor_locs=0.5, fontsize=18)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label='Radius of Gyration (nm)', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes=ElementParam(off=False, semiopen=True)

        legend = ElementParam(loc='upper right')

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def rmsdPerPeptide(cls, df, colors=None, title=None, output=None):
        colors = {
            'peptide1':'#11174b',
            'peptide2':'#3967a4',
            'peptide3':'#62bed2',
        }

        data = []
        i = 0
        for column in df.columns:
            if isinstance(colors, dict):
                d = Data(df=df, x=df.index, y=df[column], color=colors[column], label='Peptide {}'.format(i+1))
            if isinstance(colors, list):
                d = Data(df=df, x=df.index, y=df[column], color=colors[i], label='Peptide {}'.format(i+1))
            data.append(d)
            i += 1
        
        fig = None

        xticks = ElementParam(xmin=0, xmax=df.index[-1], locs=100, fontsize=18, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=1.0, locs=0.2, fontsize=18, minor_locs=0.1)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label='RMSD (nm)', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes=ElementParam(off=False, semiopen=True)

        annotations = None

        legend = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)
    
    @classmethod
    def mindistOverTime(cls, df, colors=None, title=None, legend=None, output=None):
        
        data = []
        i = 0
        _max = None
        for column in df.columns:
            column_parts = column.split('_')
            label = 'Peptide {} and Peptide {}'.format(column_parts[0][1], column_parts[1][1])
            if isinstance(colors, dict):
                d = Data(df=df, x=df.index, y=df[column], color=colors[column], label=label)
            if isinstance(colors, list):
                d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=label)
            if _max is None:
                _max = df[column].max()
            else:
                if df[column].max() > _max:
                    _max = df[column].max()
            data.append(d)
            i += 1

        fig = ElementParam(width=10, height=6)

        xticks = ElementParam(xmin=0, xmax=df.index[-1], fontsize=18, locs=100, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=1500, minor_locs=200, fontsize=18)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label='Number of Contacts', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes = axes = ElementParam(off=False, semiopen=True)

        if legend is True:
            legend = ElementParam(loc='upper right')
        else:
            legend = None
        
        annotations = [ElementParam(atype='plot', x=df.index, y=[_max]*len(df.index), linestyle='--', color='black', linewidth=1),
                        ElementParam(atype='annotate', text=str(round(_max)), xy=(df.index[-1]-40, _max+10), fontsize=18)]

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def hbondsOverTime(cls, dfs, colors=None, title=None, xlabels=None, ymax=2, output=None, block_average=100):
        data = []
        i = 0
        _ymax = None
        s = len(dfs.keys())*block_average
        for t in range(block_average, (len(dfs.keys())*block_average)+1, block_average):
            df = dfs[str(t)]
            start_time = t - block_average
            label = '{}-{} ns'.format(int(start_time), int(t))
            d = Data(df=df, x=df.index, y=df['hbonds'], yerr=df['stdev'], color=colors[i], label=label, s=s)
            i += 1
            data.append(d)
            if _ymax is None:
                _ymax = df['hbonds'].max()
            else:
                if df['hbonds'].max() > ymax:
                    _ymax = df['hbonds'].max()
            s -= block_average

        fig = None

        if ymax is None:
            ymax = round(_ymax) + 2
        
        if xlabels is not None:
            if xlabels[0] != 0:
                xlabels.insert(0,0)
        xticks = ElementParam(locs=1, fontsize=17, xlabels=xlabels)
        yticks = ElementParam(ymin=0, ymax=ymax, fontsize=17, locs=0.5, minor_locs=0.25)

        xlabel = ElementParam(label='Residue', fontsize=20)
        ylabel = ElementParam(label='Number of Hydrogen Bonds', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=False)

        legend = ElementParam(loc='upper right', markerscale=0.60)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def hbondsHeatmap(cls, df, colormap=None, title=None, name=None, vmin=0, vmax=2, output=None):
        data = Data(df=df, vmin=vmin, vmax=vmax, colormap=colormap)
        xticks = ElementParam(locs=1, xlabels=[])
        yticks = ElementParam(locs=1, ylabels=[])
        fig = None
        title = ElementParam(title=title, fontsize=24)
        legend = ElementParam(colorbar=True)
        xlabel = None
        ylabel = None

        axes = None
        annotations = None

        saveto = output
        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def mindist(cls, df, xlabels=None, ylabels=None, colormap='plasma', legend=True, title=None, vmin=0, vmax=1, panel=True, output=None):
        data = Data(df=df, vmin=vmin, vmax=vmax, colormap=colormap)
        if xlabels is None:
            xlabels = []
        else:
            if xlabels[0] != 0:
                xlabs = [0]
                for item in xlabels:
                    xlabs.append(item)
                xlabels = xlabs
        if ylabels is not None:
            # ylabels = list(reversed(ylabels))
            if ylabels[0] != 0:
                ylabs = [0]
                for item in ylabels:
                    ylabs.append(item)
                ylabels = ylabs
        else:
            ylabels = []
        if panel is True:
            offset=0.35
        else:
            offset=0.25
        xticks = ElementParam(locs=1, xlabels=xlabels, offset=offset, fontsize=16)
        yticks = ElementParam(locs=1, ylabels=ylabels, offset=offset, fontsize=16)
        fig = None
        title = ElementParam(title=title, fontsize=24)
        if legend is True:
            legend = ElementParam(colorbar=True)
        else:
            legend = None
        xlabel = ElementParam(label='Residue', fontsize=20)
        ylabel = ElementParam(label='Residue', fontsize=20)

        axes = None
        anno_data = df.round(2).to_numpy()
        annotations = ElementParam(data=anno_data, fontsize=16)

        saveto = output
        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def mindistResidueSM(cls, df, colormap=None, title=None, xlabels=None, ylabels=None, offset=False, rotation=None, colorbar=True, vmin=0, vmax=1, name=None, output=None):
        # if name != 'EPI':
        #     df = df.reindex(['ringC_carbonyl_norm', 'ringC_hydroxyl_norm', 'ringB_hydroxyl_norm', 'ringA_hydroxyl_norm', 'ringC_norm', 'ringB_norm', 'ringA_norm', '{}_norm'.format(name)])
        # else:
        #     df = df.reindex(['ringC_hydroxyl_norm', 'ringB_hydroxyl_norm', 'ringA_hydroxyl_norm', 'ringC_norm', 'ringB_norm', 'ringA_norm', '{}_norm'.format(name)])
        data = Data(df=df, vmin=vmin, vmax=vmax, colormap=colormap)
        if xlabels is None:
            xlabels = []
        else:
            if xlabels[0] != 0:
                xlabs = [0]
                for item in xlabels:
                    xlabs.append(item)
                xlabels = xlabs
        if ylabels is not None:
            ylabels = list(reversed(ylabels))
            if ylabels[0] != 0:
                ylabs = [0]
                for item in ylabels:
                    ylabs.append(item)
                ylabels = ylabs
        else:
            ylabels = []
        if (offset is True) and (xlabels != []):
            xticks = ElementParam(locs=1, xlabels=xlabels, offset=0.5, rotation=rotation, fontsize=20)
            yticks = ElementParam(locs=1, ylabels=ylabels, offset=0.5, fontsize=20)
        else:
            xticks = ElementParam(locs=1, xlabels=xlabels)
            yticks = ElementParam(locs=1, ylabels=ylabels)

            

        fig = ElementParam(width=16, height=12)

        title = ElementParam(title=title, fontsize=36)

        legend = ElementParam(colorbar=colorbar)
        
        xlabel = ElementParam()
        ylabel = ElementParam()

        axes = None
        
        anno_data = df.round(2).to_numpy()
        annotations = ElementParam(data=anno_data, fontsize=20)

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def rmsfBarPlot(cls, df, title=None, color=None, labels=None, output=None):
        if labels is None:
            data = Data(df=df, color=color, label=None, width=0.5, xpos=[i for i in range(len(df.index))], capsize=4)
        else:
            data = Data(df=df, color=color, label=None, width=0.5, xpos=[i for i in range(len(df.index))], capsize=4)

        if labels is None:
            xticks = ElementParam(labels=list(df.index), fontsize=18, locs=[i for i in range(len(df.index))])
        else:
           xticks = ElementParam(labels=labels, fontsize=18, locs=[i for i in range(len(df.index))])
    
        yticks = ElementParam(ymin=0, ymax=4, locs=0.5, minor_locs=0.25, fontsize=18)   

        fig = None

        xlabel = ElementParam(label=None, fontsize=24)
        ylabel = ElementParam(label='RMSF (nm)', fontsize=18)

        title = ElementParam(title=title, fontsize=20)

        legend = None

        axes = ElementParam(off=False, semiopen=True)

        annotations = [ElementParam(atype='autolabel', fontsize=14)]

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)


    '''
    Multi-System Class Methods
    '''

    @classmethod
    def gyrationMulti(cls, df, title=None, ecc=False, colors=None, output=None):
        
        data = []
        i = 0
        for column in df.columns:
            if isinstance(colors, dict):
                d = Data(df=df, x=df.index, y=df[column], color=colors[column], label=column)
            if isinstance(colors, list):
                d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=column)
            data.append(d)
            i += 1
        
        fig = None

        xticks = ElementParam(xmin=0, xmax=df.index[-1], locs=100, fontsize=18, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=10, locs=2, fontsize=18, minor_locs=0.5)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        
        if ecc is False:
            ylabel = ElementParam(label='Radius of Gyration (nm)', fontsize=20)
        else:
            ylabel = ElementParam(label='Eccentricity', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes=ElementParam(off=False, semiopen=True)

        legend = ElementParam(loc='upper right')

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def hbondsTotalBarMulti(cls, df, title=None, color=None, labels=None, output=None):
        if labels is None:
            data = Data(df=df, color=color, label=None, xpos=[i for i in range(len(df.index))], capsize=4)
        else:
            data = Data(df=df, color=color, label=None, xpos=[i for i in range(len(df.index))], capsize=4)
        
        if labels is None:
            xticks = ElementParam(labels=list(df.index), fontsize=18, locs=[i for i in range(len(df.index))])
        else:
           xticks = ElementParam(labels=labels, fontsize=18)
        
        fig = None

        yticks = ElementParam(ymin=0, ymax=3, locs=0.5, minor_locs=0.25, fontsize=18)

        xlabel = ElementParam(label=None, fontsize=24)
        ylabel = ElementParam(label='Average Hydrogen Bonds', fontsize=18)

        title = ElementParam(title=title, fontsize=20)

        axes = ElementParam(off=False, semiopen=True)

        legend = None

        annotations = [ElementParam(atype='autolabel', fontsize=14)]

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def contactDistributionMulti(cls, dfs, title=None, colors=None, labels=None, output=None):
        data = []
        i = 0
        width = 0.4
        for df in dfs:
            if i == 0:
                xpos = np.arange(len(df.index)) - (width/2)
            if i == 1:
                xpos = np.arange(len(df.index)) + (width/2)
            if labels is None:
                d = Data(df=df, color=colors[i], xpos=xpos, width=width, label=df.name, capsize=4)
            else:
                d = Data(df=df, color=colors[i], xpos=xpos, width=width, label=df.name, capsize=4)
            data.append(d)
            i += 1
        

        if labels is None:
            xticks = ElementParam(labels=list(df.index), fontsize=18, locs=[i for i in range(len(df.index))])
        else:
           xticks = ElementParam(labels=labels, fontsize=18, locs=[i for i in range(len(df.index))])

        fig = None

        yticks = ElementParam(ymin=0, ymax=8, locs=1, minor_locs=0.25, fontsize=18)

        xlabel = ElementParam(label=None, fontsize=24)
        ylabel = ElementParam(label='Interaction Frequency', fontsize=18)

        title = ElementParam(title=title, fontsize=20)

        axes = ElementParam(off=False, semiopen=True)

        legend = ElementParam(loc='upper right')

        annotations = [ElementParam(atype='autolabel', fontsize=14)]

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def mindistOverTimeMulti(cls, dfs, title=None, colors=None, legend=None, labels=None, output=None):
        data = []
        i = 0
        for df in dfs:
            d = Data(df=df, x=df.index, y=df['mean'], color=colors[i], label=labels[i], fill_between=df['stdev'])
            data.append(d)
            i += 1

        fig = ElementParam(width=10, height=6)

        xticks = ElementParam(xmin=0, xmax=df.index[-1], fontsize=18, locs=100, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=5, locs=1, minor_locs=0.25, fontsize=18)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label='Normalized Number of Contacts', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes = axes = ElementParam(off=False, semiopen=True)

        if legend is True:
            legend = ElementParam(loc='upper right')
        else:
            legend = None
        
        # annotations = [ElementParam(atype='plot', x=df.index, y=[_max]*len(df.index), linestyle='--', color='black', linewidth=1),
        #                 ElementParam(atype='annotate', text=str(round(_max)), xy=(df.index[-1]-40, _max+10), fontsize=18)]
        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def dsspOverTimeMulti(cls, df, title=None, colors=None, legend=None, labels=None, output=None):
        data = []
        i = 0
        for column in df.columns:
            print(column)
            print(labels)
            d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=labels[i])
            i += 1
            data.append(d)
        # coil = Data(df=df, x=df.index, y=df['coil_percent'], color=colors['coil_percent'], label='Coil')
        # bsheet = Data(df=df, x=df.index, y=df['bsheet_percent'], color=colors['bsheet_percent'], label=r'$\beta$-Strand')
        # helix = Data(df=df, x=df.index, y=df['helix_percent'], color=colors['helix_percent'], label='Helix')
        
        fig = None


        xticks = ElementParam(xmin=0, xmax=df.index[-1], locs=100, fontsize=18, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=100, locs=20, fontsize=18, minor_locs=10)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label='Percentage (%)', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=True)

        legend = ElementParam(loc='upper right')

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)
class Data:

    def __init__(self, df, **kwargs):
        self.df = df
        self.__dict__.update(kwargs)
        

class ElementParam:

    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    
    


