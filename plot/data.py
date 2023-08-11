import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from cycler import cycler
from collections.abc import Iterable
from pandas.core.algorithms import isin
import math


class Data:

    def __init__(self, df, **kwargs):
        self.df = df
        self.__dict__.update(kwargs)
        

class ElementParam:

    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

class PlotData:

    def __init__(self, plot_type, data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto, lineconnect=None, legend_data=None):
        self.plot_type = plot_type
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
        self.lineconnect = lineconnect
        self.legend_data = legend_data

    def plot(self):
        plotter = Plotter()
        if self.plot_type is None:
            raise ValueError('No plot type specified')
        if self.plot_type == 'timeseries':
            plotter.timeseries(self)

    @staticmethod
    def label(df, column, i, labels='columns'):
        if labels is None:
            return '__nolegend__'
        elif labels == 'columns':
            return column
        else:
            return labels[i]

    @staticmethod
    def color(colors, column, i):
        if isinstance(colors, dict):
            return colors[column]
        elif (isinstance(colors, list)) or (isinstance(colors, np.ndarray)):
            return colors[i]
        else:
            return colors
    
    @staticmethod
    def max(df, column, _max=None):
        if _max is None:
            _max = df[column].max()
        else:
            if df[column].max() > _max:
                _max = df[column].max()
        return _max
    
    @staticmethod
    def min(df, column, _min=None):
        if _min is None:
            _min = df[column].min()
        else:
            if df[column].max() > _min:
                _min = df[column].min()
        return _min


    @staticmethod
    def ymax_hist(x, bins, density, _ymax=None):
        counts, _ = np.histogram(x, bins=bins, density=density)
        if _ymax is None:
            return counts.max()
        if counts.max() > _ymax:
            return counts.max()
        return _ymax

    @staticmethod
    def xbounds_hist(x, bins, density, bounds=(None, None)):
        _, bins = np.histogram(x, bins=bins, density=density)
        if (None in bounds):
            return bins.min(), bins.max()
        if bounds[0] > bins.min():
            _xmin = bins.min()
        else:
            _xmin = bounds[0]
        if bounds[1] < bins.max():
            _xmax = bins.max()
        else:
            _xmax = bounds[1]
        return _xmin, _xmax
    '''
    Single-System Class Methods
    '''

    @classmethod
    def timeseries(cls, df, title=None, labels='columns', x_label=None, y_label=None, major_xticks=None, minor_xticks=None, colors=None, output=None, ymin=0.0, ymax=None, 
    legend=False, ncol=1, linewidth=2, average=False, std=False, major_yticks=None, minor_yticks=None, tick_label_fontsize=14, ax_label_fontsize=18, title_fontsize=20,
    legend_fontsize=12, semiopen=True, superlegend=False, alpha=1, grid=False, weight='regular', xtick_labels=None, xtick_label_rotation='horizontal', title_weight='bold',
    scatter=False, marker='o', s=60, xpad=0, legend_loc='best', linestyle='solid'):
        '''
        Arguments:
        df (pandas DataFrame): dataframe containing all the data you want plotted. output of postprocess.getDataFrame
        title (optional, str): title for the plot
        labels (optional, list): labels for legend 
        x_label (str): x axis label (default to value in .xvg)
        y_label (str): y axis label (default to value in .xvg) 
        major_xticks (int): interval for major x ticks
        minor_xticks (int): interval for minor x ticks
        major_yticks (int): interval for major yticks
        minor_yticks (int): interval for minor yticks
        colors (optional, list or dict): either a list of hex values or a dict mapping labels to colors
        output (optional, string): output filepath. if not specified, will not save image. 
        ymin (optional, float/int): minimum for y axis. default is 0
        ymax (optional, float): set the max for y axis. will automatically try to figure out where ymax should be if None
        legend (bool): if you want a legend. false if no legend, true if you want a legend. default: False
        ncol (int): number of columns for the legend. default 1
        linewidth (int/float): linewidth. default 2
        average (bool): if average of dataframe should be taken. default false
        std (bool): if standard deviation should be shown on plot. default false. average must be true 
        semiopen (bool): Whether to hide the axis lines on the right and top of the plot. default is True (hide)
        superlegend (bool): whether to include legend data for multiple potential subplots. default is false
        '''
        if colors is None:
            # colors = ['#1f464c', '#2a8a2d', '#8a47b0', '#9bcccc']
            colors = plt.cm.tab20b(np.linspace(0, 1,len(df.columns)))
            
        fig = None 
        # get plot data
        data = []
        i = 0
        _ymax = None
        # average and std if true
        if average is True:
            _df = df
            df = pd.DataFrame()
            df['mean'] = _df.mean(axis=1)
            df.attrs = _df.attrs

        for column in df.columns:
            # break for std
            if (std is True):
                if ('std' in df.columns):
                    if column == 'std':
                        break
                elif column.endswith('_std'):
                    continue
                else:
                    pass

            # decide labels for legend
            if labels is None:
                label = '__nolegend__'
            elif labels == 'columns':
                label = column
            else:
                if isinstance(labels, (list, tuple, set)):
                    label = labels[i]
                elif isinstance(labels, dict):
                    label = labels[column]
                else:
                    print("WARNING: Parameter labels must be: None, 'columns', list, tuple, set, or dict.\n")
                    print('Parameter labels ({}) is type {}'.format(labels, type(labels)))
                    print('Label will default to None (nolegend)')
                    label = '__nolegend__'
            # decide color
            if isinstance(colors, dict):
                color = colors[column]

            elif (isinstance(colors, list)) or (isinstance(colors, np.ndarray)):
                color = colors[i]
            else:
                color=colors

            # make Data instances and append to data list
            if std is False:
                d = Data(df=df, x=df.index, y=df[column], color=color, label=label, linewidth=linewidth, alpha=alpha, scatter=scatter, marker=marker, s=s, linestyle=linestyle)
            else:
                if 'std' in df.columns:
                    d = Data(df=df, x=df.index, y=df[column], color=color, label=label, linewidth=linewidth, fill_between=df['std'], alpha=alpha, scatter=scatter, marker=marker, s=s, linestyle=linestyle)
                elif '{}_std'.format(column) in df.columns:
                    d = Data(df=df, x=df.index, y=df[column], color=color, label=label, linewidth=linewidth, fill_between=df['{}_std'.format(column)], alpha=alpha, scatter=scatter, marker=marker, s=s, linestyle=linestyle)
                else:
                    try:
                        assert average is True
                        std_df = pd.DataFrame()
                        std_df['std'] = _df.std(axis=1)
                        d = Data(df=df, x=df.index, y=df[column], color=color, label=label, linewidth=linewidth, fill_between=std_df['std'], alpha=alpha, scatter=scatter, marker=marker, s=s, linestyle=linestyle)
                    except:
                        if i == 0:
                            print('WARNING: average must be True for standard deviation to plot. Plotting dataframe containing the following files individually: \n {}'.format(df.name))
                        d = Data(df=df, x=df.index, y=df[column], color=color, label=label, linewidth=linewidth, alpha=alpha, scatter=scatter, marker=marker, s=s, linestyle=linestyle)
            data.append(d)

            i += 1

            # calculate ymax
            if _ymax is None:
                _ymax = df[column].max()
            else:
                if df[column].max() > _ymax:
                    _ymax = df[column].max()
        # tick labels and locations
        if xtick_labels is not None:
            _xtick_labels = [0]
            for l in xtick_labels:
                _xtick_labels.append(l)
            xtick_labels = _xtick_labels
        xticks = ElementParam(xmin=df.index[0]-xpad, xmax=round(df.index[-1])+xpad, locs=major_xticks, fontsize=tick_label_fontsize, minor_locs=minor_xticks, labels=xtick_labels,
                            xtick_label_rotation=xtick_label_rotation)
        # xticks = ElementParam(xmin=df.index[0], xmax=round(1000), locs=major_xticks, fontsize=14, minor_locs=minor_xticks)

        if ymax is None:
            ymax = math.ceil(_ymax)
        yticks = ElementParam(ymin=ymin, ymax=ymax, locs=major_yticks, fontsize=tick_label_fontsize, minor_locs=minor_yticks)

        if (hasattr(df, 'attrs')) and (df.attrs != {}):
            pass
        else:
            df.attrs = {
                'title':'',
                'x_label':'',
                'y_label':'', 
                's':{}
            }
        # axis labels and locations
        if x_label is not None:
            xlabel = ElementParam(label=x_label, fontsize=ax_label_fontsize, weight=weight)
        else:
            xlabel = ElementParam(label=df.attrs['x_label'], fontsize=ax_label_fontsize, weight=weight)
        if y_label is not None:
            ylabel = ElementParam(label=y_label, fontsize=ax_label_fontsize, weight=weight)
        else:
            ylabel = ElementParam(label=df.attrs['y_label'], fontsize=ax_label_fontsize, weight=weight)

        # title
        if title is not None:
            title = ElementParam(title=title, fontsize=title_fontsize, weight=title_weight)
        else:
            title = None

        # axes
        axes=ElementParam(off=False, semiopen=semiopen, grid=grid)

        # annotations
        annotations = None

        #legend 
        if legend is not False:
            legend = ElementParam(loc=legend_loc, ncol=ncol, fontsize=legend_fontsize)
        
        # panel legend data
        if superlegend is True:
            if (labels is not None):
                if isinstance(labels, str):
                    labels = df.columns
                legend_data = []
                i = 0
                for label in labels:
                    legend_data.append((label, colors[i]))
                    i += 1
        else:
            legend_data = None

        # output
        saveto = output

        return cls('timeseries', data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto, legend_data=legend_data)


    @classmethod
    def histogram(cls, df, title='', x_label=None,y_label=None,  output=None, ymin=0, ymax=None, xmin=None, xmax=None, xpad=0.5, n_bins='auto', density=True, center=False, annotate=False, labels='columns', alpha=1, colors=None, legend=False, ncol=2, ignore=[], major_xticks=None, minor_xticks=None, major_yticks=None, minor_yticks=None, tick_label_fontsize=14, ax_label_fontsize=18, title_fontsize=20,
    legend_fontsize=12, semiopen=True, superlegend=False, grid=False, weight='regular', xtick_labels=None, xtick_label_rotation='horizontal', title_weight='bold', legend_loc='upper right'):
        if colors is None:
            # colors = ['#1f464c', '#2a8a2d', '#8a47b0', '#9bcccc']
            colors = plt.cm.tab20b(np.linspace(0, 1,len(df.columns)))

        if center is True:
            df = df.subtract(df.mean())
        data = []
        i = 0
        _ymax = None
        _xmin = None
        _xmax = None
        for column in df.columns:
            if column in ignore:
                continue

            # decide labels for legend
            label = cls.label(df, column, i, labels)
            # decide color
            # colors = cls.color(colors, column, i)
            
            #  bins
            if (n_bins == 'auto') or (isinstance(n_bins, str)):
                n_bins = math.ceil(1 + (3.322 * math.log(len(df[column]))))
            print(n_bins)
            # k=100
            x = df[column]
            print(colors, colors[i])
            d = Data(df=df, x=x, bins=n_bins, density=density, alpha=alpha, color=colors[i], label=label)
            data.append(d)
            i += 1
            # calculate ymax
            _ymax = cls.ymax_hist(x, n_bins, density)
            _xmin, _xmax = cls.xbounds_hist(x, n_bins, density, (_xmin, _xmax))
        if xmax is None:
            xmax = _xmax
        if xmin is None:
            xmin = _xmin
        if ymax is None:
            ymax = _ymax
        fig = None
        xticks = ElementParam(xmin=xmin-xpad, xmax=xmax+xpad, fonlocs=major_xticks, fontsize=tick_label_fontsize, minor_locs=minor_xticks, labels=xtick_labels,
                            xtick_label_rotation=xtick_label_rotation)
        
        yticks = ElementParam(ymin=ymin, ymax=ymax, locs=major_yticks, fontsize=tick_label_fontsize, minor_locs=minor_yticks)

        x_label = ElementParam(label=x_label, fontsize=ax_label_fontsize, weight=weight)

        if y_label is not None:
            y_label = ElementParam(label=y_label, fontsize=ax_label_fontsize, weight=weight)
        else:
            if density is True:
                y_label = ElementParam(label='Probability Density', fontsize=ax_label_fontsize, weight=weight)
            else:
                y_label = ElementParam(label='Count', fontsize=ax_label_fontsize, weight=weight)
                

        title = ElementParam(title=title, fontsize=title_fontsize, weight=title_weight)

        axes = ElementParam(off=False, semiopen=semiopen, grid=grid)

        if (legend is not False) and (legend is not None):
            legend = ElementParam(loc=legend_loc, ncol=ncol, fontsize=legend_fontsize)


        annotations = None
        if annotate is True:
            _max = df['mean'].max()
            annotations = [ElementParam(atype='plot', x=df.index, y=[_max]*len(df.index), linestyle='--', color='black', linewidth=1),
                         ElementParam(atype='annotate', text=str(round(_max)), xy=(df.index[-1]-40, _max+10), fontsize=18)]
        else:
            annotations = None

        saveto = output
        return cls('histogram', data, fig, xticks, yticks, x_label, y_label, title, axes, legend, annotations, saveto, lineconnect=None)

    @classmethod
    def markerPlot(cls, df, colors=None, title=None, marker='o', x_label='Residue', y_label=None, ymin=-0.2, ymax=None, output=None, major_xticks=None, major_yticks=None, 
                    minor_xticks=None, minor_yticks=None, rotation=0, ncol=1, labels='columns', tick_label_fontsize=14, ax_label_fontsize=18, title_fontsize=20,
                    legend_fontsize=12, semiopen=False):
        data = []
        i = 0
        s = 14900/len(df.index)
        smax = s
        _ymax = None
        markers = [marker] * len(df.columns)
        last_marker = None
        if colors is None:
            # colors = ['#1f464c', '#2a8a2d', '#8a47b0', '#9bcccc']
            colors = plt.cm.RdBu(np.linspace(0, 2.5,len(df.columns)))
        # if major_xticks is None:
        #     major_xticks = len(df.index) / 10
        # if minor_xticks is None:
        #     if major_xticks > 1:
        #         minor_xticks = major_xticks / 2
        for column in df.columns:
            try:
                if column.endswith('std'):
                    break
            except:
                pass
            # decide labels for legend
            if labels is None:
                label = '__nolegend__'
            elif labels == 'columns':
                label = column
            else:
                label = labels[i]

            # decide color
            if isinstance(colors, dict):
                color = colors[column]
            if (isinstance(colors, list)) or (isinstance(colors, np.ndarray)):
                color = colors[i]

            # assign marker size 
            if last_marker is None:
                last_marker = markers[i]
            else:
                if markers[i] != last_marker:
                    s = len(df) * 100
            std_col = column + '-std'

            if std_col in df.columns:
                yerr = df[std_col]
                if _ymax is None:
                    _ymax = df[std_col].max()
                else:
                    if df[std_col].max() > _ymax:
                        _ymax = df[std_col].max()
            else:
                yerr=None
            d = Data(df=df, x=df.index, y=df[column], yerr=yerr, capsize=smax/100, color=colors[i], label=label, 
                    s=s, marker=markers[i])

            data.append(d)
            i += 1 
            s = s * 0.75
        # tick labels and locations
        if (isinstance(df.index[-1], int)) or (isinstance(df.index[-1], float)):
            xticks = ElementParam(xmin=df.index[0], xmax=round(df.index[-1]), locs=major_xticks, fontsize=tick_label_fontsize, minor_locs=minor_xticks, rotation=rotation)
        else:
            xticks = ElementParam(xmin=0.1, xmax=len(df.index)+0.9, locs=major_xticks, fontsize=tick_label_fontsize, minor_locs=minor_xticks, rotation=rotation)
        if ymax is None:
            if _ymax is None:
                ymax = math.ceil(df.max().max())
            else:
                ymax = math.ceil(df.max().max() + _ymax)

        yticks = ElementParam(ymin=ymin, ymax=ymax, locs=major_yticks, fontsize=tick_label_fontsize, minor_locs=minor_yticks)

        if (hasattr(df, 'attrs')) and (df.attrs != {}):
            pass
        else:
            df.attrs = {
                'title':'',
                'x_label':'',
                'y_label':'', 
                's':{}
            }
        # axis labels and locations
        if x_label is not None:
            xlabel = ElementParam(label=x_label, fontsize=ax_label_fontsize)
        else:
            xlabel = ElementParam(label=df.attrs['x_label'], fontsize=ax_label_fontsize)
        if y_label is not None:
            ylabel = ElementParam(label=y_label, fontsize=ax_label_fontsize)
        else:
            ylabel = ElementParam(label=df.attrs['y_label'], fontsize=ax_label_fontsize)


        fig = ElementParam(width=8, height=6, layout='tight')


        xlabel = ElementParam(label=x_label, fontsize=16)
        ylabel = ElementParam(label=y_label, fontsize=16)

        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=semiopen)

        legend = ElementParam(loc='upper right', markerscale=0.60, ncol=ncol, fontsize=16)

        annotations = None

        saveto = output

        return cls('marker', data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)
    @classmethod
    def dsspOverTime(cls, df, title=None, structure='all', annotations=None, output=None, legend=True, colors=None, alpha=1, std=False, linewidth=2, gridline=None):
        if colors is None:
            colors = plt.cm.tab20b(np.linspace(0, 1,len(df.columns)))
            colors = {
                'coil':colors[0],
                'bsheet':colors[1],
                'helix':colors[2]
            }
        
        if structure == 'all':
            if std is True:
                coil = Data(df=df, x=df.index, y=df['coil'], color=colors['coil'], label='Coil', fill_between=df['coil_std'],linewidth=linewidth, alpha=alpha, scatter=False)
                bsheet = Data(df=df, x=df.index, y=df['bsheet'], color=colors['bsheet'], label=r'$\beta$-Strand', fill_between=df['bsheet_std'], alpha=alpha,linewidth=linewidth, scatter=False)
                helix = Data(df=df, x=df.index, y=df['helix'], color=colors['helix'], label='Helix', fill_between=df['helix_std'],linewidth=linewidth, alpha=alpha, scatter=False)
            else:
                coil = Data(df=df, x=df.index, y=df['coil'], color=colors['coil'], label='Coil',linewidth=linewidth, alpha=alpha, scatter=False)
                bsheet = Data(df=df, x=df.index, y=df['bsheet'], color=colors['bsheet'], label=r'$\beta$-Strand',linewidth=linewidth, alpha=alpha, scatter=False)
                helix = Data(df=df, x=df.index, y=df['helix'], color=colors['helix'], label='Helix',linewidth=linewidth, alpha=alpha, scatter=False)
        else:
            if structure == 'coil':
                if std is True:
                    coil = Data(df=df, x=df.index, y=df['coil'], color=colors['coil'], label='Coil', fill_between=df['coil_std'],linewidth=linewidth, alpha=alpha, scatter=False)
                    bsheet = None
                    helix = None
                else:
                    coil = Data(df=df, x=df.index, y=df['coil'], color=colors['coil'], linewidth=linewidth, alpha=alpha, scatter=False)
                    bsheet = None
                    helix = None
            if structure == 'bsheet':
                if std is True:
                    coil = None
                    bsheet = Data(df=df, x=df.index, y=df['bsheet'], color=colors['bsheet'], label=r'$\beta$-Strand',  fill_between=df['bsheet_std'],linewidth=linewidth, alpha=alpha, scatter=False)
                    helix = None
                else:
                    coil = None
                    bsheet = Data(df=df, x=df.index, y=df['bsheet'], color=colors['bsheet'], label=r'$\beta$-Strand', linewidth=linewidth, alpha=alpha, scatter=False)
                    helix = None
            if structure == 'helix':
                if std is True:
                    coil = None
                    bsheet = None
                    helix = Data(df=df, x=df.index, y=df['helix'], color=colors['helix'], label='Helix', fill_between=df['helix_std'],linewidth=linewidth, alpha=alpha, scatter=False)
                else:
                    coil = None
                    bsheet = None
                    helix = Data(df=df, x=df.index, y=df['helix'], color=colors['helix'], label='Helix', linewidth=linewidth, alpha=alpha, scatter=False)

        fig = None

        data = [coil, bsheet, helix]



        xticks = ElementParam(xmin=0, xmax=round(df.index[-1]), locs=100, fontsize=16, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=100, locs=20, fontsize=20, minor_locs=5)

        xlabel = ElementParam(label='Time (ns)', fontsize=22, weight='bold')
        ylabel = ElementParam(label='Percentage (%)', fontsize=22, weight='bold')

        if title is not None:
            title = ElementParam(title=title, fontsize=22, weight='bold')
        else:
            title = None

        if gridline == True:
            axes = ElementParam(off=False, semiopen=False, grid=True)
        elif gridline == False:
            axes = ElementParam(off=False, semiopen=True, grid=False)
        else:
            axes = ElementParam(off=False, semiopen=True, grid=False)

        if (legend is not None) and (legend is not False):
            legend = ElementParam(loc='upper right', fontsize=14)

        annotations = annotations

        saveto = output

        return cls('timeseries', data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def dsspPerResidue(cls, dfs, colors, title=None, xlabels=None, output=None):
        data = []
        s = 500
        i = 0
        labels = {
            'alpha':'Helix',
            'beta':r'$\beta$-Strand',
            'coil':'Coil'
        }
        for df in dfs:
            label = labels[df.name]
            d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label=label, s=s, lineconnect=True)
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

        #legend = ElementParam(loc='upper right', markerscale=0.60, ncol=3, individual=False)
        legend = None

        lineconnect_data = []
        linestyles = ['solid', 'dashed', 'dotted']
        i = 0
        lineconnect_colors = ['#000000'] * len(dfs)
        for df in dfs:
            lbl = labels[df.name]
            d = Data(df=df, x=df.index, y=df['mean'], color=lineconnect_colors[i], linestyle=linestyles[i], label=lbl)
            lineconnect_data.append(d)
            i += 1

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto, lineconnect_data)

    @classmethod
    def gyration(cls, df, title=None, colors=None, label=None, output=None, ymin=0, ymax=10, ncol=1, legend=True):
        data = []
        i = 0
        for column in df.columns:
            if column != 'std':
                if isinstance(colors, dict):
                    if 'std' in df.columns:
                        d = Data(df=df, x=df.index, y=df[column], color=colors[column], label='Replicate {}'.format(i+1), fill_between=df['std'])
                    else:
                        d = Data(df=df, x=df.index, y=df[column], color=colors[column], label='Replicate {}'.format(i+1))
                if isinstance(colors, list):
                    if 'std' in df.columns:
                        d = Data(df=df, x=df.index, y=df[column], color=colors[i], label='Replicate {}'.format(i+1))
                    else:
                        d = Data(df=df, x=df.index, y=df[column], color=colors[i], label='Replicate {}'.format(i+1), fill_between=df['std'])
                if isinstance(colors, str):
                    d = Data(df=df, x=df.index, y=df[column], color=colors, label='Replicate {}'.format(i+1), fill_between=df['std'])
                data.append(d)
            i += 1

        fig = None

        xticks = ElementParam(xmin=0, xmax=round(df.index[-1]), locs=100, fontsize=18, minor_locs=20)
        yticks = ElementParam(ymin=ymin, ymax=ymax, locs=5, minor_locs=1, fontsize=18)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label=r'Radius of Gyration ($\AA$)', fontsize=20)

        title = ElementParam(title=title, fontsize=24)

        axes=ElementParam(off=False, semiopen=True)

        if (legend is not False) and (legend is not None):
            legend = ElementParam(loc='upper right', ncol=ncol, fontsize=16)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def rmsdPerPeptide(cls, df, colors=None, title=None, output=None, ymax=1, legend=False, ncol=1):
        # colors = {
        #     'peptide1':'#11174b',
        #     'peptide2':'#3967a4',
        #     'peptide3':'#62bed2',
        #     'peptide4':'#4e3b6c',
        #     'peptide5':'#705c93',
        #     'peptide6':'#be98fc'
        # }
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

        xticks = ElementParam(xmin=0, xmax=round(df.index[-1]), locs=100, fontsize=16, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=ymax, locs=1, fontsize=14, minor_locs=0.1)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label='|$\mu$| (D)', fontsize=20)

        if title is not None:
            title = ElementParam(title=title, fontsize=24)
        else:
            title=None

        axes=ElementParam(off=False, semiopen=True)

        annotations = None
        
        if legend is not False:
            legend = ElementParam(loc='upper right', ncol=ncol)

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def peptideBondDipoleMoments(cls, df, output, fill_between=False, labels=None, title=None, colors=None, ymin=None, ymax=None, legend=False, ncol=2):
        data = []
        i = 0
        _ymax = None
        _ymin = None
        for column in df.columns:
            if 'stdev' in column:
                continue
            if labels is None:
                label = column
            else:
                label = labels[i]
            if fill_between is True:
                if isinstance(colors, dict):
                    d = Data(df=df, x=df.index, y=df[column], color=colors[column], label=label, fill_between=df['stdev'])
                if isinstance(colors, list):
                    d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=label, fill_between=df['stdev'])
            else:
                if isinstance(colors, dict):
                    d = Data(df=df, x=df.index, y=df[column], color=colors[column], label=label)
                if isinstance(colors, list):
                    d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=label)
            data.append(d)
            i += 1
            if _ymax is None:
                _ymax = df[column].max()
            else:
                if df[column].max() > _ymax:
                    _ymax = df[column].max()
            if _ymin is None:
                _ymin = df[column].min()
            else:
                if df[column].min() < _ymin:
                    _ymin = df[column].min()
        if ymax is None:
            ymax = _ymax + 0.5
        if ymin is None:
            ymin = _ymin - 0.5
        fig = ElementParam(width=10, height=6)

        xticks = ElementParam(xmin=0, xmax=round(df.index[-1]), locs=100, fontsize=16, minor_locs=20)
        yticks = ElementParam(ymin=ymin, ymax=ymax, locs=0.5, fontsize=14, minor_locs=0.1)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label=r'|$\mu$| (D)', fontsize=20)

        if title is not None:
            title = ElementParam(title=title, fontsize=22)
        else:
            title=None

        axes=ElementParam(off=False, semiopen=False, grid=True)

        annotations = None
        
        if (legend is not False) and (legend is not None):
            legend = ElementParam(loc='best', ncol=ncol, fontsize=14)

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)

    @classmethod
    def dipoleMomentsOverTime(cls, df, labels=None, colors=None, title=None, output=None, ymin=0, ymax=None, legend=False, ncol=1):
        # colors = {
        #     'peptide1':'#11174b',
        #     'peptide2':'#3967a4',
        #     'peptide3':'#62bed2',
        #     'peptide4':'#4e3b6c',
        #     'peptide5':'#705c93',
        #     'peptide6':'#be98fc'
        # }
        data = []
        i = 0
        _ymax = None
        for column in df.columns:
            if labels is None:
                label = column
            else:
                label = labels[i]
            if isinstance(colors, dict):
                d = Data(df=df, x=df.index, y=df[column], color=colors[column], label=label)
            if isinstance(colors, list):
                d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=label)
            data.append(d)
            i += 1
            if _ymax is None:
                _ymax = df.mean(axis=1).max()
            else:
                if df.mean(axis=1).max() > _ymax:
                    _ymax = df.mean(axis=1).max()
        fig = None

        if ymax is None:
            ymax = round(_ymax) + 10
        if ymax <=1:
            locs=0.2
            minor_locs=0.1
        elif (ymax <= 2) and (ymax>1):
            locs=0.5
            minor_locs=0.1
        elif (ymax > 2):
            locs=1
            minor_locs=0.2
        
        xticks = ElementParam(xmin=0, xmax=round(df.index[-1])+1, locs=100, fontsize=16, minor_locs=20)
        yticks = ElementParam(ymin=ymin, ymax=ymax, locs=locs, fontsize=14, minor_locs=minor_locs)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label='Dipole Moment (D)', fontsize=20)

        if title is not None:
            title = ElementParam(title=title, fontsize=22)
        else:
            title=None

        axes=ElementParam(off=False, semiopen=True)

        annotations = None
        
        if legend is not False:
            legend = ElementParam(loc='best', ncol=ncol, fontsize=14)

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
    def dipoleMomentsPerResidue(cls, dfs, labels=None, lineconnect_dfs=None, lineconnect_labels=None, lineconnect_colors=None, colors=None, title=None, markers=None, xlabels=None, ymin=0, ymax=10, output=None, rotation=None, xlabelfontsize=17, ncol=1):
        data = []
        i = 0
        s = len(dfs)*100
        _ymax = None
        if markers is None:
            markers = ['o'] * len(dfs)
        last_marker = None
        for df in dfs:
            if last_marker is None:
                last_marker = markers[i]
            else:
                if markers[i] != last_marker:
                    s = len(dfs) * 100
            if labels is None:
                if (lineconnect_dfs is None):
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label='Replicate {}'.format(i+1), 
                            s=s, lineconnect=True, marker=markers[i])
                else:
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label='Replicate {}'.format(i+1), 
                            s=s, lineconnect=False, marker=markers[i])
            else:
                if (lineconnect_dfs is None):
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label=labels[i], 
                            s=s, lineconnect=True, marker=markers[i])
                else:
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label=labels[i], 
                            s=s, lineconnect=False, marker=markers[i]) 
            data.append(d)
            i += 1 
            s -= 100
            if _ymax is None:
                _ymax = df['mean'].max()
            else:
                if df['mean'].max() > _ymax:
                    _ymax = df['mean'].max()
        
        lineconnect_data = None
        if lineconnect_dfs is not None:
            lineconnect_data = []
            linestyles = ['solid', 'dashed']
            i = 0
            if lineconnect_colors is None:
                lineconnect_colors = ['#000000'] * len(lineconnect_dfs)
            for df in dfs:
                if lineconnect_labels is None:
                    lbl = '__nolegend__'
                else:
                    lbl = labels[i]
                d = Data(df=df, x=df.index, y=df['mean'], color=lineconnect_colors[i], linestyle=linestyles[i], label=lbl)
                lineconnect_data.append(d)
                i += 1


        fig = ElementParam(width=8, height=6, layout='tight')

        if ymax is None:
            ymax = round(_ymax) + 2

        if xlabels is not None:
            if xlabels[0] != 0:
                xlabels.insert(0,0)

        xticks = ElementParam(locs=1, fontsize=xlabelfontsize, xlabels=xlabels, rotation=rotation)
        yticks = ElementParam(ymin=ymin, ymax=ymax, fontsize=13, locs=5, minor_locs=1)

        xlabel = ElementParam(label='Residue', fontsize=16)
        ylabel = ElementParam(label='|$\mu$| (D)', fontsize=16)

        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=False)

        legend = ElementParam(loc='upper right', markerscale=0.60, ncol=ncol, fontsize=16)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto, lineconnect=lineconnect_data)

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
    def mindistResidueSM(cls, df, colormap='plasma', title=None, xlabels=None, ylabels=None, offset=False, rotation=None, colorbar=True, vmin=0, vmax=1, name=None, annotate=True, output=None):
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
            xticks = ElementParam(locs=1, xlabels=xlabels, offset=0.5, rotation=rotation, fontsize=16)
            yticks = ElementParam(locs=1, ylabels=ylabels, offset=0.5, fontsize=16)
        else:
            xticks = ElementParam(locs=1, xlabels=xlabels)
            yticks = ElementParam(locs=1, ylabels=ylabels)

            

        fig = ElementParam(width=16, height=12)

        title = ElementParam(title=title, fontsize=36)

        legend = ElementParam(colorbar=colorbar)
        
        xlabel = ElementParam()
        ylabel = ElementParam()

        axes = None
        
        if annotate is True:
            anno_data = df.round(2).to_numpy()
            annotations = ElementParam(data=anno_data, fontsize=20)
        else:
            annotations = None

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

    @classmethod
    def sasa(cls, dfs, colors=None, title=None, markers=None, xlabels=None, ymin=0, ymax=300, output=None, rotation=None, xlabelfontsize=17, ncol=1, labels=None, lineconnect_dfs=None, lineconnect_labels=None, lineconnect_colors=None):
        data = []
        i = 0
        s = len(dfs)*100
        _ymax = None
        if markers is None:
            markers = ['o'] * len(dfs)
        last_marker = None
        for df in dfs:
            if last_marker is None:
                last_marker = markers[i]
            else:
                if markers[i] != last_marker:
                    s = len(dfs) * 100
            if labels is None:
                if (lineconnect_dfs is None):
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label='Replicate {}'.format(i+1), 
                            s=s, lineconnect=True, marker=markers[i])
                else:
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label='Replicate {}'.format(i+1), 
                            s=s, lineconnect=False, marker=markers[i])
            else:
                if (lineconnect_dfs is None):
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label=labels[i], 
                            s=s, lineconnect=True, marker=markers[i])
                else:
                    d = Data(df=df, x=df.index, y=df['mean'], yerr=df['stdev'], color=colors[i], label=labels[i], 
                            s=s, lineconnect=False, marker=markers[i]) 
            data.append(d)
            i += 1 
            s -= 100
            if _ymax is None:
                _ymax = df['mean'].max()
            else:
                if df['mean'].max() > _ymax:
                    _ymax = df['mean'].max()
        
        lineconnect_data = None
        if lineconnect_dfs is not None:
            lineconnect_data = []
            linestyles = ['solid', 'dashed']
            i = 0
            if lineconnect_colors is None:
                lineconnect_colors = ['#000000'] * len(lineconnect_dfs)
            for df in lineconnect_dfs:
                if lineconnect_labels is None:
                    lbl = '__nolegend__'
                else:
                    lbl = lineconnect_labels[i]
                d = Data(df=df, x=df.index, y=df['mean'], color=lineconnect_colors[i], linestyle=linestyles[i], label=lbl)
                lineconnect_data.append(d)
                i += 1


        fig = ElementParam(width=8, height=6, layout='tight')

        if ymax is None:
            ymax = round(_ymax) + 2

        if xlabels is not None:
            if xlabels[0] != 0:
                xlabels.insert(0,0)

        xticks = ElementParam(locs=1, fontsize=xlabelfontsize, xlabels=xlabels, rotation=rotation)
        yticks = ElementParam(ymin=ymin, ymax=ymax, fontsize=13, locs=100, minor_locs=20)

        xlabel = ElementParam(label='Residue', fontsize=16)
        ylabel = ElementParam(label='Solvent Accessible Surface Area ($\AA^2$)', fontsize=16)

        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=False)

        legend = ElementParam(loc='best', markerscale=0.60, ncol=ncol, fontsize=16)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto, lineconnect=lineconnect_data)
    
    @classmethod
    def sasaOverTime(cls, df, title=None, colors=None, locs=5, minor_locs=1, label=None, output=None, ymin=0, ymax=10, ncol=1, legend=True):
        data = []
        i = 0
        for column in df.columns:
            if label == 'columns':
                label = column
            else:
                label = 'Replicate {}'.format(i+1)
            if column != 'std':
                if isinstance(colors, dict):
                    d = Data(df=df, x=df.index, y=df[column], color=colors[column], label=label)
                if isinstance(colors, list):
                    d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=label)
                if isinstance(colors, str):
                    d = Data(df=df, x=df.index, y=df[column], color=colors, label=label)
                data.append(d)
            i += 1

        fig = None

        xticks = ElementParam(xmin=0, xmax=1000, locs=100, fontsize=18, minor_locs=20)
        yticks = ElementParam(ymin=ymin, ymax=ymax, locs=locs, minor_locs=minor_locs, fontsize=18)

        xlabel = ElementParam(label='Time (ns)', fontsize=20)
        ylabel = ElementParam(label=r'Solvent Accessible Surface Area ($\AA^2$)', fontsize=20)

        title = ElementParam(title=title, fontsize=22)

        axes=ElementParam(off=False, semiopen=True)

        if (legend is not False) and (legend is not None):
            legend = ElementParam(loc='best', ncol=ncol, fontsize=12)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)
    

    @classmethod
    def heatmap(cls, df, title='', vmin=0, vmax=1, colormap='viridis', xtick_labels=None, ytick_labels=None, legend=True, annotate=False, output=None, major_xticks=None, minor_xticks=None, major_yticks=None, minor_yticks=None, tick_label_fontsize=12, xtick_label_rotation=0, title_fontsize=20, ax_label_fontsize=18, x_label='', y_label='', weight='regular', title_weight='bold', ytick_label_rotation=0, ytick_offset=0, xtick_offset=0, annotation_fontsize=16, hide_yticks=False, hide_xticks=False):
        if colormap is None:
            colormap='viridis'
        data = Data(df=df, vmin=vmin, vmax=vmax, colormap=colormap)
        if xtick_labels is None:
            xtick_labels = df.columns
        xlabs = [0]
        for item in xtick_labels:
            xlabs.append(item)
        xtick_labels = xlabs
        if ytick_labels is None:
            ytick_labels = df.index
            # ylabels = list(reversed(ylabels))
        ylabs = [0]
        for item in ytick_labels:
            ylabs.append(item)
        ytick_labels = ylabs


        xticks = ElementParam(xmin=0, xmax=len(df.columns), locs=major_xticks, fontsize=tick_label_fontsize, minor_locs=minor_xticks, labels=xtick_labels,
                            xtick_label_rotation=xtick_label_rotation, offset=xtick_offset, hide=hide_xticks)

        yticks = ElementParam(ymin=0, ymax=len(df.index), locs=major_yticks, labels=ytick_labels, fontsize=tick_label_fontsize, minor_locs=minor_yticks, ytick_label_rotation=ytick_label_rotation, offset=ytick_offset, hide=hide_yticks)

        fig = None
        title = ElementParam(title=title, fontsize=title_fontsize, weight=title_weight)

        legend = ElementParam(colorbar=legend)
        xlabel = ElementParam(label=x_label, fontsize=ax_label_fontsize, weight=weight)
        ylabel = ElementParam(label=y_label, fontsize=ax_label_fontsize, weight=weight)

        axes = ElementParam(off=False, semiopen=False, grid=False)
        if annotate is True:
            anno_data = df.round(2).to_numpy()
            annotations = ElementParam(data=anno_data, fontsize=annotation_fontsize, atype='heatmap')
        else:
            annotations = None

        saveto = output
        return cls('heatmap',data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)
    
    @classmethod
    def bar(cls, df, x='default', width=1, tick_label=None, color='#9bcccc', alpha=1, title=None, xlabel=None, title_fontsize=20, title_weight='bold', tick_label_fontsize=14, label_fontsize=16, weight='bold', ylabel=None, ymin=0, ymax=None, xmin=None, xmax=None, locs=20, minor_locs=5, xlim=None, capsize=1, capthick=0, elinewidth=1, ecolor=None, orientation='vertical', annotations=True, std=pd.DataFrame(), output=None, bar_labels=None, bar_padding=3, bar_fmt='', bar_fontsize=12):
        if x == 'default':
            x = [i for i in range(len(df.index))]
        elif x == 'index':
            x = df.index
        else:
            if not isinstance(x, Iterable):
                raise ValueError("x must be 'default', 'index', or an iterable")
        if isinstance(df, pd.Series):
            values = np.reshape(df.values, -1)
        elif (isinstance(df, pd.DataFrame)) or (isinstance(df, np.ndarray)):
            if len(df.shape) > 1:
                values = np.reshape(df.values, -1)
            else:
                values = df.values
        if len(values) != len(x):
            raise ValueError('bar coordinates length ({}) does not match the length of values to plot ({})'.format(len(x), len(np.reshape(df.values, -1))))
        error_kw = {
            'capsize':capsize,
            'capthick':capthick,
            'elinewidth':elinewidth,
            'ecolor':ecolor
        }
        if std.empty:
            err = None
        else:
            err = np.reshape(std.values, -1)
        data = [Data(df=df, x=x, values=values, color=color, tick_label=tick_label, width=width, err=err, error_kw=error_kw, orientation=orientation, alpha=alpha)]

        if xmin is None:
            xmin = min(x)
        if xmax is None:
            xmax = max(x)
        if ymin is None:
            ymin = min(values)
        if ymax is None:
            ymin = max(values)

        if orientation == 'vertical':
            _xlocs = None
            _xlocs_minor = None
            _ylocs = locs
            _ylocs_minor = minor_locs
        else:
            _xlocs = locs
            _xlocs_minor = minor_locs
            _ylocs = None
            _ylocs_minor = None
        xticks = ElementParam(labels=tick_label, locs=_xlocs, minor_locs=_xlocs_minor, fontsize=tick_label_fontsize, xmin=xmin,xmax=xmax)
    
        yticks = ElementParam(labels=tick_label, ymin=ymin, ymax=ymax, locs=_ylocs, minor_locs=_ylocs_minor, fontsize=tick_label_fontsize)   

        fig = None

        xlabel = ElementParam(label=xlabel, fontsize=label_fontsize, weight=weight)
        ylabel = ElementParam(label=ylabel, fontsize=label_fontsize, weight=weight)

        if title is not None:
            title = ElementParam(title=title, fontsize=title_fontsize, weight=title_weight)
        else:
            title = ElementParam(title='', fontsize=2, weight=title_weight)

        legend = None

        axes = ElementParam(off=False, semiopen=True)
        if annotations:
            annotations = [ElementParam(atype='autolabel', labels=bar_labels, fmt=bar_fmt, padding=bar_padding, fontsize=bar_fontsize)]
        else:
            annotations = None

        saveto = output

        return cls('bar', data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)
   
   
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
        if colors is None:
            colors = ['#000080', '#a00000']
        data = []
        i = 0
        for column in df.columns:
            if 'std' in column:
                continue
            try:
                d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=labels[i], fill_between='{}_std'.format(column))
            except:
                d = Data(df=df, x=df.index, y=df[column], color=colors[i], label=labels[i])
            i += 1
            data.append(d)
        # coil = Data(df=df, x=df.index, y=df['coil'], color=colors['coil'], label='Coil')
        # bsheet = Data(df=df, x=df.index, y=df['bsheet'], color=colors['bsheet'], label=r'$\beta$-Strand')
        # helix = Data(df=df, x=df.index, y=df['helix'], color=colors['helix'], label='Helix')
        
        fig = None

        xticks = ElementParam(xmin=0, xmax=df.index[-1], locs=100, fontsize=16, minor_locs=20)
        yticks = ElementParam(ymin=0, ymax=100, locs=20, fontsize=20, minor_locs=5)

        xlabel = ElementParam(label='Time (ns)', fontsize=22)
        ylabel = ElementParam(label='Percentage (%)', fontsize=22)

        title = ElementParam(title=title, fontsize=24)

        axes = ElementParam(off=False, semiopen=False, grid=True)

        legend = ElementParam(loc='upper right', fontsize=14)

        annotations = None

        saveto = output

        return cls(data, fig, xticks, yticks, xlabel, ylabel, title, axes, legend, annotations, saveto)


    
    


