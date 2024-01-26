import os
from pymd.plot import PlotData
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import sys
import socket
from collections.abc import Iterable

if os.name == 'posix':
    matplotlib.use('Agg')

class Plotter():

    def __init__(self, nrows=1, ncols=1, sharex=False, sharey=False, w=8, h=6, font='default'):
        if not font == 'default':
            plt.rcParams['font.family'] = font
        self.fig, self.axes = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey)
        self.fig.set_size_inches(w,h)
        self.pdata = PlotData(*[None]*14)
        self.nrows = nrows
        self.ncols = ncols 
        self.sharex=sharex
        self.sharey=sharey
        self.w = w
        self.h = h


    def graph(self, pdata, ax, container=None, **kwargs):
        '''
        Controls the specs of the graph. 
        Arguments:
        - pdata (pymd.plot.PlotData): plot data object
        - ax (matplotlib.axes.Axes): current axis
        - container (matplotlib.axes)
        - **kwargs: mostly just for msc things. such as passing img from pcolor to draw a colorbar. 
        '''
        # labels and axes
        ax.set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
        ax.set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

        if (hasattr(pdata.xticks, 'locs')) and (pdata.xticks.locs is not None):
            ax.xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
        if (hasattr(pdata.xticks, 'labels')) and (pdata.xticks.labels is not None):
            ax.set_xticklabels(pdata.xticks.labels, rotation=pdata.xticks.xtick_label_rotation)
            if hasattr(pdata.xticks, 'offset'):
                if pdata.xticks.offset is not None:
                    dx = pdata.xticks.offset
                    dy = 0
                    offset = transforms.ScaledTranslation(dx, dy, self.fig.dpi_scale_trans)
                    for label in ax.xaxis.get_majorticklabels():
                        label.set_transform(label.get_transform() + offset)
        if (hasattr(pdata.yticks, 'locs')) and (pdata.yticks.locs is not None):
            ax.yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
        if (hasattr(pdata.yticks, 'labels')) and (pdata.yticks.labels is not None):
            ax.set_yticklabels(pdata.yticks.labels, rotation=pdata.yticks.ytick_label_rotation)
            if hasattr(pdata.yticks, 'offset'):
                if pdata.yticks.offset is not None:
                    dx = 0
                    dy = pdata.yticks.offset
                    offset = transforms.ScaledTranslation(dx, dy, self.fig.dpi_scale_trans)
                    for label in ax.yaxis.get_majorticklabels():
                        label.set_transform(label.get_transform() + offset)
        ax.tick_params(axis='both', which='major', labelsize=pdata.xticks.fontsize)
        if (hasattr(pdata.yticks, 'hide')):
            if pdata.yticks.hide:
                ax.set_yticks([])
        if hasattr(pdata.xticks, 'hide'):
            if pdata.xticks.hide:
                ax.set_xticks([])


        if (hasattr(pdata.xticks, 'minor_locs')) and (pdata.xticks.minor_locs is not None):
            ax.xaxis.set_minor_locator(MultipleLocator(pdata.xticks.minor_locs))
        if (hasattr(pdata.yticks, 'minor_locs')) and (pdata.yticks.minor_locs is not None):
            ax.yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))

        # show/hide axes
        if pdata.axes.semiopen is True:
            ax.spines['right'].set_visible(False) # hide right axis
            ax.spines['top'].set_visible(False) # hide top axis
            ax.spines['right'].set_linewidth(20)
            ax.spines['top'].set_linewidth(20)

        # grid
        if (hasattr(pdata.axes, 'grid')):
            if pdata.axes.grid is True:
                ax.grid(b=True, which='major', axis='both', c='black', alpha=0.2)

        # title
        if pdata.title is not None:
            ax.set_title(pdata.title.title, fontsize=pdata.title.fontsize, weight=pdata.title.weight)

        # legend
        if (pdata.legend is not None) and (pdata.legend is not False):
            if pdata.plot_type == 'heatmap':
                if pdata.legend.colorbar:
                    if 'img' in kwargs:
                        cbar = self.fig.colorbar(kwargs['img'])
            else:
                if hasattr(pdata.legend, 'ncol') and hasattr(pdata.legend, 'fontsize'):
                    ax.legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize, ncol=pdata.legend.ncol)
                elif hasattr(pdata.legend, 'fontsize'):
                    ax.legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize)
                elif hasattr(pdata.legend, 'ncol'):
                    ax.legend(loc=pdata.legend.loc, ncol=pdata.legend.ncol)
                else:
                    ax.legend(loc=pdata.legend.loc)

        ax.set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize, weight=pdata.xlabel.weight)
        ax.set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize, weight=pdata.xlabel.weight)


        # annotations
        if pdata.annotations is not None:
            if isinstance(pdata.annotations, Iterable):
                for annotation in pdata.annotations:
                    if annotation.atype == 'plot':
                        ax.plot(annotation.x, annotation.y, linestyle=annotation.linestyle, color=annotation.color, linewidth=annotation.linewidth)
                    if annotation.atype == 'autolabel':
                        ax.bar_label(container, labels=annotation.labels, fmt=annotation.fmt, padding=annotation.padding)
            else:
                annotation = pdata.annotations
                if annotation.atype == 'annotate':
                    plt.annotate(s=annotation.text, xy=annotation.xy, fontsize=annotation.fontsize)
                elif annotation.atype == 'heatmap':
                    for i in range(0, len(pdata.data.df.index)):
                        for j in range(0, len(pdata.data.df.columns)):
                            anno = pdata.annotations.data[i,j]
                            ax.text(j+0.5, i+0.5, anno, ha="center", va="center", color='w', fontsize=pdata.annotations.fontsize)
        else:
            pass
        return ax
    
    def timeseries(self, df, out=None, panel=False, ax=None, show=False, titles=[], suptitle=None, **kwargs):
        if panel is True:
            if isinstance(df, pd.DataFrame):
                for (i, col) in enumerate(df.columns):
                    if titles != []:
                        kwargs['title'] = titles[i]
                    data = pd.DataFrame()
                    data[col] = df[col]
                    self.axes.flat[i] = self.timeseries(data, panel=False, ax=self.axes.flat[i], **kwargs)
            elif isinstance(df, (list, tuple, np.ndarray)):
                for (i, data) in enumerate(df):
                    if titles != []:
                        kwargs['title'] = titles[i]
                    self.axes.flat[i] = self.timeseries(data, panel=False, ax=self.axes.flat[i], **kwargs)
            self._fix_labels()
            if suptitle is not None:
                plt.suptitle(suptitle, weight='bold', fontsize=18)
            plt.tight_layout()
            if out is not None:
                if out.endswith('png'):
                    plt.savefig(out, dpi=300)
                if out.endswith('svg'):
                    plt.savefig(out, dpi=300)
                print('Plotted {}'.format(out))
            if show:
                plt.show()
            return self.axes
        else:
            if ax is None:
                ax = self.axes
            if not isinstance(df, (list, tuple, set)):
                dfs = [df]
            else:
                dfs = df
            for df in dfs:
                pdata = PlotData.timeseries(df, **kwargs)
                self.pdata = pdata
                for d in pdata.data:
                    ax.plot(d.x, d.y, color=d.color, label=d.label, linestyle=d.linestyle, linewidth=d.linewidth, alpha=d.alpha)
                    if hasattr(d, 'fill_between'):
                        ax.fill_between(d.x, d.y-d.fill_between, d.y+d.fill_between, alpha=0.2, color=d.color)
                    if d.scatter:
                        ax.scatter(d.x, d.y, color=d.color, marker=d.marker, s=d.s)
                    ax = self.graph(pdata, ax)
        if panel is False:
            plt.tight_layout()
            if out is not None:
                if out.endswith('png'):
                    plt.savefig(out, dpi=300)
                if out.endswith('svg'):
                    plt.savefig(out, dpi=300)
                print('Plotted {}'.format(out))
                if show is True:
                    plt.show()
            elif show is False:
                pass
            elif show is True:
                plt.show()
            else:
                plt.show()
                plt.close()
        return ax

    def heatmap(self, df, panel=False, ax=None, show=False, **kwargs):
        if panel is True:
            for (i, col) in enumerate(df.columns):
                data = pd.DataFrame()
                data[col] = df[col]
                self.axes.flat[i] = self.heatmap(data, panel=False, ax=self.axes.flat[i], **kwargs)
            self._fix_labels()
            return self.axes
        else:
            if ax is None:
                ax = self.axes
            if not isinstance(df, (list, tuple, set)):
                dfs = [df]
            else:
                dfs = df
            for df in dfs:
                pdata = PlotData.heatmap(df, **kwargs)
                self.pdata = pdata
                img = ax.pcolor(pdata.data.df, vmin=pdata.data.vmin, vmax=pdata.data.vmax, cmap=pdata.data.colormap)
                ax = self.graph(pdata, ax, img=img)
        if panel is False:
            plt.tight_layout()
            if pdata.saveto is not None:
                if pdata.saveto.endswith('png'):
                    plt.savefig(pdata.saveto, dpi=300)
                if pdata.saveto.endswith('svg'):
                    plt.savefig(pdata.saveto, dpi=300)
                print('Plotted {}'.format(pdata.saveto))
                if show is True:
                    plt.show()
            elif show is False:
                pass
            elif show is True:
                plt.show()
            else:
                plt.show()
                plt.close()
        return ax

    def histogram(self, df, out=None, panel=False, show=True, ax=None, suptitle=None, **kwargs):
        if panel is True:
            for (i, col) in enumerate(df.columns):
                data = pd.DataFrame()
                data[col] = df[col]
                self.axes.flat[i] = self.histogram(data, panel=False, ax=self.axes.flat[i], **kwargs)
            self._fix_labels()
            if suptitle is not None:
                plt.suptitle(suptitle, weight='bold', fontsize=18)
            plt.tight_layout()
            if out is not None:
                if out.endswith('png'):
                    plt.savefig(out, dpi=300)
                if out.endswith('svg'):
                    plt.savefig(out, dpi=300)
                print('Plotted {}'.format(out))
            if show:
                plt.show()
            return self.axes
        else:
            if ax is None:
                ax = self.axes
            pdata = PlotData.histogram(df, **kwargs)
            self.pdata = pdata
            for d in pdata.data:
                ax.hist(d.x, bins=d.bins, density=d.density, alpha=d.alpha, color=d.color, label=d.label)
                ax = self.graph(pdata, ax)
            if suptitle is not None:
                plt.suptitle(suptitle, weight='bold', fontsize=18)
            plt.tight_layout()
            if pdata.saveto is not None:
                if pdata.saveto.endswith('png'):
                    plt.savefig(pdata.saveto, dpi=300)
                if pdata.saveto.endswith('svg'):
                    plt.savefig(pdata.saveto, dpi=300)
                print('Plotted {}'.format(pdata.saveto))
                if show is True:
                    plt.show()
            elif show is False:
                pass
            elif show is True:
                plt.show()
            else:
                plt.show()
                plt.close()
        return ax

    def bar(self, df, panel=False, ax=None, show=False, **kwargs):
        if panel is True:
            for (i, col) in enumerate(df.columns):
                data = pd.DataFrame()
                data[col] = df[col]
                self.axes.flat[i] = self.bar(data, panel=False, ax=self.axes.flat[i], **kwargs)
            self._fix_labels()
            return self.axes
        else:
            if ax is None:
                ax = self.axes
            if not isinstance(df, (list, tuple, set)):
                dfs = [df]
            else:
                dfs = df
            for df in dfs:
                pdata = PlotData.bar(df, **kwargs)
                self.pdata = pdata
                for d in pdata.data:
                    if d.orientation == 'vertical':
                        rects = ax.bar(d.x, d.values, color=d.color, tick_label=d.tick_label, width=d.width,alpha=d.alpha, yerr=d.err, error_kw=d.error_kw)
                    elif d.orientation == 'horizontal':
                        rects = ax.barh(d.x, d.values, color=d.color, tick_label=d.tick_label, width=d.width,alpha=d.alpha, xerr=d.err, error_kw=d.error_kw)
                    ax = self.graph(pdata, ax, container=rects)
        if panel is False:
            plt.tight_layout()
            if pdata.saveto is not None:
                if pdata.saveto.endswith('png'):
                    plt.savefig(pdata.saveto, dpi=300)
                if pdata.saveto.endswith('svg'):
                    plt.savefig(pdata.saveto, dpi=300)
                print('Plotted {}'.format(pdata.saveto))
                if show is True:
                    plt.show()
            elif show is False:
                pass
            elif show is True:
                plt.show()
            else:
                plt.show()
                plt.close()
        return ax
    def panel(self, _type, df, ax, **kwargs):
        pass

    def _fix_labels(self):
        row = 0
        for x in range(0, len(self.axes.flat)):
            if x % self.ncols == 0:
                row += 1
            else:
                if self.sharey:
                    self.axes.flat[x].set_ylabel(None)
            if (row != self.nrows) and (self.sharex):
                self.axes.flat[x].set_xlabel(None)
                # pdata.xlabel.label = None
        plt.tight_layout()

    def show(self):
        if (os.name == 'posix'):
            raise OSError('A Linux operating system has been detected! Plot cannot be shown.')
        plt.show()
        plt.close()
        # reset self
        self.fig, self.axes = plt.subplots(self.nrows, self.ncols, sharex=self.sharex, sharey=self.sharey)
        self.fig.set_size_inches(self.w,self.h)
        self.pdata = PlotData(*[None]*14)
    
    def save(self, output=None, dpi=300):
        if output is None:
            if self.pdata.saveto is not None:
                output = self.pdata.saveto
            else:
                raise ValueError('No output file specified')
        if output.endswith('png'):
            plt.savefig(output, dpi=dpi)
        if output.endswith('svg'):
            plt.savefig(output)
        print('Plotted {}'.format(output))
        plt.close()
        # reset self
        self.fig, self.axes = plt.subplots(self.nrows, self.ncols, sharex=self.sharex, sharey=self.sharey)
        self.fig.set_size_inches(self.w,self.h)
        self.pdata = PlotData(*[None]*14)

