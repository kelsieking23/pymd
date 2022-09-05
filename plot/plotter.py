import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import sys
import socket

if os.name == 'posix':
    matplotlib.use('Agg')

# from pymd.mdanalysis.postprocess import PostProcess
# from pymd.plot.data import PlotData

# sys.path.append('/work/cascades/kelsieking23/iapp_analysis/scripts/python')
# from mdanalysis.postprocess import PostProcess
# from plot.data import PlotData


class Plotter:
    def __init__(self, system=None):
        self.system = system


    def timeseries(self, pdata, fig=None, ax=None, show=True):
        # init plot
        if (show is True) and (os.name == 'posix'):
            print('A linux operating system is detected, plot output will not be shown')
            show = False
        panel = False
        if fig is None:
            fig, ax = plt.subplots()

            if pdata.fig is None:
                fig.set_size_inches(8,6)
        else:
            if pdata.fig is not None:
                fig.set_size_inches(pdata.fig.width, pdata.fig.height)
            panel = True
        
        p = 0
        for d in pdata.data:
            if d is None:
                continue
            if hasattr(d, 'linestyle'):
                linestyle=d.linestyle
            else:
                linestyle='solid'
            if hasattr(d, 'linewidth'):
                ax.plot(d.x, d.y, color=d.color, label=d.label, linestyle=linestyle, linewidth=d.linewidth)
            else:
                ax.plot(d.x, d.y, color=d.color, label=d.label, linestyle=linestyle)
            if hasattr(d, 'fill_between'):
                ax.fill_between(d.x, d.y-d.fill_between, d.y+d.fill_between, alpha=0.2, color=d.color)
            p += 1

        # labels and axes
        ax.set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
        ax.set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

        if (hasattr(pdata.xticks, 'locs')) and (pdata.xticks.locs is not None):
            ax.xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
        # plt.xticks(fontsize=pdata.xticks.fontsize)
        if (hasattr(pdata.yticks, 'locs')) and (pdata.yticks.locs is not None):
            ax.yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
        # plt.yticks(fontsize=pdata.yticks.fontsize)
        ax.tick_params(axis='both', which='major', labelsize=pdata.xticks.fontsize)

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
            ax.set_title(pdata.title.title, fontsize=pdata.title.fontsize)

        # legend
        if (pdata.legend is not None) and (pdata.legend is not False):
            if hasattr(pdata.legend, 'ncol') and hasattr(pdata.legend, 'fontsize'):
                ax.legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize, ncol=pdata.legend.ncol)
            elif hasattr(pdata.legend, 'fontsize'):
                ax.legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize)
            elif hasattr(pdata.legend, 'ncol'):
                ax.legend(loc=pdata.legend.loc, ncol=pdata.legend.ncol)
            else:
                ax.legend(loc=pdata.legend.loc)

        ax.set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
        ax.set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

        # annotations
        if pdata.annotations is not None:
            for annotation in pdata.annotations:
                if annotation.atype == 'plot':
                    ax.plot(annotation.x, annotation.y, linestyle=annotation.linestyle, color=annotation.color, linewidth=annotation.linewidth)
                if annotation.atype == 'annotate':
                    plt.annotate(s=annotation.text, xy=annotation.xy, fontsize=annotation.fontsize)

        # plt.show()
        # saving
        
        if panel is False:
            plt.tight_layout()
            if pdata.saveto is not None:
                if pdata.saveto.endswith('png'):
                    plt.savefig(pdata.saveto, dpi=300)
                if pdata.saveto.endswith('svg'):
                    plt.savefig(pdata.saveto, dpi=300)
                print('Plotted {}'.format(pdata.saveto))
                plt.close()
            elif show is False:
                pass
            elif show is True:
                plt.show()
            else:
                plt.show()
                plt.close()
        
        # plt.show()
        return ax
    
    @staticmethod
    def minmax(xlim, ylim, _min_x, _max_x, _min_y, _max_y):
        # print(xlim, ylim)
        # print(_min_x, _max_x, _min_y, _max_y)

        if _min_x is None:
            _min_x = xlim[0]
        elif _min_x > xlim[0]:
            _min_x = xlim[0]
        else:
            pass

        if _max_x is None:
            _max_x = xlim[1]
        elif _max_x < xlim[1]:
            _max_x = xlim[1]
        else:
            pass

        if _min_y is None:
            _min_y = ylim[0]
        elif _min_y > ylim[0]:
            _min_y = ylim[0]
        else:
            pass

        if _max_y is None:
            _max_y = ylim[1]
        elif _max_y < ylim[1]:
            _max_y = ylim[1]
        else:
            pass
        # print(_min_x, _max_x, _min_y, _max_y)
        # print('**************')
        return _min_x, _max_x, _min_y, _max_y

    def timeseriesPanel(self, data, ncols, nrows, title=None, axes=[], sharex=True, sharey=True, output=None, legend=True, legend_data=None, show=True):
        fig, ax = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey)
        # fig = plt.figure() 
        normal_w = 8
        normal_h = 6
        fig_h = (normal_h * nrows) 
        fig_w = (normal_w * ncols) + 2
        fig.set_size_inches(fig_w, fig_h)
        # gs = GridSpec(ncols=ncols,nrows=nrows)
        _axes = []
        ax_index = 0
        _min_x = None
        _max_x = None
        _min_y = None
        _max_y = None
        _xlabel = None
        _ylabel = None
        rows = []
        i = 0
        row_counter = -1
        if (nrows > 1) and (ncols > 1):
            for x in range(0, nrows):
                row_counter += 1
                rows.append([])
                for y in range(0, ncols):
                    try:
                        pdata = data[i]
                        if sharex:
                            _xlabel = pdata.xlabel.label
                            pdata.xlabel.label = None
                        if sharey:
                            _ylabel = pdata.ylabel.label
                            pdata.ylabel.label = None
                    except:
                        break
                    if axes == []:
                        ax[x,y] = self.timeseries(pdata, fig=fig, ax=ax[x,y])
                        _min_x, _max_x, _min_y, _max_y = self.minmax(ax[x,y].get_xlim(), ax[x,y].get_ylim(), _min_x, _max_x, _min_y, _max_y)
                        _axes.append(ax[x,y])
                        rows[row_counter].append(ax[x,y])
                        i += 1
                    else:
                        ax[x,y] = axes[ax_index]
                        _axes.append(axes[ax_index])
                        i += 1
                        ax_index += 1
        else:
            if nrows > 1:
                indeces = nrows
            else:
                indeces = ncols
            for x in range(0, indeces):
                pdata = data[i]
                if sharex:
                    _xlabel = pdata.xlabel.label
                    pdata.xlabel.label = None
                if sharey:
                    _ylabel = pdata.ylabel.label
                    pdata.ylabel.label = None
                if axes == []:
                    ax[x] = self.timeseries(pdata, fig=fig, ax=ax[x], show=False)
                    _min_x, _max_x, _min_y, _max_y = self.minmax(ax[x].get_xlim(), ax[x].get_ylim(), _min_x, _max_x, _min_y, _max_y)
                    _axes.append(ax[x])
                    i += 1
                else:
                    ax[x] = axes[ax_index]
                    _axes.append(axes[ax_index])
                    i += 1
                    ax_index += 1

        import matplotlib.lines as mlines
        if (legend_data is not None):
            if (pdata.legend_data is not None) and (legend_data is None):
                legend_data = pdata.legend_data
            if legend is True:
                handles = []
                labels = []
                for label, color in legend_data:
                    if label == 'bstrand':
                        label = r'$\beta$-Strand'
                    handles.append(mlines.Line2D([], [], color=color, label=label))
                    labels.append(label)
                plt.figlegend(handles, labels, loc='lower center', ncol=4, fontsize=18)
            # plt.legend(handles=handles, bbox_to_anchor=(-0.81, -0.15, 1.6, .102), loc='lower left',
            # ncol=3, mode="expand", borderaxespad=0.)

        # fix axes?
        for ax_ in _axes:
            ax_.set_xlim(_min_x, _max_x)
            ax_.set_ylim(_min_y, _max_y)
        
        if sharex:
            for _ax in rows[-1]:
                _ax.set_xlabel(_xlabel, fontsize=pdata.xlabel.fontsize)
        if sharey:
            for row in rows:
                print(_ylabel)
                row[0].set_ylabel(_ylabel, fontsize=pdata.ylabel.fontsize)
        # if (len(data) % 2) != 0:
        #     ax[h-1,w-1].remove()
        # layout
        if title is not None:
            fig.suptitle(title, fontsize=20)
        plt.tight_layout()
        # fig.subplots_adjust(top=0.85)
        #saving
        if output is not None:
            if output.endswith('png'):
                plt.savefig(output, dpi=300)
            if output.endswith('svg'):
                plt.savefig(output)
            print('Plotted {}'.format(output))
            if show is True:
                plt.show()
        else:
            plt.show()
    
    def markerPlot(self, pdata):

        # init plot
        fig, ax = plt.subplots()
        fig.set_size_inches(8,6)
        
        z = 1
        for d in pdata.data:
            errorbars = ax.errorbar(x=d.x, y=d.y, xerr=None, yerr=d.yerr, marker='o',capsize=12, elinewidth=1, capthick=1.5, ms=7, fillstyle='full', fmt='none', color=d.color, label=None, zorder=z)
            z += 1
            if (hasattr(d, 'lineconnect')):
                if d.lineconnect is True:
                    lineconnect = ax.plot(d.x, d.y, color=d.color, label='_nolegend_', zorder=z)
                    z += 1
            points = ax.scatter(d.x, d.y, s=d.s, color=d.color, label=d.label, zorder=z)
            z += 1
        if hasattr(pdata, 'lineconnect'):
            if d.lineconnect is False:
                z = 0
                for d in pdata.lineconnect:
                    ax.plot(d.x, d.y, color=d.color, label=d.label, linestyle=d.linestyle, zorder=z)

        # labels and axes
        if (hasattr(pdata.xticks, 'xmin')) and (hasattr(pdata.xticks, 'xmax')):
            ax.set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
        if (hasattr(pdata.yticks, 'ymin')) and (hasattr(pdata.yticks, 'ymax')):
            ax.set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

        if hasattr(pdata.xticks, 'locs'):
            ax.xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
        if hasattr(pdata.yticks, 'locs'):
            ax.yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
        
        if hasattr(pdata.xticks, 'minor_locs'):
            ax.xaxis.set_minor_locator(MultipleLocator(pdata.xticks.minor_locs))
        if hasattr(pdata.yticks, 'minor_locs'):
            ax.yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))

        if hasattr(pdata.xticks, 'xlabels'):
            ax.set_xticklabels(pdata.xticks.xlabels)
            if hasattr(pdata.xticks, 'rotation'):
                plt.setp(ax.xaxis.get_majorticklabels(), rotation=pdata.xticks.rotation)
            if hasattr(pdata.xticks, 'pad'):
                ax.xaxis.labelpad=pdata.xticks.pad
            if hasattr(pdata.xticks, 'offset'):
                if pdata.xticks.offset is not None:
                    dx = pdata.xticks.offset
                    dy = 0
                    offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                    for label in ax.xaxis.get_majorticklabels():
                        label.set_transform(label.get_transform() + offset)
        if hasattr(pdata.yticks, 'ylabels'):
            ax.set_yticklabels(pdata.yticks.ylabels)

        plt.xticks(fontsize=pdata.xticks.fontsize)
        plt.yticks(fontsize=pdata.yticks.fontsize)

        # axis paramters
        if pdata.axes.semiopen is True:
            ax.spines['right'].set_visible(False) # hide right axis
            ax.spines['top'].set_visible(False) # hide top axis
            ax.spines['right'].set_linewidth(10)
            ax.spines['top'].set_linewidth(10)

        # titles
        if hasattr(pdata.title, 'title'):
            ax.set_title(pdata.title.title, fontsize=pdata.title.fontsize)
        ax.set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
        ax.set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

        #legend
        if pdata.legend is not None:
            if (hasattr(pdata.legend, 'individual')):
                if pdata.legend.individual is True:
                    ax.legend(loc=pdata.legend.loc, markerscale=pdata.legend.markerscale, fontsize=8, ncol=pdata.legend.ncol)
                else:
                    pass
            else:
                ax.legend(loc=pdata.legend.loc, markerscale=pdata.legend.markerscale, fontsize=10, ncol=pdata.legend.ncol)

        if hasattr(pdata.fig, 'layout'):
            if pdata.fig.layout == 'tight':
                plt.tight_layout()
        # saving
        if pdata.saveto.endswith('png'):
            plt.savefig(pdata.saveto, dpi=300)
        else:
            plt.savefig(pdata.saveto)
        plt.close()
        print('Plotted {}'.format(pdata.saveto))
        return pdata.saveto

    def markerPlotPanel(self, data, output, title, legend_data=None):
        # init plot
        h = int(int(len(data)) / 2)
        if len(data) % 2 != 0:
            h = h + 1
        w = 2
        fig, ax = plt.subplots(h, w)
        normal_w = 8
        normal_h = 6
        fig_h = (normal_h * h) 
        fig_w = (normal_w * w) + 2
        fig.set_size_inches(fig_w, fig_h)
        plotnum = 0
        i = 0
        if len(data) > 2:
            for x in range(0, h):
                for y in range(0, w):
                    try:
                        pdata = data[plotnum]
                    except:
                        continue
                    z = 0
                    for d in pdata.data:
                        print(d.x, '\n', d.y, '\n', d.yerr)
                        errorbars = ax[x,y].errorbar(x=d.x, y=d.y, xerr=None, yerr=d.yerr, marker='o',capsize=12, elinewidth=1, capthick=1.5, ms=7, fillstyle='full', fmt='none', color=d.color, label=None, zorder=z)
                        z += 1
                        if (hasattr(d, 'lineconnect')):
                            print('plotting...')
                            if d.lineconnect is True:
                                print('its true')
                                lineconnect = ax[x,y].plot(d.x, d.y, color=d.color, label='_nolegend_', zorder=z)
                                z += 1
                        points = ax[x,y].scatter(d.x, d.y, s=d.s, color=d.color, label=d.label, zorder=z)
                        z += 1
                    if hasattr(pdata, 'lineconnect'):
                        if d.lineconnect is False:
                            z = 0
                            for d in pdata.lineconnect:
                                ax[x,y].plot(d.x, d.y, color=d.color, label=d.label, linestyle=d.linestyle, zorder=z)

                    # labels and axes
                    if (hasattr(pdata.xticks, 'xmin')) and (hasattr(pdata.xticks, 'xmax')):
                        ax[x,y].set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
                    if (hasattr(pdata.yticks, 'ymin')) and (hasattr(pdata.yticks, 'ymax')):
                        ax[x,y].set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

                    if hasattr(pdata.xticks, 'locs'):
                        ax[x,y].xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
                    if hasattr(pdata.yticks, 'locs'):
                        ax[x,y].yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
                    
                    if hasattr(pdata.xticks, 'minor_locs'):
                        ax[x,y].xaxis.set_minor_locator(MultipleLocator(pdata.xticks.minor_locs))
                    if hasattr(pdata.yticks, 'minor_locs'):
                        ax[x,y].yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))

                    if hasattr(pdata.xticks, 'xlabels'):
                        ax[x,y].set_xticklabels(pdata.xticks.xlabels)
                    if hasattr(pdata.yticks, 'ylabels'):
                        ax[x,y].set_yticklabels(pdata.yticks.ylabels)

                    ax[x,y].tick_params(axis='x', which='major', labelsize=pdata.xticks.fontsize)
                    ax[x,y].tick_params(axis='y', which='major', labelsize=pdata.yticks.fontsize)

                    # axis paramters
                    if pdata.axes.semiopen is True:
                        ax[x,y].spines['right'].set_visible(False) # hide right axis
                        ax[x,y].spines['top'].set_visible(False) # hide top axis
                        ax[x,y].spines['right'].set_linewidth(10)
                        ax[x,y].spines['top'].set_linewidth(10)

                    # titles
                    if hasattr(pdata.title, 'title'):
                        ax[x,y].set_title(pdata.title.title, fontsize=pdata.title.fontsize)
                    ax[x,y].set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
                    ax[x,y].set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

                    # legend
                    if pdata.legend is not None:
                        if (hasattr(pdata.legend, 'individual')):
                            if pdata.legend.individual is True:
                                ax[x,y].legend(loc=pdata.legend.loc, markerscale=pdata.legend.markerscale, fontsize=14)
                            else:
                                pass
                        else:
                            ax[x,y].legend(loc=pdata.legend.loc, markerscale=pdata.legend.markerscale, fontsize=14)


                    plotnum += 1
        else:
            for x in range(0, w):
                try:
                    pdata = data[plotnum]
                except:
                    continue
                z = 0
                for d in pdata.data:
                    print(d.x, '\n', d.y, '\n', d.yerr)
                    errorbars = ax[x].errorbar(x=d.x, y=d.y, xerr=None, yerr=d.yerr, marker='o',capsize=12, elinewidth=1, capthick=1.5, ms=7, fillstyle='full', fmt='none', color=d.color, label=None, zorder=z)
                    z += 1
                    if (hasattr(d, 'lineconnect')):
                        print('plotting...')
                        if d.lineconnect is True:
                            print('its true')
                            lineconnect = ax[x].plot(d.x, d.y, color=d.color, label='_nolegend_', zorder=z)
                            z += 1
                    points = ax[x].scatter(d.x, d.y, s=d.s, color=d.color, label=d.label, zorder=z)
                    z += 1
                if hasattr(pdata, 'lineconnect'):
                    if d.lineconnect is False:
                        z = 0
                        for d in pdata.lineconnect:
                            ax[x].plot(d.x, d.y, color=d.color, label=d.label, linestyle=d.linestyle, zorder=z)

                # labels and axes
                if (hasattr(pdata.xticks, 'xmin')) and (hasattr(pdata.xticks, 'xmax')):
                    ax[x].set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
                if (hasattr(pdata.yticks, 'ymin')) and (hasattr(pdata.yticks, 'ymax')):
                    ax[x].set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

                if hasattr(pdata.xticks, 'locs'):
                    ax[x].xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
                if hasattr(pdata.yticks, 'locs'):
                    ax[x].yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
                
                if hasattr(pdata.xticks, 'minor_locs'):
                    ax[x].xaxis.set_minor_locator(MultipleLocator(pdata.xticks.minor_locs))
                if hasattr(pdata.yticks, 'minor_locs'):
                    ax[x].yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))

                if hasattr(pdata.xticks, 'xlabels'):
                    ax[x].set_xticklabels(pdata.xticks.xlabels)
                if hasattr(pdata.yticks, 'ylabels'):
                    ax[x].set_yticklabels(pdata.yticks.ylabels)

                ax[x].tick_params(axis='x', which='major', labelsize=pdata.xticks.fontsize)
                ax[x].tick_params(axis='y', which='major', labelsize=pdata.yticks.fontsize)

                # axis paramters
                if pdata.axes.semiopen is True:
                    ax[x].spines['right'].set_visible(False) # hide right axis
                    ax[x].spines['top'].set_visible(False) # hide top axis
                    ax[x].spines['right'].set_linewidth(10)
                    ax[x].spines['top'].set_linewidth(10)

                # titles
                if hasattr(pdata.title, 'title'):
                    ax[x].set_title(pdata.title.title, fontsize=pdata.title.fontsize)
                ax[x].set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
                ax[x].set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

                # legend
                if pdata.legend is not None:
                    if (hasattr(pdata.legend, 'individual')):
                        if pdata.legend.individual is True:
                            ax[x].legend(loc=pdata.legend.loc, markerscale=pdata.legend.markerscale, fontsize=14)
                        else:
                            pass
                    else:
                        ax[x].legend(loc=pdata.legend.loc, markerscale=pdata.legend.markerscale, fontsize=14)


                plotnum += 1
        import matplotlib.lines as mlines
        if legend_data is not None:
            handles = []
            labels = []
            for label, color, size in legend_data:
                handles.append(mlines.Line2D([], [], color=color, label=label, linestyle=None, marker='o', markersize=size))
                labels.append(label)
            plt.figlegend(handles, labels, loc='lower center', ncol=4, fontsize=18)

        # layout
        fig.suptitle(title, fontsize=28)
        plt.tight_layout(pad=5, h_pad=2, w_pad=2)
        if '\n' in title:
            plt.subplots_adjust(top=0.85)
        # saving
        if output.endswith('png'):
            plt.savefig(output, dpi=300)
        else:
            plt.savefig(output)
        plt.close()
        print('Plotted {}'.format(output))

    def heatmap(self, pdata):
        fig, ax = plt.subplots()
        if pdata.fig is None:
            fig.set_size_inches(16,12)
        else:
            fig.set_size_inches(pdata.fig.width, pdata.fig.height)
        img = ax.pcolor(pdata.data.df, vmin=pdata.data.vmin, vmax=pdata.data.vmax, cmap=pdata.data.colormap)
        # img = ax.imshow(pdata.data.df, vmin=pdata.data.vmin, vmax=pdata.data.vmax, cmap=pdata.data.colormap)
        if hasattr(pdata.legend, 'colorbar'):
            if pdata.legend.colorbar is True:
                cbar = fig.colorbar(img)
        if hasattr(pdata.xticks, 'locs'):
            ax.xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
        if hasattr(pdata.yticks, 'locs'):
            ax.yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
        if hasattr(pdata.xticks, 'xlabels'):
            ax.xaxis.set_ticklabels(pdata.xticks.xlabels)
            if hasattr(pdata.xticks, 'rotation'):
                plt.setp(ax.xaxis.get_majorticklabels(), rotation=pdata.xticks.rotation)
            if hasattr(pdata.xticks, 'offset'):
                if pdata.xticks.offset is not None:
                    dx = pdata.xticks.offset
                    dy = 0
                    offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                    for label in ax.xaxis.get_majorticklabels():
                        label.set_transform(label.get_transform() + offset)
            if hasattr(pdata.xticks, 'fontsize'):
                ax.tick_params(axis='x', which='major', labelsize=pdata.xticks.fontsize)
        if (hasattr(pdata.yticks, 'ylabels')) and (pdata.yticks.ylabels is not None):
            ax.yaxis.set_ticklabels(pdata.yticks.ylabels)
            if hasattr(pdata.yticks, 'offset'):
                if pdata.yticks.offset is not None:
                    dx = 0
                    dy = pdata.yticks.offset
                    offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                    for label in ax.yaxis.get_majorticklabels():
                        label.set_transform(label.get_transform() + offset)
            if hasattr(pdata.yticks, 'fontsize'):
                ax.tick_params(axis='y', which='major', labelsize=pdata.yticks.fontsize)

        # annotations
        if pdata.annotations is not None:
            for i in range(0, len(pdata.data.df.index)):
                for j in range(0, len(pdata.data.df.columns)):
                    anno = pdata.annotations.data[i,j]
                    ax.text(j+0.5, i+0.5, anno, ha="center", va="center", color='w', fontsize=pdata.annotations.fontsize)

        #titling
        if hasattr(pdata.title, 'title'):
            ax.set_title(pdata.title.title, fontsize=pdata.title.fontsize)
        plt.tick_params(
            axis='both',
            which='both',
            bottom=False,
            left=False
        )

        # if hasattr(pdata.xlabel, 'label'):
        #     ax.set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
        # if hasattr(pdata.ylabel, 'label'):
        #     ax.set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)
        ax.set_aspect('equal')
        if pdata.saveto is not None:
            if pdata.saveto.endswith('png'):
                plt.savefig(pdata.saveto, dpi=300)
            if pdata.saveto.endswith('svg'):
                plt.savefig(pdata.saveto, dpi=300)
            plt.close()
            print('Plotted {}'.format(pdata.saveto))
        else:
            # plt.show()
            pass
        return fig, ax
        # return pdata.saveto
        # if os.name == 'nt':
        #     plt.show()

    def heatmapPanel(self, data, output, title, legend_data=None):
        h = int(int(len(data)) / 2)
        w = 2
        fig, ax = plt.subplots(h, w)
        normal_h = 8
        normal_w = 6
        fig_h = (normal_h * h) 
        fig_w = (normal_w * w) + 6
        fig.set_size_inches(fig_w, fig_h)
        plotnum = 0
        if h > 1:
            for x in range(0, h):
                print(0,h)
                for y in range(0, w):
                    pdata = data[plotnum]
                    img = ax[x,y].pcolor(pdata.data.df, vmin=pdata.data.vmin, vmax=pdata.data.vmax, cmap=pdata.data.colormap)
                    
                    if hasattr(pdata.xticks, 'locs'):
                        ax[x,y].xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
                    if hasattr(pdata.yticks, 'locs'):
                        ax[x,y].yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
                    if hasattr(pdata.xticks, 'xlabels'):
                        ax[x,y].xaxis.set_ticklabels(pdata.xticks.xlabels)
                        if hasattr(pdata.xticks, 'rotation'):
                            plt.setp(ax[x,y].xaxis.get_majorticklabels(), rotation=pdata.xticks.rotation)
                        if hasattr(pdata.xticks, 'offset'):
                            if pdata.xticks.offset is not None:
                                dx = pdata.xticks.offset
                                dy = 0
                                offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                                for label in ax[x,y].xaxis.get_majorticklabels():
                                    label.set_transform(label.get_transform() + offset)
                        if hasattr(pdata.xticks, 'fontsize'):
                            ax[x,y].tick_params(axis='x', which='major', labelsize=pdata.xticks.fontsize)
                    if (hasattr(pdata.yticks, 'ylabels')) and (pdata.yticks.ylabels is not None):
                        ax[x,y].yaxis.set_ticklabels(pdata.yticks.ylabels)
                        if hasattr(pdata.yticks, 'offset'):
                            if pdata.yticks.offset is not None:
                                dx = 0
                                dy = pdata.yticks.offset
                                offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                                for label in ax[x,y].yaxis.get_majorticklabels():
                                    label.set_transform(label.get_transform() + offset)
                        if hasattr(pdata.yticks, 'fontsize'):
                            ax[x,y].tick_params(axis='y', which='major', labelsize=pdata.yticks.fontsize)

                    # annotations
                    if pdata.annotations is not None:
                        print(pdata.data.df)
                        for i in range(0, len(pdata.data.df.index)):
                            for j in range(0, len(pdata.data.df.columns)):
                                anno = pdata.annotations.data[i,j]
                                ax[x,y].text(j+0.5, i+0.5, anno, ha="center", va="center", color='w', fontsize=pdata.annotations.fontsize)

                    # titling
                    if hasattr(pdata.title, 'title'):
                        ax[x,y].set_title(pdata.title.title, fontsize=pdata.title.fontsize)

                    if hasattr(pdata.xlabel, 'label'):
                        ax[x,y].set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
                    if hasattr(pdata.ylabel, 'label'):
                        ax[x,y].set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

                    plotnum += 1
        else:
            for x in range(0, w):
                pdata = data[plotnum]
                img = ax[x].pcolor(pdata.data.df, vmin=pdata.data.vmin, vmax=pdata.data.vmax, cmap=pdata.data.colormap)
                
                if hasattr(pdata.xticks, 'locs'):
                    ax[x].xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
                if hasattr(pdata.yticks, 'locs'):
                    ax[x].yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
                if hasattr(pdata.xticks, 'xlabels'):
                    ax[x].xaxis.set_ticklabels(pdata.xticks.xlabels)
                    if hasattr(pdata.xticks, 'rotation'):
                        plt.setp(ax[x].xaxis.get_majorticklabels(), rotation=pdata.xticks.rotation)
                    if hasattr(pdata.xticks, 'offset'):
                        if pdata.xticks.offset is not None:
                            dx = pdata.xticks.offset
                            dy = 0
                            offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                            for label in ax[x].xaxis.get_majorticklabels():
                                label.set_transform(label.get_transform() + offset)
                    if hasattr(pdata.xticks, 'fontsize'):
                        ax[x].tick_params(axis='x', which='major', labelsize=pdata.xticks.fontsize)
                if (hasattr(pdata.yticks, 'ylabels')) and (pdata.yticks.ylabels is not None):
                    ax[x].yaxis.set_ticklabels(pdata.yticks.ylabels)
                    if hasattr(pdata.yticks, 'offset'):
                        if pdata.yticks.offset is not None:
                            dx = 0
                            dy = pdata.yticks.offset
                            offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                            for label in ax[x].yaxis.get_majorticklabels():
                                label.set_transform(label.get_transform() + offset)
                    if hasattr(pdata.yticks, 'fontsize'):
                        ax[x].tick_params(axis='y', which='major', labelsize=pdata.yticks.fontsize)

                # annotations
                if pdata.annotations is not None:
                    print(pdata.data.df)
                    for i in range(0, len(pdata.data.df.index)):
                        for j in range(0, len(pdata.data.df.columns)):
                            anno = pdata.annotations.data[i,j]
                            ax[x].text(j+0.5, i+0.5, anno, ha="center", va="center", color='w', fontsize=pdata.annotations.fontsize)

                # titling
                if hasattr(pdata.title, 'title'):
                    ax[x].set_title(pdata.title.title, fontsize=pdata.title.fontsize)

                if hasattr(pdata.xlabel, 'label'):
                    ax[x].set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
                if hasattr(pdata.ylabel, 'label'):
                    ax[x].set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

                plotnum += 1
        # layout
        fig.suptitle(title, fontsize=24)
        plt.tight_layout(pad=2, h_pad=2, w_pad=2)

        if hasattr(pdata.legend, 'colorbar'):
            print('colorbar!!!!!')
            if pdata.legend.colorbar is True:
                cbar = fig.colorbar(img, ax=ax, orientation='horizontal')
                cbar.ax.tick_params(labelsize=10)
        #saving
        if output.endswith('png'):
            plt.savefig(output, dpi=300)
        if output.endswith('svg'):
            plt.savefig(output)
        print('Plotted {}'.format(output))

    def histogram(self, pdata, show=True):
        fig, ax = plt.subplots()

        if pdata.fig is None:
            fig.set_size_inches(8,6)
        else:
            fig.set_size_inches(pdata.fig.width, pdata.fig.height)

        for d in pdata.data:
            if (hasattr(d, 'alpha')) and (hasattr(d, 'color')):
                ax.hist(d.x, bins=d.bins, density=d.density, alpha=d.alpha, color=d.color, label=d.label)
            elif hasattr(d, 'alpha'):
                ax.hist(d.x, bins=d.bins, density=d.density, alpha=d.alpha,label=d.label)
            elif hasattr(d, 'color'):
                ax.hist(d.x, bins=d.bins, density=d.density, color=d.color, label=d.label)
            else:
                ax.hist(d.x, bins=d.bins, density=d.density, label=d.label)
            

        # labels and axes
        ax.set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
        ax.set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

        if hasattr(pdata.xticks, 'locs'):
            ax.xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
        # plt.xticks(fontsize=pdata.xticks.fontsize)
        if hasattr(pdata.yticks, 'locs'):
            ax.yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
        # plt.yticks(fontsize=pdata.yticks.fontsize)
        ax.tick_params(axis='both', which='major', labelsize=pdata.xticks.fontsize)

        if hasattr(pdata.xticks, 'minor_locs'):
            ax.xaxis.set_minor_locator(MultipleLocator(pdata.xticks.minor_locs))
        if hasattr(pdata.yticks, 'minor_locs'):
            ax.yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))


        if pdata.axes.semiopen is True:
            ax.spines['right'].set_visible(False) # hide right axis
            ax.spines['top'].set_visible(False) # hide top axis
            ax.spines['right'].set_linewidth(10)
            ax.spines['top'].set_linewidth(10)

        if (hasattr(pdata.axes, 'grid')):
            if pdata.axes.grid is True:
                ax.grid(b=True, which='major', axis='both', c='black', alpha=0.2)

        ax.set_title(pdata.title.title, fontsize=pdata.title.fontsize)

        ax.set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
        ax.set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

        if (pdata.legend is not None) and (pdata.legend is not False):
            if hasattr(pdata.legend, 'ncol') and hasattr(pdata.legend, 'fontsize'):
                ax.legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize, ncol=pdata.legend.ncol)
            elif hasattr(pdata.legend, 'fontsize'):
                ax.legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize)
            elif hasattr(pdata.legend, 'ncol'):
                ax.legend(loc=pdata.legend.loc, ncol=pdata.legend.ncol)
            else:
                ax.legend(loc=pdata.legend.loc)

        # annotations
        if pdata.annotations is not None:
            for annotation in pdata.annotations:
                if annotation.atype == 'plot':
                    ax.plot(annotation.x, annotation.y, linestyle=annotation.linestyle, color=annotation.color, linewidth=annotation.linewidth)
                if annotation.atype == 'annotate':
                    plt.annotate(s=annotation.text, xy=annotation.xy, fontsize=annotation.fontsize)

        if pdata.saveto is not None:
            print(pdata.saveto)
            if pdata.saveto.endswith('png'):
                plt.savefig(pdata.saveto, dpi=300)
            if pdata.saveto.endswith('svg'):
                plt.savefig(pdata.saveto, dpi=300)
            print('Plotted {}'.format(pdata.saveto))
        else:
            print('noshow')
            print(pdata.saveto)
            plt.show()
        if show:
            plt.show()
        plt.close()

    def histogramPanel(self, data, output, title):
        if (len(data) % 2) == 0:
            h = int(int(len(data)) / 2)
        else:
            h = int(int(len(data)+1) / 2)
        w = 2
        fig, ax = plt.subplots(h, w)
        normal_w = 8
        normal_h = 6
        fig_h = (normal_h * h) 
        fig_w = (normal_w * w) + 2
        fig.set_size_inches(fig_w, fig_h)
        i = 0
        if h > 1:
            for x in range(0, h):
                for y in range(0, w):
                    try:
                        pdata = data[i]
                    except:
                        break
                    for d in pdata.data:
                        if (hasattr(d, 'alpha')) and (hasattr(d, 'color')):
                            ax[x,y].hist(d.x, bins=d.bins, density=d.density, alpha=d.alpha, color=d.color, label=d.label)
                        elif hasattr(d, 'alpha'):
                            ax[x,y].hist(d.x, bins=d.bins, density=d.density, alpha=d.alpha, label=d.label)
                        elif hasattr(d, 'color'):
                            ax[x,y].hist(d.x, bins=d.bins, density=d.density, color=d.color, label=d.label)
                        else:
                            ax[x,y].hist(d.x, bins=d.bins, density=d.density, label=d.label)
                      

                    # labels and axes
                    ax[x,y].set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
                    ax[x,y].set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

                    if hasattr(pdata.xticks, 'locs'):
                        ax[x,y].xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
                    # plt.xticks(fontsize=pdata.xticks.fontsize)
                    if hasattr(pdata.yticks, 'locs'):
                        ax[x,y].yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
                    # plt.yticks(fontsize=pdata.yticks.fontsize)
                    ax[x,y].tick_params(axis='both', which='major', labelsize=pdata.xticks.fontsize)

                    if hasattr(pdata.xticks, 'minor_locs'):
                        ax[x,y].xaxis.set_minor_locator(MultipleLocator(pdata.xticks.minor_locs))
                    if hasattr(pdata.yticks, 'minor_locs'):
                        ax[x,y].yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))


                    if pdata.axes.semiopen is True:
                        ax[x,y].spines['right'].set_visible(False) # hide right axis
                        ax[x,y].spines['top'].set_visible(False) # hide top axis
                        ax[x,y].spines['right'].set_linewidth(10)
                        ax[x,y].spines['top'].set_linewidth(10)

                    if (hasattr(pdata.axes, 'grid')):
                        if pdata.axes.grid is True:
                            ax[x,y].grid(b=True, which='major', axis='both', c='black', alpha=0.2)

                    ax[x,y].set_title(pdata.title.title, fontsize=pdata.title.fontsize)

                    ax[x,y].set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
                    ax[x,y].set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

                    if (pdata.legend is not None) and (pdata.legend is not False):
                        if hasattr(pdata.legend, 'ncol') and hasattr(pdata.legend, 'fontsize'):
                            ax[x,y].legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize, ncol=pdata.legend.ncol)
                        elif hasattr(pdata.legend, 'fontsize'):
                            ax[x,y].legend(loc=pdata.legend.loc, fontsize=pdata.legend.fontsize)
                        elif hasattr(pdata.legend, 'ncol'):
                            ax[x,y].legend(loc=pdata.legend.loc, ncol=pdata.legend.ncol)
                        else:
                           ax[x,y].legend(loc=pdata.legend.loc)
                    # annotations
                    if pdata.annotations is not None:
                        for annotation in pdata.annotations:
                            if annotation.atype == 'plot':
                                ax[x,y].plot(annotation.x, annotation.y, linestyle=annotation.linestyle, color=annotation.color, linewidth=annotation.linewidth)
                            if annotation.atype == 'annotate':
                                plt.annotate(s=annotation.text, xy=annotation.xy, fontsize=annotation.fontsize)
                    i += 1
        else:
            for x in range(0, w):
                pdata = data[i]

                for d in pdata.data:
                    if (hasattr(d, 'alpha')) and (hasattr(d, 'color')):
                        ax[x].hist(d.x, bins=d.bins, density=d.density, alpha=d.alpha, color=d.color)
                    elif hasattr(d, 'alpha'):
                        ax[x].hist(d.x, bins=d.bins, density=d.density, alpha=d.alpha)
                    elif hasattr(d, 'color'):
                        ax[x].hist(d.x, bins=d.bins, density=d.density, color=d.color)
                    else:
                        ax[x].hist(d.x, bins=d.bins, density=d.density)
                  
                # labels and axes
                ax[x].set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
                ax[x].set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)

                if hasattr(pdata.xticks, 'locs'):
                    ax[x].xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
                # plt.xticks(fontsize=pdata.xticks.fontsize)
                if hasattr(pdata.yticks, 'locs'):
                    ax[x].yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
                # plt.yticks(fontsize=pdata.yticks.fontsize)
                ax[x].tick_params(axis='both', which='major', labelsize=pdata.xticks.fontsize)

                if hasattr(pdata.xticks, 'minor_locs'):
                    ax[x].xaxis.set_minor_locator(MultipleLocator(pdata.xticks.minor_locs))
                if hasattr(pdata.yticks, 'minor_locs'):
                    ax[x].yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))


                if pdata.axes.semiopen is True:
                    ax[x].spines['right'].set_visible(False) # hide right axis
                    ax[x].spines['top'].set_visible(False) # hide top axis
                    ax[x].spines['right'].set_linewidth(10)
                    ax[x].spines['top'].set_linewidth(10)

                ax[x].set_title(pdata.title.title, fontsize=pdata.title.fontsize)


                ax[x].set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize)
                ax[x].set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize)

                # annotations
                if pdata.annotations is not None:
                    for annotation in pdata.annotations:
                        if annotation.atype == 'plot':
                            ax[x].plot(annotation.x, annotation.y, linestyle=annotation.linestyle, color=annotation.color, linewidth=annotation.linewidth)
                        if annotation.atype == 'annotate':
                            plt.annotate(s=annotation.text, xy=annotation.xy, fontsize=annotation.fontsize)
                i += 1

        # remove last empty axis if any
        if (len(data) % 2) != 0:
            ax[h-1,w-1].remove()
        # layout
        fig.suptitle(title, fontsize=20)
        plt.tight_layout(pad=7, h_pad=2, w_pad=2)
        fig.subplots_adjust(top=0.9)
        #saving
        if output is not None:
            if output.endswith('png'):
                plt.savefig(output, dpi=300)
            if output.endswith('svg'):
                plt.savefig(output)
            print('Plotted {}'.format(output))   
        else:
            plt.show()
        plt.close()    

    def barPlot(self, pdata):

        # initialize fig, ax
        fig, ax = plt.subplots()
        if pdata.fig is None:
            fig.set_size_inches(8,6)
        else:
            fig.set_size_inches(pdata.fig.width, pdata.fig.height)
        
        # plot data
        rects = []
        if not isinstance(pdata.data, list):
            if pdata.data.color is None:
                if 'stdev' in pdata.data.df.columns:
                    rect = ax.bar(pdata.data.xpos, pdata.data.df['mean'], yerr=pdata.data.df['stdev'], width=pdata.data.width, label=pdata.data.label, ecolor='#5c5c5c', capsize=pdata.data.capsize, color='#a0b7d2')
                else:
                    rect = ax.bar(pdata.data.xpos, pdata.data.df['mean'], width=pdata.data.width, label=pdata.data.label, ecolor='#5c5c5c', capsize=pdata.data.capsize, color='#a0b7d2')
            else:
                if 'stdev' in pdata.data.df.columns:
                    rect = ax.bar(pdata.data.xpos, pdata.data.df['mean'], yerr=pdata.data.df['stdev'], width=pdata.data.width, label=pdata.data.label, ecolor='#5c5c5c', capsize=pdata.data.capsize, color=pdata.data.color)
                else:
                    rect = ax.bar(pdata.data.xpos, pdata.data.df['mean'], width=pdata.data.width, label=pdata.data.label, ecolor='#5c5c5c', color=pdata.data.color)
                # rect = ax.bar(pdata.data.xpos, pdata.data.df['mean'], width=pdata.data.width, label=pdata.data.label, ecolor='#5c5c5c', capsize=pdata.data.capsize, color=pdata.data.color)
            rects.append(rect)
        else:
            for data in pdata.data:
                print(data)
                if data.color is None:
                    rect = ax.bar(data.xpos, data.df['mean'], yerr=data.df['stdev'], width=data.width, label=data.label, ecolor='#5c5c5c', capsize=data.capsize, color='#a0b7d2')
                else:
                    rect = ax.bar(data.xpos, data.df['mean'], yerr=data.df['stdev'], width=data.width, label=data.label, ecolor='#5c5c5c', capsize=data.capsize, color=data.color)
                rects.append(rect)
        ax.axhline(color='black', linewidth=0.8)
        # annotations
        if pdata.annotations is not None:
            for annotation in pdata.annotations:
                if annotation.atype == 'autolabel':
                    for rect in rects:
                        caplines = rect.errorbar.lines[1][1]._y
                        self.autolabel(rect, ax, caplines=caplines, fontsize=annotation.fontsize)
        

        # labels and ticks
        if (hasattr(pdata.xticks, 'xmin')) and (hasattr(pdata.xticks, 'xmax')): 
            ax.set_xlim(pdata.xticks.xmin, pdata.xticks.xmax)
        if (hasattr(pdata.yticks, 'ymin')) and (hasattr(pdata.yticks, 'ymax')):
            ax.set_ylim(pdata.yticks.ymin, pdata.yticks.ymax)


        # ax.set_xticks(MultipleLocator(pdata.xticks.locs))
        # ax.set_xticklabels(pdata.xticks.labels)
        ax.xaxis.set_major_locator(MultipleLocator(pdata.xticks.locs))
        plt.xticks(fontsize=pdata.xticks.fontsize)
        if hasattr(pdata.yticks, 'locs'):
            ax.yaxis.set_major_locator(MultipleLocator(pdata.yticks.locs))
        if hasattr(pdata.yticks, 'minor_locs'):
            ax.yaxis.set_minor_locator(MultipleLocator(pdata.yticks.minor_locs))
        plt.yticks(fontsize=pdata.yticks.fontsize)

        if pdata.xlabel.label is not None:
            ax.set_xlabel(pdata.xlabel.label, fontsize=pdata.xlabel.fontsize, weight='bold')
        if pdata.ylabel.label is not None:
            ax.set_ylabel(pdata.ylabel.label, fontsize=pdata.ylabel.fontsize, weight='bold')

        # axes
        if pdata.axes.semiopen is True:
            ax.spines['right'].set_visible(False) # hide right axis
            ax.spines['top'].set_visible(False) # hide top axis
            ax.spines['right'].set_linewidth(10)
            ax.spines['top'].set_linewidth(10)

        # title
        if pdata.title.title is not None:
            ax.set_title(pdata.title.title, fontsize=pdata.title.fontsize)

        # legend
        if pdata.legend is not None:
            ax.legend(loc=pdata.legend.loc)

        # saving
        if pdata.saveto is not None:
            if pdata.saveto.endswith('png'):
                plt.savefig(pdata.saveto, dpi=600)
            if pdata.saveto.endswith('svg'):
                plt.savefig(pdata.saveto, dpi=300)
            plt.close()
            print('Plotted {}'.format(pdata.saveto))
        else:
            if os.name == 'nt':
                plt.show()
        return ax

    def barPlotPanel(self, data, title, output=None):
        h = int(int(len(data)) / 2)
        w = 2
        fig, ax = plt.subplots(h, w)
        normal_w = 8
        normal_h = 6
        fig_h = (normal_h * h) 
        fig_w = (normal_w * w) + 2
        fig.set_size_inches(fig_w, fig_h)
        i = 0
        if h > 1:
            for x in range(0, h):
                for y in range(0, w):
                    pdata = data[i]
                    ax[x,y] = self.barPlot(pdata)
                    i += 1
        else:
            for x in range(0, w):
                pdata = data[i]
                ax[x] = self.barPlot(pdata)
                i += 1

        fig.suptitle(title, fontsize=30)
        plt.tight_layout(pad=7, h_pad=2, w_pad=2)
        fig.subplots_adjust(top=0.85)
        #saving
        if output is not None:
            if output.endswith('png'):
                plt.savefig(output, dpi=300)
            if output.endswith('svg'):
                plt.savefig(output)
            print('Plotted {}'.format(output))
        else:
            if os.name == 'nt':
                plt.show()
            
    def rmsf(self, marker=True, lineplot=False):
        aggregate = []
        for rep in range(1, self.system.reps+1):
            files = []
            root = self.system.directory[rep]['root']
            if self.system.peptides is not None:
                for pep in range(1, self.system.peptides+1):
                    filename = os.path.join(root, 'rmsf_{}.xvg'.format(str(pep)))
                    files.append(filename)
            else:
                filename = os.path.join(root, 'rmsf.xvg')
                files.append(filename)
            color_cycler = ['#11174b', '#2d5e9e', '#62bed2']
            fig, ax = plt.subplots()
            i = 0
            dfs = []
            ms = 13
            for filename in files:
                # get data
                df = pd.read_csv(filename, delim_whitespace=True, skiprows=17, header=None)
                df.columns=['residue', 'rmsf']
                # fix for glycines
                resnums = []
                k = 0
                for item in df['residue']:
                    resnum = int(item)
                    if len(resnums) == 0:
                        resnums.append(resnum)
                    else:
                        if resnum != (resnums[-1] + 1):
                            temp = pd.DataFrame({'residue':(resnum-1), 'rmsf':float(0)}, index=[k])
                            df = pd.concat([df.iloc[0:k], temp, df.iloc[k:]]).reset_index(drop=True)
                        resnums.append(resnum)
                    k += 1
                k = 0
                for item in df['residue']:
                    if item == 19:
                        df.drop(k, inplace=True)
                    if item == 30:
                        df.drop(k, inplace=True)
                    k += 1
                df.reset_index(inplace=True)

                # plot 
                dfs.append(df)
                aggregate.append(df)
                if marker == True:
                    label = 'Peptide {}'.format(str(i + 1))
                    lines = {'linestyle': 'None'}
                    plt.rc('lines', **lines)
                    ax.errorbar(df['residue'], df['rmsf'], xerr=None, yerr=None, marker='o', color=color_cycler[i], label=label, ms=ms)
                else:
                    ax.plot(df['residue'], df['rmsf'], color=color_cycler[i])
                i += 1
                ms -= 3
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.set_ylim(None,1)
            ax.set_xlabel('Residue')
            ax.set_ylabel('RMSF (nm)')
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            saveto = os.path.join(self.system.directory[rep]['root'], 'rmsf.png')
            plt.savefig(saveto, dpi=300)
            plt.close()
            
            # plot by peptide if specified
            if self.system.peptides is not None:
                stdevs = []
                for i in range(len(df['rmsf'])):
                    values = []
                    k = 0
                    for df in dfs:
                        values.append(df['rmsf'][i])
                    stdevs.append(np.std(values))
                concat = pd.concat(dfs)
                by_row_index = concat.groupby(concat.index)
                mean = by_row_index.mean()
                fig, ax = plt.subplots()
                if marker == True:
                    lines = {'linestyle': 'None'}
                    plt.rc('lines', **lines)
                    ax.errorbar(mean['residue'], mean['rmsf'], xerr=None, yerr=stdevs, marker='o',capsize=5, elinewidth=0.5, capthick=0.5, color=color_cycler[1], ms=7, fillstyle='full')
                else:
                    ax.plot(df['residue'], df['rmsf'], color=color_cycler[i])
                ax.set_ylim(None,1)
                ax.set_xlabel('Residue')
                ax.set_ylabel('RMSF (nm)')
                ax.yaxis.set_minor_locator(MultipleLocator(0.1))
                ax.xaxis.set_major_locator(MultipleLocator(1))
                saveto = os.path.join(self.system.directory[rep]['root'], 'rmsf_average.png')
                plt.savefig(saveto, dpi=300)
                plt.close()
        
        stdevs = []
        for i in range(len(df['rmsf'])):
            values = []
            k = 0
            for df in aggregate:
                values.append(df['rmsf'][i])
            stdevs.append(np.std(values))
        concat = pd.concat(aggregate)
        concat = pd.concat(dfs)
        by_row_index = concat.groupby(concat.index)
        mean = by_row_index.mean()
        fig, ax = plt.subplots()
        if marker == True:
            lines = {'linestyle': 'None'}
            plt.rc('lines', **lines)
            ax.errorbar(mean['residue'], mean['rmsf'], xerr=None, yerr=stdevs, marker='o',capsize=5, capthick=0.5, elinewidth=0.5, color=color_cycler[1], ms=7, fillstyle='full')
        else:
            ax.plot(df['residue'], df['rmsf'], color=color_cycler[i])
        ax.set_ylim(None,1)
        ax.set_xlabel('Residue')
        ax.set_ylabel('RMSF (nm)')
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.xaxis.set_major_locator(MultipleLocator(1))
        saveto = os.path.join(self.system.root, 'rmsf_average.png')
        plt.savefig(saveto, dpi=300)
        plt.close()

    def mindist(self, matrix, dtype, rep=None, title=None, average_title=None, output=None, testing=False, outfiletype='svg'):
        fig, ax = plt.subplots()
        fig.set_size_inches(8,6)
        img = ax.pcolor(matrix, vmin=0, vmax=20, cmap='plasma')
        cbar = fig.colorbar(img)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

        # titling
        if rep is not None:
            if (isinstance(rep, int)):
                ax.set_title('Replicate {}'.format(rep), fontsize=24)
            else:
                if title is not None:
                    ax.set_title(title, fontsize=24)
                    print(title)

        if output is None:
            if isinstance(rep, int):
                saveto = os.path.join(self.system.directory['images']['mindist'], 'mindist_{}_{}_vmax_20.{}'.format(dtype, rep, outfiletype))
            else:
                saveto = os.path.join(self.system.directory['images']['mindist'], 'mindist_{}_average_vmax_20.{}'.format(dtype, outfiletype))
        plt.savefig(saveto, dpi=300, format=outfiletype)
        plt.close()
        print('Plotted {}'.format(saveto))
        return saveto
    
    def dsspPerResidue(self, df, plot, xpm, rep=None, peptides=None, average=False):
        color_cycler = ['#11174b', '#4479bf', '#62bed2']
        fig, ax = plt.subplots()
        fig.set_size_inches(8,6)
        if plot == 'line':
            i = 2
            for column in df.columns:
                ax.plot(df[column], color=color_cycler[i], label=column)
                i -= 1
            # locs = [float(i) for i in range(1, len(df.index)+1)]
            # plt.xticks(locs)
            plt.xticks(fontsize=8)
            plt.xlim(1, len(df.index))
            plt.xlabel('Residue')
            plt.ylabel('Secondary Sructure Probability')
            # plt.legend()
            # plt.show()
        if plot == 'bar':
            print('i am terrified')
            for index in df.index:
                alpha = df.loc[index, 'alpha']
                beta = df.loc[index, 'beta']
                coil = df.loc[index, 'coil']
                tup = [('alpha', alpha), ('beta', beta), ('coil', coil)]
                lst = len(tup)  
                for i in range(0, lst):  
                    
                    for j in range(0, lst-i-1):  
                        if (tup[j][1] > tup[j + 1][1]):  
                            temp = tup[j]  
                            tup[j]= tup[j + 1]  
                            tup[j + 1]= temp  
                tup.reverse()
                if index == 1:
                    if tup[0][0] == 'alpha':
                        top_color = color_cycler[1]
                    if tup[0][0] == 'beta':
                        top_color = color_cycler[0]
                    if tup[0][0] == 'coil':
                        top_color = color_cycler[2]
                    if tup[1][0] == 'alpha':
                        mid_color = color_cycler[1]
                    if tup[1][0] == 'beta':
                        mid_color = color_cycler[0]
                    if tup[1][0] == 'coil':
                        mid_color = color_cycler[2]
                    if tup[2][0] == 'alpha':
                        low_color = color_cycler[1]
                    if tup[2][0] == 'beta':
                        low_color = color_cycler[0]
                    if tup[2][0] == 'coil':
                        low_color = color_cycler[2]
                    ax.bar(x=index, height=tup[0][1], label=tup[0][0], color=top_color)
                    ax.bar(x=index, height=tup[1][1], label=tup[1][0], color=mid_color)
                    ax.bar(x=index, height=tup[2][1], label=tup[2][0], color=low_color)
                else:
                    if tup[0][0] == 'alpha':
                        top_color = color_cycler[1]
                    if tup[0][0] == 'beta':
                        top_color = color_cycler[0]
                    if tup[0][0] == 'coil':
                        top_color = color_cycler[2]
                    if tup[1][0] == 'alpha':
                        mid_color = color_cycler[1]
                    if tup[1][0] == 'beta':
                        mid_color = color_cycler[0]
                    if tup[1][0] == 'coil':
                        mid_color = color_cycler[2]
                    if tup[2][0] == 'alpha':
                        low_color = color_cycler[1]
                    if tup[2][0] == 'beta':
                        low_color = color_cycler[0]
                    if tup[2][0] == 'coil':
                        low_color = color_cycler[2]
                    ax.bar(x=index, height=tup[0][1], color=top_color)
                    ax.bar(x=index, height=tup[1][1], color=mid_color)
                    ax.bar(x=index, height=tup[2][1], color=low_color)

            locs = [float(i) for i in range(1, len(df.index)+1)]
            plt.xticks(locs)
            plt.xticks(fontsize=8)
            plt.xlabel('Residue')
            plt.ylabel('Secondary Sructure Probability')
        # if rep is None:
        #     xpm_base = xpm[:-3] + 'png'
        #     saveto = xpm_base
        if average is False:
            xpm_base = xpm.split('.')[0]
            saveto = os.path.join(self.system.directory[rep]['dssp']['root'], '{}.png'.format(xpm_base))
        else:
            xpm_base = xpm.split('.')[0] + '_average.png'
            saveto = os.path.join(self.system.root, xpm_base)
        plt.savefig(saveto, dpi=300)
        plt.close()
        print('Plotted {}'.format(saveto))
        # plt.legend()
        # plt.show()
    
    def dsspOverTime(self, df, xvg, rep=None, output=None, legend=True, title=None, average_title=None):
        '''
        Plot DSSP over time for given data.
        Arguments:
        ** df (DataFrame): dataframe containing secondary structure information (output of PostProcess.dsspOverTime())
        ** xvg (str): .xvg filename
        ** rep (int, optional): replicate number
        ** output (str, optional): output file path
        ** legend (bool, default True): if true will show legend on upper right corner
        '''
        colors = {
            'coil_percent':'#11174b',
            'bsheet_percent':'#62bed2',
            'helix_percent':'#3967a4'
        }

        # init plot
        fig, ax = plt.subplots()
        fig.set_size_inches(8,6)

        ax.plot(df.index, df['coil_percent'], color=colors['coil_percent'], label='Coil')
        ax.plot(df.index, df['bsheet_percent'], color=colors['bsheet_percent'], label=r'$\beta$-Strand')
        ax.plot(df.index, df['helix_percent'], color=colors['helix_percent'], label='Helix')

        # labels and axes
        ticks = ax.get_xticks().tolist()
        ax.set_xlim(ticks[1], ticks[-2])
        ax.set_ylim(0,100)

        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)

        ax.spines['right'].set_visible(False) # hide right axis
        ax.spines['top'].set_visible(False) # hide top axis
        ax.spines['right'].set_linewidth(10)
        ax.spines['top'].set_linewidth(10)


        # titling
        if rep is not None:
            if (isinstance(rep, int)):
                ax.set_title('Replicate {}'.format(rep), fontsize=24)
        if average_title is not None:
            ax.set_title(average_title, fontsize=24)


        ax.set_xlabel('Time (ns)', fontsize=20)
        ax.set_ylabel('Percentage (%)', fontsize=20)

        ax.xaxis.set_minor_locator(MultipleLocator(25000))

        # # legend

        # leg = ax.legend()
        # if 'legend.png' not in os.listdir(self.system.directory['images']['dssp']):
        #     self.legend(leg, 'dssp')
        # if not legend:
        #     leg.remove()

        # saving
        xvg_base = xvg.split('.')[0]
        if output is None:
            try:
                saveto = os.path.join(self.system.directory['images']['dssp'], '{}_{}.png'.format(xvg_base, rep))
            except:
                if (isistance(rep, int)):
                    saveto = os.path.join(os.getcwd(), '{}_{}.png'.format(xvg_base, rep))
                elif (isinstance(rep, str)):
                    saveto = os.path.join(os.getcwd(), '{}_average.png'.format(xvg_base))
                else:
                    saveto = 'dssp.png'
        else:
            saveto = output
        plt.savefig(saveto, dpi=600)
        plt.close()
        print('Plotted {}'.format(saveto))
        return saveto

    def legendIMG(self, legend, location):
        fig  = legend.figure
        fig.canvas.draw()
        bbox  = legend.get_window_extent()
        expand=[-5,-5,5,5]
        bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
        bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
        output = os.path.join(self.system.directory['images'][location], 'legend.png')
        fig.savefig(output, dpi=300, bbox_inches=bbox)

    @staticmethod
    def legend(data, handle_type):
        import matplotlib.lines as mlines
        import matplotlib.patches as mpatches 

        handles = []
        for label, color in data.values():
            if handle_type == 'line2d':
                handles.append(mlines.Line2D([], [], color=color, label=label))
        
        legend = plt.legend(handles=handles)
        return legend

    @staticmethod
    def autolabel(rects, ax, caplines=None, stdev=None, xpos='center', fontsize=None):
        """
        Attach a text label above each bar in *rects*, displaying its height.

        *xpos* indicates which side to place the text w.r.t. the center of
        the bar. It can be one of the following {'center', 'right', 'left'}.
        """

        xpos = xpos.lower()  # normalize the case of the parameter
        ha = {'center': 'center', 'right': 'left', 'left': 'right'}
        offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

        i = 0
        for rect in rects:
            height = rect.get_height()
            _height = rect.get_height()
            if caplines is not None:
                _height = caplines[i]
            if stdev is None:
                if _height > 0:
                    ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.02*_height,
                            '{}'.format(round(height,2)), ha=ha[xpos], va='bottom', fontsize=fontsize)
                else:
                    ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.02*_height,
                            '{}'.format(round(height,2)), ha=ha[xpos], va='top', fontsize=fontsize)

            else:
                std = stdev[i]
                txt = '{} '.format(round(height, 2)) + r'$\pm$' + ' {}'.format(round(stdev[i],2))
                ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.02*height,
                        txt, ha=ha[xpos], va='bottom', fontsize=fontsize)
            i += 1
        return ax


    def barPlotOld(self,csv,labels, output, title=None, mix=False, color=None, ab_color=None, iapp_color=None):
        colors = ['#000080']
        fig, ax = plt.subplots()
        if len(labels) == 10:
            fig.set_size_inches(7,3)
        elif len(labels) == 7:
            fig.set_size_inches(5,3)
        else:
            fig.set_size_inches(12,3)
        df = pd.read_csv(csv, index_col=0)
        print(df)
        x_pos = [i for i in range(1, len(df.index)+1)]
        if mix is True:
            ab_params = {
                'xpos':[],
                'height':[],
                'stdev':[]
            }
            iapp_params = {
                'xpos':[],
                'height':[],
                'stdev':[]
            }

            ab = True
            iapp = False
            for i in range(1, len(labels)+1):
                if labels[i-1] == 'SER20':
                    iapp = True
                    ab = False
                    iapp_params['xpos'].append(i)
                elif ab is True:
                    ab_params['xpos'].append(i)
                else:
                    iapp_params['xpos'].append(i)
                
            ab_params['height'] = df['hydrogen bonds'][:7]
            ab_params['stdev'] = df['stdev'][:7]
            iapp_params['height'] = df['hydrogen bonds'][7:]
            iapp_params['stdev'] = df['stdev'][7:]
            rects1 = ax.bar(ab_params['xpos'], ab_params['height'], yerr=ab_params['stdev'], alpha=0.5, ecolor='#5c5c5c', capsize=4, color=ab_color)
            rects2 = ax.bar(iapp_params['xpos'], iapp_params['height'], yerr=iapp_params['stdev'], alpha=0.5, ecolor='#5c5c5c', capsize=4, color=iapp_color)
            caplines = rects1.errorbar.lines[1]
            caplines = caplines[1]._y
            self.autolabel(rects1, ax, caplines)
            caplines = rects2.errorbar.lines[1]
            caplines = caplines[1]._y
            self.autolabel(rects2, ax, caplines)
        else:
            rects = ax.bar(x_pos, df['hydrogen bonds'], yerr=df['stdev'], alpha=0.5, ecolor='#5c5c5c', capsize=4, color=color)
            caplines = rects.errorbar.lines[1]
            caplines = caplines[1]._y
            self.autolabel(rects, ax, caplines)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels)
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.set_ylim(0,10)
        ax.set_ylabel('Average Number of Hydrogen Bonds')
        fontdict = {
            'fontsize':20
        }
        if title is not None:
            ax.set_title(title, fontdict)
        fig.tight_layout()
        plt.savefig(output, dpi=300)



class SystemPlot:

    def __init__(self, system=None):
        if system is not None:
            self.root = system.root
            self.reps = system.reps
            self.name = system.name
            self.source = system.source
            self.peptides = system.peptides
            self.ligands = system.ligands
            self.cascades = system.cascades
            self.directory = system.directory
            self.protein = system.protein
            self.system = system
        else:
            self.system = None

        self.plotter = Plotter(system=self.system)

    def plotRMSF(self, marker=True, lineplot=False):
        plotter = Plotter(self)
        plotter.rmsf(marker, lineplot)
    
    def mindist(self, group, title, average_title, **kwargs):
        matrices = []
        for rep in range(1, self.reps+1):
            matrix = self.system.post.mindist(dtype=group, rep=rep)
            pdata = PlotData.mindist(df=matrix, xlabels=xlabels, ylabels=ylabels, title=title)
            matrix = self.system.post.mindist(dtype=group, rep=rep)
            matrices.append(matrix)
            self.plotter.mindist(matrix=matrix, dtype=group, rep=rep, title=title)
            print('\n')
        # combined = pd.concat(matrices)
        # average = combined.groupby(level=0).mean()
        # self.plotter.mindist(matrix=average, dtype=group, rep='average', title=average_title)
        # print('\n\n')

    def mindistResidueSM(self, group, title, average_title, xlabels=None, ylabels=None, offset=False, rotation=False, colormap='plasma', normalize=True, ndx=None):
        dfs = []
        files = []
        for rep in range(1, self.system.reps+1):
            df = self.system.post.mindistResidueSM(group=group, rep=rep, normalize=normalize, ndx=ndx)
            rep_title = 'Replicate {}'.format(rep)
            output = os.path.join(self.system.directory['images']['mindist'], 'mindist_{}_{}.svg'.format(group, rep))
            pdata = PlotData.mindistResidueSM(df=df, colormap=colormap, title=rep_title, xlabels=xlabels, ylabels=ylabels, offset=offset, colorbar=True, vmin=0, vmax=1, name=self.system.name, output=output)
            fname = self.plotter.heatmap(pdata)
            dfs.append(df)
            files.append(fname)
            # sys.exit(0)
        combined = pd.concat(dfs)
        average = combined.groupby(level=0).mean()
        fname = '{}_mindist_residue_sm.csv'.format(self.system.name)
        average.to_csv(fname)
        if group == 'sm_sm':
            average = average.reindex(df.index)
            average = average.T.reindex(df.index).T
            average.to_csv(os.path.join(self.system.directory['images']['mindist'], 'average.csv'))
        output = os.path.join(self.system.directory['images']['mindist'], 'mindist_{}_average.svg'.format(group))
        pdata = PlotData.mindistResidueSM(df=average, colormap=colormap, title=average_title, xlabels=xlabels, ylabels=ylabels, offset=offset, colorbar=True, vmin=0, vmax=1, name=self.system.name, output=output)
        self.plotter.heatmap(pdata)

    def mindistOverTime(self, group, title, average_title, colors, block_average=100, legend=True, panel=True, width=2):
        dfs = []
        files = []
        for rep in range(1, self.system.reps+1):
            df = self.system.post.mindistOverTime(group=group, rep=rep, block_average=block_average)
            dfs.append(df)
            output = os.path.join(self.system.directory['images']['mindist'], 'mindist_over_time_{}_{}.png'.format(group, rep))
            pdata = PlotData.mindistOverTime(df=df, colors=colors, title='Replicate {}'.format(rep), legend=legend, output=output)
            fname = self.plotter.timeseries(pdata)
            files.append(fname)
        combined = pd.concat(dfs)
        average = combined.groupby(level=0).mean()
        output = os.path.join(self.system.directory['images']['mindist'], 'mindist_over_time_{}_average.png'.format(group))
        pdata = PlotData.mindistOverTime(df=average, colors=colors, title=average_title, legend=legend, output=output)
        self.plotter.timeseries(pdata)
        if panel is True:
            self.panel(files=files, width=width, output=os.path.join(self.system.directory['images']['mindist'], 'mindist_over_time_panel.png'), title=title, legend_data=None)

    def hbondsOverTime(self, group, title, average_title, colors, labels, block_average=100, ymax=2, panel=True, width=2):
        dfs = []
        pngs = []
        svgs = []
        for rep in range(1, self.system.reps+1):
            df_dict = self.system.post.hbondsOverTime(dtype=group, rep=rep, block_average=block_average)
            output = os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_{}.png'.format(group, rep))
            pdata = PlotData.hbondsOverTime(dfs=df_dict, colors=colors, title='Replicate {}'.format(rep), xlabels=labels, ymax=ymax, output=output, block_average=100)
            fname = self.plotter.markerPlot(pdata)
            pngs.append(fname)
            output = os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_{}.svg'.format(group, rep))
            pdata = PlotData.hbondsOverTime(dfs=df_dict, colors=colors, title='Replicate {}'.format(rep), xlabels=labels, ymax=ymax, output=output, block_average=100)
            fname = self.plotter.markerPlot(pdata)
            svgs.append(fname)
            dfs.append(df_dict)
        averages = {}
        for dic in dfs:
            for i in range(block_average, (len(df_dict.keys())*block_average)+1, block_average):
                if i not in averages.keys():
                    averages[str(i)] = []
                averages[str(i)].append(dic[str(i)])
        for t, dfs in averages.items():
            combined = pd.concat(dfs)
            average = combined.groupby(level=0).mean()
            averages[t] = average
        output = os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_average.{}'.format(group, 'png'))
        pdata = PlotData.hbondsOverTime(dfs=averages, colors=colors, title=average_title, xlabels=labels, output=output, block_average=100)
        self.plotter.markerPlot(pdata)
        output = os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_average.{}'.format(group, 'svg'))
        pdata = PlotData.hbondsOverTime(dfs=averages, colors=colors, title=average_title, xlabels=labels, output=output, block_average=100)
        self.plotter.markerPlot(pdata)

        if panel is True:
            self.panel(files=pngs, width=width, output=os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_panel.png'.format(group)), title=title, legend_data=None)
            # self.panel(files=svgs, width=width, output=os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_panel.svg'.format(group)), title=title, legend_data=None)

    def hbondsHeatmap(self, group, title, colormap='viridis', period=(400,600), panel=True, width=2):
        dfs = []
        files = []
        for rep in range(1, self.system.reps+1):
            df = self.system.post.hbondsHeatmap(dtype=group, period=period, rep=rep)
            dfs.append(df)
            rep_title = 'Replicate {}'.format(rep)
            output = os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_{}.png'.format(group, rep))
            pdata = PlotData.hbondsHeatmap(df=df, colormap=colormap, title=rep_title, name=self.system.name, output=output)
            fname  = self.plotter.heatmap(pdata=pdata)
            files.append(fname)
        combined = pd.concat(dfs)
        average = combined.groupby(level=0).mean()
        output = os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_average.png'.format(group))
        pdata = PlotData.hbondsHeatmap(df=average, colormap=colormap, title=title, name=self.system.name, output=output)
        self.plotter.heatmap(pdata)

        if panel is True:
            self.panel(files=files, width=width, output=os.path.join(self.system.directory['images']['hbonds'], 'hbonds_{}_panel.png'.format(group)), title=title)

    def dsspPerResidue(self, plot, xpm):
        dfs = []
        for rep in range(1, self.system.reps+1):
            xpm_file = os.path.join(self.system.directory[rep]['dssp']['root'], xpm)
            post = PostProcess(system=self.system)
            df = post.dsspPerResidue(gro=self.system.gro, xpm=xpm_file, rep=rep)
            dfs.append(df)
            plotter = Plotter(self.system)
            plotter.dsspPerResidue(df=df, plot=plot, xpm=xpm, rep=rep)
        combined = pd.concat(dfs)
        average = combined.groupby(level=0).mean()
        plotter.dsspPerResidue(df=average, plot=plot, xpm=xpm, average=True)

    def gyrationOverTime(self, group, title, color, block_average=100):
        dfs = []
        files = []
        for rep in range(1, self.system.reps+1):
            df = self.system.post.gyration(group=group, rep=rep, block_average=block_average)
            dfs.append(df)
            output = os.path.join(self.system.directory['images']['gyration'], 'gyration_{}_{}.png'.format(group, rep))
            pdata = PlotData.gyration(df=df, color=color, title='Replicate {}'.format(rep), label=self.system.name, output=output)
            fname = self.plotter.timeseries(pdata)
            files.append(fname)
        combined = pd.concat(dfs)
        average = combined.groupby(level=0).mean()
        output = os.path.join(self.system.directory['images']['gyration'], 'gyration_{}_average.png'.format(group))
        pdata = PlotData.gyration(df=average, color=color, label=system.name, title=title, output=output)
        self.plotter.timeseries(pdata)
    
    def timeseries(self, type, block_average=100, panel=True, width=2, legend=True, title=None, average_title=None, colors=None):
        pass 
        #TODO do this 

    def dsspOverTime(self, xvg, block_average=100, panel=True, width=2, legend=True, title=None, average_title=None):
        '''
        Plot DSSP over time for System object.
        Arguments:
        ** xvg (string): path to .xvg file to plot
        ** block_average (int, optional, default 100): number of frames to block average. If none, no block averaging will be performed. 
        ** panel (bool, optional, default True): if True and plotting from System object, will create panel image with each plot from replicate graphs
        ** width (int, optional, default 2): if panel is True, will create a panel with specified width.  
        ** legend (bool, default True): if true will show legend
        ** title (str, optional): title for panel image
        '''
        dfs = []
        files = []
        plotter = Plotter(self)
        for rep in range(1, self.system.reps+1):
            xvg_file = os.path.join(self.system.directory[rep]['dssp']['root'], xvg)

            post = PostProcess(system=self.system)
            df = post.dsspOverTime(xvg=xvg_file, block_average=block_average)
            dfs.append(df)

            _title = 'Replicate {}'.format(rep)

            xvg_base = xvg.split('.')[0]
            output = os.path.join(self.system.directory['images']['dssp'], '{}_{}.png'.format(xvg_base, rep))
            data = PlotData.dsspOverTime(df, title=_title, output=output)
            files.append(plotter.timeseries(data))

        combined = pd.concat(dfs)
        average = combined.groupby(level=0).mean()
        coil = []
        bsheet = []
        helix = []
        for i in average.index:
            coil.append(np.std([df['coil_percent'][i] for df in dfs]))
            bsheet.append(np.std([df['bsheet_percent'][i] for df in dfs]))
            helix.append(np.std([df['helix_percent'][i] for df in dfs]))
        average['coil_std'] = coil
        average['bsheet_std'] = bsheet
        average['helix_std'] = helix
        output = os.path.join(self.system.directory['images']['dssp'], '{}_average.png'.format(xvg_base))
        data = PlotData.dsspOverTime(average, title=average_title, output=output)
        plotter.timeseries(data)

        # if panel is True:
        #     colors = [('Coil','#11174b'), (r'$\beta$-Strand','#62bed2'), ('Helix','#3967a4')]
        #     self.panel(files=files, width=width, output=os.path.join(self.system.directory['images']['dssp'], 'dssp_panel.png'), title=title, legend_data=colors)

    def rmsdPerPeptide(self, panel=True, width=2, legend=True, title=None, average_title=None, colors=None):
        dfs = []
        data = []
        plotter = Plotter(self)
        for rep in range(1, self.system.reps+1):
            xvg_files = self.system.directory[rep]['rmsd']['data']
            post = PostProcess(system=self.system)
            df = post.rmsdPerPeptide(files=xvg_files, peptides=self.system.peptides)
            dfs.append(df)

            t = 'Replicate {}'.format(rep)

            output = os.path.join(self.system.directory['images']['rmsd'], 'rmsd_{}.png'.format(rep))
            pdata = PlotData.rmsdPerPeptide(df, colors=colors, title=t, output=output)
            # self.plotter.timeseries(pdata)
            data.append(pdata)

        combined = pd.concat(dfs)
        average = combined.groupby(level=0).mean()
        output = os.path.join(self.system.directory['images']['rmsd'], 'rmsd_average.png')
        pdata = PlotData.rmsdPerPeptide(average, title=average_title, output=output)
        plotter.timeseries(pdata)

        if panel is True:
            for d in data:
                d.xlabel.fontsize = 14
                d.ylabel.fontsize = 14
                d.title.fontsize = 18
                d.xticks.fontsize = 12
            colors = [('Peptide 1','#11174b'), ('Peptide 2','#3967a4'), ('Peptide 3','#62bed2')]
            self.plotter.timeseriesPanel(data=data, output=os.path.join(self.system.directory['images']['rmsd'], 'rmsd_panel.png'), title=title, legend_data=colors)
            self.plotter.timeseriesPanel(data=data, output=os.path.join(self.system.directory['images']['rmsd'], 'rmsd_panel.svg'), title=title, legend_data=colors)

  
    def panelNew(self, data, ptype, width, output, title, legend_data=None):
        h = int(int(len(data)) / int(width))
        w = width
        fig, axes = plt.subplots(h, w)
        fig.set_size_inches(10,12)
        i = 0
        for x in range(0, h):
            for y in range(0, w):
                print(axes[x,y])
                pdata = data[i]
                pdata.saveto = None
                axes[x,y] = self.plotter.timeseries(pdata)
                print(axes[x,y])
                i += 1

        fig.suptitle(title, y=0.90, fontsize=16)

        import matplotlib.lines as mlines
        if legend_data is not None:
            handles = []
            for label, color in legend_data:
                handles.append(mlines.Line2D([], [], color=color, label=label))
            
            plt.legend(handles=handles, bbox_to_anchor=(-0.81, -0.15, 1.6, .102), loc='lower left',
            ncol=3, mode="expand", borderaxespad=0.)

        plt.show()
        # if output.endswith('png'):
        #     plt.savefig(output, dpi=300)
        # else:
        #     plt.savefig(output)
        # print('Plotted {}'.format(output))
        plt.close()  


    @staticmethod
    def panel(files, width, output, title, legend_data=None):
        fig, axs = plt.subplots(int(len(files)/width), width)
        fig.set_size_inches(10,12)
        i = 0
        for row in axs:
            for ax in row:
                img = plt.imread(files[i])
                ax.imshow(img)
                ax.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    left=False,
                    labelbottom=False,
                    labelleft=False,)
                ax.axis('off')
                f = ax
                i += 1
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.suptitle(title, y=0.90, fontsize=16)

        # legend

        import matplotlib.lines as mlines
        if legend_data is not None:
            handles = []
            for label, color in legend_data:
                handles.append(mlines.Line2D([], [], color=color, label=label))
            
            plt.legend(handles=handles, bbox_to_anchor=(-0.81, -0.15, 1.6, .102), loc='lower left',
            ncol=3, mode="expand", borderaxespad=0.)
        
        if output.endswith('png'):
            plt.savefig(output, dpi=300)
        else:
            plt.savefig(output)
        print('Plotted {}'.format(output))
        plt.close()

    @staticmethod
    def panelMultiSystem(systems, png, width, output, title, legend_data=None):
        fig, axs = plt.subplots(int(len(systems)/width), width)
        fig.set_size_inches(10,12)
        i = 0
        for row in axs:
            for ax in row:
                sys = systems[i]
                location = png.split('_')[0]
                f = os.path.join(sys.directory['images'][location], png)
                img = plt.imread(f)
                ax.imshow(img)
                ax.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    left=False,
                    labelbottom=False,
                    labelleft=False,)
                ax.axis('off')
                f = ax
                i += 1
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.suptitle(title, y=0.90, fontsize=16)

        # legend

        import matplotlib.lines as mlines
        if legend_data is not None:
            handles = []
            for label, color in legend_data:
                handles.append(mlines.Line2D([], [], color=color, label=label))
            
            plt.legend(handles=handles, bbox_to_anchor=(-0.81, -0.15, 1.6, .102), loc='lower left',
            ncol=3, mode="expand", borderaxespad=0.)

        plt.savefig(output, dpi=300)
        print('Plotted {}'.format(output))
        plt.close()


# colors = {
#     'Coil':'#11174b',
#     r'$\beta$-Strand':'#62bed2',
#     'Helix':'#3967a4'
# }
# colors = [('Coil','#11174b'), (r'$\beta$-Strand','#62bed2'), ('Helix','#3967a4')]
# files = ['images/dssp_0_600_1.png', 'images/dssp_0_600_2.png', 'images/dssp_0_600_3.png', 'images/dssp_0_600_4.png', 'images/dssp_0_600_5.png', 'images/dssp_0_600_6.png']
# SystemPlot.panel(files, width=2, output='images/testtitle.png', title='IAPP$_{(20-29)}$ Secondary Structure', legend_data=colors)
        
        

class FreePlot:

    @staticmethod
    def dsspPerResidue(xpm, gro, plot, peptides=None):
        post = PostProcess(system=None)
        df = post.dsspPerResidue(gro=gro, xpm=xpm, rep=None)
        plotter = Plotter(system=None)
        plotter.dsspPerResidue(df=df, plot=plot, xpm=xpm, rep=None, peptides=peptides)

    @staticmethod
    def dsspOverTime(xvg, block_average=100, panel=True, width=2, rep=None, output=None, average=False, legend=True, title=None):
        '''
        Plot DSSP over time.
        Arguments:
        ** xvg (string, list): path (str) or paths (list) to .xvg file(s) to plot
        ** block_average (int, optional, default 100): number of frames to block average. If none, no block averaging will be performed. 
        ** panel (bool, optional, default True): if True and plotting from System object, will create panel image with each plot from replicate graphs
        ** width (int, optional, default 2): if panel is True, will create a panel with specified width.  
        ** rep (int, optional): replicate number (for titling)
        ** output (str, optional): output for single graph or panel image
        ** average (bool, default False): if xvg is list, create average plot
        ** legend (bool, default True): if True, will show legend on plot
        '''
        
        post = PostProcess(system=None)
        plotter = Plotter(system=None)
        if isinstance(xvg, str):
            df = post.dsspOverTime(xvg, block_average=block_average)
            plotter.dsspOverTime(df, xvg, rep, output, legend=legend)
        if isinstance(xvg, list):
            files = []
            dfs = []
            for x in xvg:
                df = post.dsspOverTime(xvg, block_average=block_average)
                dfs.append(df)
                files.append(plotter.dsspOverTime(df, xvg, rep, output, title=title))
            if average is True:
                combined = pd.concat(dfs)
                average = combined.groupby(level=0).mean()
                plotter.dsspOverTime(average, xvg, rep='average')
                plotter.dsspOverTime(average, xvg, rep='average', output=output, title=title)
            if panel is True:
                plotter.panel(files, width, output)


