import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import os
import pandas as pd
import sys

sys.path.append(os.getcwd())

from pymd.mdanalysis.postprocess import PostProcess


class Plotter:
    def __init__(self, system=None):
        self.system = system

    def rmsd(self, xvgs):
        for xvg in xvgs:
            df = pd.read_csv(xvg, delim_whitespace=True, skiprows=17, header=None)
            df.columns=['residue', 'rmsf']
        

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


    def mindist(self, matrix, group, rep=None, testing=False, kwargs=None):
        fig, ax = plt.subplots()
        fig.set_size_inches(8,6)
        img = ax.pcolor(matrix, vmin=0, vmax=30, cmap='plasma')
        cbar = fig.colorbar(img)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        if testing == False:
            if rep is not None:
                saveto = os.path.join(self.system.directory[rep]['mindist'][group], 'mindist_{}.png'.format(group))
            else:
                saveto = os.path.join(self.system.root, 'mindist_{}_average.png'.format(group))
            plt.savefig(saveto, dpi=300)
            plt.close()
            print('Plotted {}'.format(saveto))
        else:
            plt.show()
    
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
        print(xpm)
        if rep is None:
            xpm_base = xpm[:-3] + 'png'
            saveto = xpm_base
        elif average is False:
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
    




class SystemPlot:

    def __init__(self, system):
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

    def plotRMSF(self, marker=True, lineplot=False):
        plotter = Plotter(self)
        plotter.rmsf(marker, lineplot)
    
    def plotMindist(self, group, **kwargs):
        matrices = []
        for rep in range(1, self.reps+1):
            post = PostProcess(self)
            path = self.directory[rep]['mindist'][group]
            matrix = post.mindist(path)
            matrices.append(matrix)
            plotter = Plotter(self)
            plotter.mindist(matrix=matrix, group=group, rep=rep)
            print('\n')
        combined = pd.concat(matrices)
        average = combined.groupby(level=0).mean()
        plotter.mindist(matrix=average, group=group)
        print('\n\n')


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
        

class FreePlot:

    def dsspPerResidue(self, xpm, gro, plot, peptides=None):
        post = PostProcess(system=None)
        df = post.dsspPerResidue(gro=gro, xpm=xpm, rep=None)
        plotter = Plotter(system=None)
        plotter.dsspPerResidue(df=df, plot=plot, xpm=xpm, rep=None, peptides=peptides)