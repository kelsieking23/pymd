import os

import matplotlib
if os.name != 'nt':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
import sys

if os.name == 'nt':
    sys.path.append(os.getcwd())
    sys.path.append('D:/Work/iapp')
    from pymd.mdanalysis.postprocess import PostProcess
    from pymd.plot.data import PlotData
    from pymd.plot.plotter import Plotter
else:
    sys.path.append('/work/cascades/kelsieking23/iapp_analysis/scripts/python/')
    from mdanalysis.postprocess import PostProcess
    from plot.data import PlotData
    from plot.plotter import Plotter


class MultiSystemPlot:

    def __init__(self, systems):
        self.systems = systems
        self.plotter = Plotter()

        self.save = None
        self.last_ran = None

    def mindistResidueSM(self, group, title, ndxfiles, colormap, xlabels, ylabels, output):
        i = 0
        md = pd.DataFrame()
        for system in self.systems:
            dfs = []
            for rep in range(1, system.reps+1):
                print(ndxfiles[i])
                print(system.name)
                df = system.post.mindistResidueSM(group=group, rep=rep, normalize=True, ndx=ndxfiles[i])
                dfs.append(df)
            combined = pd.concat(dfs)
            average = combined.groupby(level=0).mean()
            md[average.iloc[0,:].name] = average.iloc[0,:]
            i += 1
        df = md.T

        pdata = PlotData.mindistResidueSM(df=df, colormap=colormap, title=title, xlabels=xlabels, ylabels=ylabels, offset=True, rotation=None, colorbar=True, output=output)
        self.plotter.heatmap(pdata)

    def gyration(self, group, title=None, colors=None, output=None, ecc=False, block_average=200, tabulate=(400,600)):
        data = {}
        names = []
        for system in self.systems:
            dfs = []
            names.append(system.name)
            if ecc is False:
                df = system.post.gyration(group=group, block_average=block_average, average_replicates=True)
            else:
                df = system.post.eccentricity(group=group, block_average=block_average, average_replicates=True)
            index = df['time']
            if system.alias is not None:
                data[system.alias] = df['gyr']
            else:
                data[system.name] = df['gyr']
        df = pd.DataFrame(data)
        df.index = index
        df = df.loc[400:]
        df.to_csv('gyr_all.csv')
        # if ecc is False:
        #     pdata = PlotData.gyrationMulti(df=df, title=title, colors=colors, output=output)
        # else:
        #     pdata = PlotData.gyrationMulti(df=df, title=title, ecc=True, colors=colors, output=output)
        # self.plotter.timeseries(pdata)

        # if tabulate is not None:
        #     data = {}
        #     for system in self.systems:
        #         dfs = []
        #         for rep in range(1, system.reps+1):
        #             if ecc is False:
        #                 df = system.post.gyration(group=group, rep=rep, block_average=1)
        #             else:
        #                 df = system.post.eccentricity(group=group, rep=rep, block_average=1)
        #             dfs.append(df)
        #         combined = pd.concat(dfs)
        #         average = combined.groupby(level=0).mean()
        #         data[system.name] = average['gyr']
        #     df = pd.DataFrame(data)
        #     df.index = average['time']
        #     data = {
        #         'average':{},
        #         'stdev':{}
        #     }
        #     for column in df.columns:
        #         mean = sum(df.loc[tabulate[0]:tabulate[1], column])/len(df.loc[tabulate[0]:tabulate[1], column])
        #         stdev = np.std(df.loc[tabulate[0]:tabulate[1], column])
        #         data['average'][column] = round(mean,1)
        #         data['stdev'][column] = round(stdev,1)
        #     df = pd.DataFrame(data)
        #     df.reindex(names)
        #     print(df)
    
    def hbondsBarPlot(self, group, title, output, period=(400,600), labels=None, color=None):
        hbonds = {
            'mean':{},
            'stdev':{}
        }
        names = []
        for system in self.systems:
            total_ = 0
            frames_ = 0
            means = []
            names.append(system.name)
            for rep in range(1, system.reps+1):
                total, frames = system.post.hbondsTotal(rep=rep, dtype=group, period=period)
                total_ += total
                frames_ += frames
                m = total / frames
                means.append(m)
            mean = total_ / frames
            stdev = np.std(means) 
            hbonds['mean'][system.name] = round(mean,2)
            hbonds['stdev'][system.name] = round(stdev,2)

        df = pd.DataFrame(hbonds).reindex(names)
        
        print(df)
        pdata = PlotData.hbondsTotalBarMulti(df=df, title=title, color=None, labels=None, output=output)
        self.plotter.barPlot(pdata)

    def mindistOverTime(self, group, title, output, labels=None, colors=None, tabulate=(400,600), block_average=100):
        values = {}
        first = True
        k = 0
        dfs = []
        for system in self.systems:
            mindist = pd.DataFrame()
            for rep in range(1, system.reps+1):
                df = system.post.mindistOverTime(group=group, rep=rep, block_average=block_average)
                mindist[rep] = df.mean(axis=1)
            # average = combined.groupby(level=0).mean().mean(axis=1) # average across peptides
            # mindist[system.name] = average
            mindist['mean'] = mindist.mean(axis=1)
            mindist['stdev'] = mindist.std
            df = pd.DataFrame()
            df['mean'] = mindist.mean(axis=1)
            df['stdev'] = mindist.std(axis=1)
            dfs.append(df)
            k += 1
        pdata = PlotData.mindistOverTimeMulti(dfs=dfs, colors=colors, title=title, legend=True, output=output, labels=labels)
        self.plotter.timeseries(pdata)

        if tabulate is not None:
            md = pd.DataFrame()
            i = 0
            for df in dfs:
                md[i] = df['mean']
                i += 1
            
            md = md.loc[tabulate[0]:tabulate[1], :]
            for column in md.columns:
                for col in md.columns:
                    if col == column:
                        continue
                    l1 = list(md[column])
                    l2 = list(md[col])
                    tstat, pval = stats.ttest_ind(l1,l2)
                    print('{} and {}: {}'.format(column, col, round(pval, 2)))
            md.to_csv('mindist_overtime_avg_ctrl_mut.csv')

    def dsspOverTime(self, xvg, suptitle, titles, output, block_average=200, tabulate=(400,600), colors=None, legend=True, labels=None):
        # dssp = pd.DataFrame()
        data = []
        k = 0
        ala = pd.DataFrame()
        alb = pd.DataFrame()
        for system in self.systems:
            dfs = []
            print(system.name, system.alias)
            for rep in range(1, system.reps+1):
                if system.name != 'Dihydroquercetin':
                    xvg_file = os.path.join(system.directory[rep]['dssp']['root'], xvg) 
                else:
                    xvg_file = os.path.join(system.directory[rep]['dssp']['root'], 'dssp.xvg') 
                df = system.post.dsspOverTime(xvg=xvg_file, block_average=block_average)
                dfs.append(df)
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
            alb[system.alias] = average['bsheet_percent'].loc[400:]
            ala[system.alias] = average['coil_percent'].loc[400:]
            pdata = PlotData.dsspOverTime(df=average, title=titles[k], output=None)
            data.append(pdata)
            k += 1
        for d in data:
            d.xlabel.fontsize = 14
            d.ylabel.fontsize = 14
            d.title.fontsize = 18
            d.xticks.fontsize = 12
        colors = {
            'Coil':'#11174b',
            r'$\beta$-Strand':'#62bed2',
            'Helix':'#3967a4'
        }
        colors = [('Coil', '#11174b'), (r'$\beta$-Strand','#62bed2'),('Helix','#3967a4')]
        # self.plotter.timeseriesPanel(data=data, output=output, title=suptitle, legend_data=colors)
        ala.to_csv('coil_avgs_wmut.csv')
        alb.to_csv('bsheet_avgs_wmut.csv')


        
        # if labels is None:
        #     labels = list(dssp.columns)
        # pdata = PlotData.dsspOverTimeMulti(df=dssp, title=title, colors=colors, legend=True, labels=labels, output=output)
        # self.plotter.timeseries(pdata)

        # if tabulate is not None:
        #     pairs = []
        #     ds = dssp.loc[tabulate[0]:tabulate[1], :]
        #     mean = ds.mean(axis=0)
        #     stdev = ds.std(axis=0)
        #     df = pd.DataFrame({'mean':mean, 'stdev':stdev}).round(2)
        #     print(df)
        #     ds.to_csv('average_ss.csv')
        #     for column in ds.columns:
        #         for col in ds.columns:
        #             if col == column:
        #                 continue
        #             if ((col, column) in pairs) or ((column, col) in pairs):
        #                 continue
        #             l1 = list(ds[column])
        #             l2 = list(ds[col])
        #             tstat, pval = stats.ttest_ind(l1,l2)
        #             if pval < 0.05:
        #                 print('{} and {}: {}'.format(column, col, pval))
        #             pairs.append((column, col))
        #             pairs.append((col, column))
    
    def contactDistribution(self, title, output, colors=None, legend=True, labels=None):
        polar = {
            'mean':{},
            'stdev':{}
        }
        hydrophobic = {
            'mean':{},
            'stdev':{}
        }
        totals_polar = {}
        totals_hydrophobic = {}
        names = []
        for system in self.systems:
            ndxt = '{}_ndx_template.ndxt'.format(system.ligands.lower())
            df, temp = system.post.contactTypes(ndxt=ndxt)
            polar['mean'][system.name] = df['polar']['mean']
            polar['stdev'][system.name] = df['polar']['stdev']
            hydrophobic['mean'][system.name] = df['hydrophobic']['mean']
            hydrophobic['stdev'][system.name] = df['hydrophobic']['stdev']
            names.append(system.name)
            totals_hydrophobic[system.name] = temp['hydrophobic']
            totals_polar[system.name] = temp['polar']
        polar = pd.DataFrame(polar).round(2)
        polar = polar.reindex(names)
        polar.name = 'Polar'
        hydrophobic = pd.DataFrame(hydrophobic).round(2)
        hydrophobic = hydrophobic.reindex(names)
        hydrophobic.name = 'Hydrophobic'
        pdata = PlotData.contactDistributionMulti(dfs=[polar, hydrophobic], title=title, colors=colors, labels=labels, output=output)
        self.plotter.barPlot(pdata)

        pairs = []
        print('Polar T-Test')
        print('******')
        for key, value in totals_polar.items():
            for k, v in totals_polar.items():
                if k == key:
                    continue
                if ((key,k) in pairs) or ((k,key) in pairs):
                    continue
                tstat, pval = stats.ttest_ind(value, v)
                if pval < 0.05:
                    print('{} and {}: {}'.format(key,k,pval))
                pairs.append((key,k))
                pairs.append((k,key))
        pairs = []
        print('\n')
        print('Hydrophobic T-Test')
        print('******')
        for key, value in totals_hydrophobic.items():
            for k, v in totals_hydrophobic.items():
                if k == key:
                    continue
                if ((key,k) in pairs) or ((k,key) in pairs):
                    continue
                tstat, pval = stats.ttest_ind(value, v)
                if pval < 0.05:
                    print('{} and {}: {}'.format(key,k,pval))
                pairs.append((key,k))
                pairs.append((k,key))
                
    def mindist(self, group, suptitle, titles, output, xlabels=None, ylabels=None, colormap='plasma'):
        data = []
        k = 0
        for system in self.systems:
            print(system.name)
            print(titles[k])
            dfs = system.post.mindist(dtype=group, normalize=True, average_peptides=True, testing=False)
            combined = pd.concat(dfs)
            average = combined.groupby(level=0).mean()
            print(xlabels)
            if isinstance(xlabels[0], str):
                xl = xlabels
            else:
                xl = xlabels[k]
            if isinstance(ylabels[0], str):
                yl = ylabels
            else:
                yl = ylabels[k]
            pdata = PlotData.mindist(df=average, xlabels=xl, ylabels=yl, colormap=colormap, legend=False, 
                    title=titles[k], output=None)
            data.append(pdata)
            k += 1
        self.plotter.heatmapPanel(data=data, output=output, title=suptitle, legend_data=None)
        
    def hbondsOverTime(self, group, suptitle, output, titles, xlabels, colors, block_average=100):
        data = []
        k = 0
        for system in self.systems:
            dfs = system.post.hbondsOverTime(dtype=group, block_average=100)
            averages = {}
            for dic in dfs:
                for i in range(block_average, (len(dic.keys())*block_average)+1, block_average):
                    if i not in averages.keys():
                        averages[str(i)] = []
                    averages[str(i)].append(dic[str(i)])
            for t, dfs in averages.items():
                combined = pd.concat(dfs)
                average = combined.groupby(level=0).mean()
                averages[t] = average
            pdata = PlotData.hbondsOverTime(dfs=averages, colors=colors, title=titles[k], xlabels=xlabels)
            data.append(pdata)
            k += 1
        self.plotter.markerPlotPanel(data, output=output, title=suptitle)

    def dsspPerResidue(self, xpm, suptitle, output, titles, xlabels, colors):
        data = []
        i = 0
        for system in self.systems:
            dfs = system.post.dsspPerResidueAverage(gro=system.gro, xpm=xpm, peptides=system.peptides)
            pdata = PlotData.dsspPerResidue(dfs=dfs, colors=colors, title=titles[i], xlabels=xlabels, output=None)
            i += 1
            data.append(pdata)
        legend_data = []
        i = 0
        s = 20
        for df in dfs:
            l = (df.name, colors[i], s)
            legend_data.append(l)
            i += 1
            s -= 3
        self.plotter.markerPlotPanel(data, output=output, title=suptitle, legend_data=legend_data)

    def rmsf(self, group, title, labels, output, color=None, start=None, stop=None):
        data = {}
        for system in self.systems:
            data[system.name] = {}
            df = pd.DataFrame()
            if group == 'sm':
                grp = system.ligands
            else:
                grp = group
            dfs = system.post.rmsf(grp, start=400, stop=600)
            means = []
            for d in dfs:
                mean = d['rmsf'].mean(axis=0)
                means.append(mean)
            mean = sum(means) / len(means)
            stdev = np.std(means)
            data[system.name]['mean'] = mean
            data[system.name]['stdev'] = stdev
        df = pd.DataFrame(data).T
        print(df)
        pdata = PlotData.rmsfBarPlot(df=df, title=title, color=color, labels=labels, output=output)
        self.plotter.barPlot(pdata)
            
    def linearRegressionResBsheet(self, xvg, group='sidechain', res1='PHE23', res2='PHE23', period=(400,600)):
        
        skip = False
        ids = self.systems[0].protein.ids
        y = []
        for res_id in ids:
            x = []
            for system in self.systems:
                if skip is False:
                # get bsheet
                    if system.name != 'Dihydroquercetin':
                        xvg_file = xvg
                    else:
                        xvg_file = 'dssp.xvg'
                    dssp = system.post.dsspOverTime(xvg=xvg_file, block_average=1, average_replicates=True)
                    bsheet = dssp['bsheet_percent'].loc[period[0]:period[1]].mean()
                    y.append(bsheet)

                # get res int freq
                freq = system.post.residueContactFrequency(group=group, res1=res_id, res2=res_id)
                x.append(freq)
            skip = True
        
            x = np.array(x).reshape((-1,1))
            y = np.array(y)

            model = LinearRegression().fit(x,y)
            r_sq = model.score(x,y)


            yield res_id, r_sq
        

            
        

            
                
            
    


            
            

                








        