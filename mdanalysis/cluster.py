from cProfile import run
import pandas as pd
import os
import sys
import threading
import time
import mdtraj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.metrics import pairwise_distances_argmin_min
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestCentroid
import multiprocessing as mp
import gc
import mdtraj
import warnings
warnings.filterwarnings('ignore')
# sys.path.append('D:/Work/')
from mdanalysis.rmsd_old import RMSD
# from pymd.mdanalysis.ndx import Ndx
from pymd.utilities.rewritepdb import writePDB
from pymd.structure.protein import Protein
from datetime import datetime
import random

class Cluster:

    def __init__(self, xtc, top, root=None, job_name='cluster', system_name='system', selection='backbone', stride=0, b=0, e=-1, cutoff=0.2, n_clusters=None, verbose=False, job_params = {}):
        self.xtc = xtc
        self.top = top
        self.selection = selection
        self.n_clusters = n_clusters
        self.clusters = []
        self.cutoff = cutoff
        self._traj = None
        self.df = None
        self._X = None
        self.X = None
        self.y = None
        self.job_name = job_name
        self.stride = stride
        self.b = b
        self.e = e
        self.analyze = None
        self._parent = None
        self._root = root
        self.argmins = None
        self.centroids = None
        self._job_params = job_params
        self._log = None
        self.verbose = verbose
        self.method = None
        self.fig = None
        self.ax = None
        self._root = None
        self.system_name = system_name

    @staticmethod
    def now():
        return datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    @property
    def parent(self):
        return self._parent
    
    @parent.setter
    def parent(self, parent):
        self._parent = parent
    
    @property
    def job_params(self):
        return self._job_params
    
    @job_params.setter
    def job_params(self, params):
        self._job_params = params

    @property
    def root(self):
        if self.parent is not None:
            root = os.path.join(self.parent.root, self.job_name)
        else:
            if os.path.dirname(self.xtc) == '':
                root =  os.path.join(os.getcwd(), self.job_name)
            else:
                root = os.path.join(os.path.dirname(self.xtc), self.job_name)
        if not os.path.isdir(root):
            os.mkdir(root)
        return root
    
    @root.setter
    def root(self, root):
        if not os.path.isdir(root):
            os.mkdir(root)
        self._root = root
    
    @property
    def log(self):
        return os.path.join(self.root, 'cluster.log')
    

    def load(self, filename):
        df = pd.read_csv(filename, index_col = 0, header = 0)
        self.df = df
        # self.analyze = ClusterData(df, self.traj, self.n_clusters)

    def loadTrajectory(self, stride=None):
        print('Loading {}'.format(self.xtc))
        print('Using {} as topology'.format(self.top))
        if stride is not None:
            if stride != 0:
                traj = mdtraj.load(self.xtc, top=self.top, stride=stride)
            else:
                traj = mdtraj.load(self.xtc, top=self.top)
        else:
            if self.stride != 0:
                traj = mdtraj.load(self.xtc, top=self.top, stride=self.stride)
            else:
                traj = mdtraj.load(self.xtc, top=self.top)
        traj = traj.superpose(traj)

        self._traj = traj
        sele = traj.top.select(self.selection)
        traj = traj.atom_slice(sele)
        if (self.e == -1):
            self.traj = traj.center_coordinates()[self.b:]
        else:
            self.traj = traj.center_coordinates()[self.b:self.e]
        self.frames = self.traj._xyz
        print(self._traj, self._traj._xyz.shape)

    def getPartitions(self, nprocs):
        if nprocs == 'auto':
            nprocs = mp.cpu_count() / 2
        nframes, _, _ = self.traj._xyz.shape
        print(self.traj.xyz.shape)
        interval = int(nframes // nprocs)
        print(nframes, nprocs)
        partitions = []
        procid=1
        for i in range(0, nframes, interval):
            data = {
                'b':i,
                'e':i+interval,
                'procid':procid,
            }
            partitions.append(data)
            procid+=1
            if ((i + interval + interval) > nframes) and (i+interval != nframes):
                data = {
                    'b':i+interval,
                    'e':nframes,
                    'procid':procid
                }
                partitions.append(data)
                break
        return partitions
    
    def run(self, nprocs, output=None, matrix=None):
        log = open(self.log, 'w')
        log.write('Created on: {}\n'.format(self.now()))
        log.write('Created by: pymd.mdanalysis.cluster\n\n')
        log.write('Cluster\n')
        log.write('**********\n')
        log.write('Job Parameters:\n')
        log.write('Trajectory: {}\n'.format(self.xtc))
        log.write('Topology: {}\n'.format(self.top))
        log.write('Selection: {}\n'.format(self.job_params['selection']))
        log.write('Start Index: {}\n'.format(self.job_params['b']))
        log.write('Stop Index: {}\n'.format(self.job_params['e']))
        log.write('Stride: {}\n'.format(self.job_params['stride']))
        if self.job_params['n_clusters'] is None:
            log.write('N Clusters: None (will determine via distortions)\n')
        else:
            log.write('N Clusters: {}\n'.format(self.job_params['n_clusters']))
        log.write('Method: {}\n'.format(self.job_params['method']))
        if self.job_params['method'] == 'agglomerative':
            log.write('Linkage: {}\n'.format(self.job_params['linkage']))
            if self.job_params['distance_threshold'] is not None:
                log.write('Distance Threshold: {}\n'.format(self.job_params['distance_threshold']))
            else:
                log.write('Distance Threshold: Not specified, will go by N-clusters\n')
        log.write('Dimensionality Reduction: {}\n'.format(self.job_params['reduce']))
        log.write('Scale: {}\n'.format(self.job_params['scale']))
        log.write('**********\n\n')
        log.write('Loading trajectory...\n')
        if self.verbose:
            print('Loading trajectory...')
        self.loadTrajectory()
        if matrix is None:
            log.write('Beginning RMSD matrix calculation for {}...\n'.format(self.system_name))
            if self.verbose:
                print('Beginning RMSD matrix calculation for {}...'.format(self.system_name))
            log.write('Performing {} processes\n'.format(nprocs))
            partitions = self.getPartitions(nprocs)
            # for partition in partitions:
            #     df = self.gromosV2(partition)
            try:
                pool = mp.Pool(processes=nprocs)
                results = pool.map(self.gromos, partitions)
            except Exception as e:
                log.write('Process failed.\n{}\n'.format(e))
                log.write('Terminating.')
                log.close()
                print(e)
                sys.exit(1)
            log.write('{} processes complete\n'.format(nprocs))
            if self.verbose:
                print('{} processes complete for {}'.format(nprocs, self.system_name))
            log.write('Concatenating data...\n')
            df = pd.concat(results, axis=1)
            df = df.fillna(df.T)
            log.write('Concatenation complete\n')
        else:
            log.write('Loading RMSD matrix...\nFilepath {}\n'.format(matrix))
            if self.verbose:
                print('Loading RMSD matrix...\nFilepath {}'.format(matrix))
            df = pd.read_csv(matrix, index_col=0)
            log.write('Loaded RMSD matrix successfully\n')
            if self.verbose:
                print('Loaded RMSD matrix successfully')
        # print(df)
        self.df = df
        # self.analyze = ClusterData(parent=None, dataframe=self.df, n_clusters=None)
        if output is not None:
            log.write('Writing .csv to {}...\n'.format(output))
            # print('writing csv to {}'.format(output))
            df.to_csv(output)
            log.write('Writing .csv complete\n')
        if matrix is None:
            log.write('RMSD matrix calculation complete. End time {}\n\n'.format(self.now()))
        log.close()
        # print('Done')
        return df

    def gromos(self, data):
        start_index = data['b']
        stop_index = data['e']
        procid = data['procid']
        # print('Starting procid {} - start {} stop {}'.format(procid, start_index, stop_index))
        df = pd.DataFrame()
        i = 1
        rms = RMSD(traj = self.traj, cutoff=self.cutoff, binary=False)
        nframes, _, _ = self.traj._xyz.shape
        # print(nframes)
        # for reference_index in range(start_index, stop_index): 
        # for reference_index in range(len(frames)):
        for reference_index in range(start_index, stop_index): 
            reference = self.frames[reference_index]
            rms.reference = reference
            rmsd = rms.rmsdOverTimeOLD(frames=self.frames[reference_index:])
            fill = [None]*(reference_index)
            if fill != []:
                df[reference_index] = fill + rmsd
            else:
                df[reference_index] = rmsd
            i += 1
        stop = time.time()
        # df.index = [i for i in range(slice[0], interval[1]-1)]
        # df.columns =  [i for i in range(interval[0], interval[1])]
        # print('Completed procid {}'.format(procid))
        # print('******** Procid {} ********'.format(procid))
        # print(df)
        # print('******** /end Procid {} ********'.format(procid))
        df.index = self.traj.time
        df.columns = self.traj.time[start_index:stop_index]
        self.df = df
        return df

    def pca(self, n_components=2, scale=True):
        if 'cluster.log' not in os.listdir(self.root):
            log = open(self.log, 'w')
        else:
            log = open(self.log, 'a')
        log.write('Performing PCA...\n')
        if self.verbose:
            print('Performing PCA for {}'.format(self.system_name))
        p = PCA(n_components=n_components)
        if scale is True:
            scaler = StandardScaler()
            self._X = scaler.fit(self.df).transform(self.df)
        self._X = p.fit_transform(self._X)
        self.X = self._X
        self.df = pd.DataFrame(self.X)
        self.df.to_csv(os.path.join(self.root, 'pca.csv'))
        log.write('PCA complete\nWrote {}\n\n'.format(os.path.join(self.root, 'pca.csv')))
        if self.verbose:
            print('PCA complete for {}'.format(self.system_name))
        log.close()
        return self.df

    def kmeans(self, n_clusters=None, output=None, title='K-Means', init='k-means++'):
        if 'cluster.log' not in os.listdir(self.root):
            log = open(self.log, 'w')
        else:
            log = open(self.log, 'a')
        log.write('Performing K-Means Clustering...\n')
        if self.verbose:
            print('Performing K-means for {}'.format(self.system_name))
        if self.n_clusters is not None:
            n_clusters = self.n_clusters
        elif n_clusters is not None:
            if self.n_clusters is None:
                self.n_clusters = n_clusters
        else:
            log.write('Number of clusters not defined\n')
            log.write('\tFinding K via distortions...\n')
            n_clusters = self.findK()
            self.n_clusters = n_clusters
            log.write('\tFinding K via distortions complete.\n')
        log.write('Number of clusters: {}\n'.format(n_clusters))
        if self.verbose:
            print('N clusters: {}'.format(n_clusters))
        km = KMeans(
            n_clusters=n_clusters, init=init,
            n_init=10, max_iter=300,
            tol=1e-04, random_state=0
        )
        self.y = km.fit_predict(self.X)
        self.centroids = km.cluster_centers_
        if output is not None:
            self.plotClusters(title='{}\nK = {}'.format(self.system_name, n_clusters), output=output)
        else:
            output = os.path.join(self.root, 'kmeans.png')
            self.plotClusters(title='{}\nK = {}'.format(self.system_name, n_clusters), output=output)
        self.argmins, _ = pairwise_distances_argmin_min(km.cluster_centers_, self.X)
        log.write('K-means clustering complete.\n\n')
        if self.verbose:
            print('K-means clustering for {} complete'.format(self.system_name))
        log.close()
        self.method = 'kmeans'

    def agglomerative(self, n_clusters=None, linkage='ward', distance_threshold=None):
        # if self.n_clusters is None:
        #     if n_clusters is not None:
        #         self.n_clusters = n_clusters
        wa = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage, distance_threshold=distance_threshold)
        self.y = wa.fit_predict(self.X)
        self.n_clusters = wa.n_clusters_
        self.centroids = self.getCentroids()
        if n_clusters is not None:
            output = os.path.join(self.root, 'agglomerative.n{}.{}.png'.format(n_clusters, linkage))
        else:
            output = os.path.join(self.root, 'agglomerative.dt{}.{}.png'.format(distance_threshold, linkage))
        if n_clusters is None:
            n = len(self.centroids.keys())
        else:
            n = n_clusters
        print('N Clusters: {}'.format(n))
        self.plotClusters(output=output, title='Agglomerative Clustering - {} Linkage\nN = {}'.format(linkage.capitalize(), n), legend=False)
        self.method = 'agglomerative'
        # print(self.n_clusters)
        # plt.scatter(self.X[:,0],self.X[:,1], c=wa.labels_, cmap='rainbow')
        # plt.savefig(output)


    def dbscan(self, eps=0.5, min_samples=5, metric='euclidean'):
        d = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
        self.X = self.df.to_numpy()
        self.y = d.fit_predict(self.X)
        self.centroids = self.getCentroids()
        output = os.path.join(self.root, 'dbscan.png')
        title='{}\nN = {}'.format(self.parent.name, len(self.centroids.keys()))
        output = os.path.join(self.root, 'dbscan.png')
        self.plotClusters(title=title, output=output)
        self.method = 'dbscan'

    def gromosCluster(self, filepath=None, prefix=None, cutoff=0.15):
        if prefix is not None:
            filename = os.path.join(self.root, '{}.cluster.csv'.format(prefix))
        elif filepath is not None:
            filename = filepath
        else:
            filename = os.path.join(self.root, 'cluster.csv')
        df = pd.read_csv(filename, index_col=0, header=0)
        self.X = df
        central_clusters = []
        n_clusters = []
        last_df = None
        i = -1
        while not df.empty:
            i += 1
            df, _max_index, _max = self.neighbor_search(df, cutoff)
            central_clusters.append(_max_index)
            n_clusters.append(_max)
            last_df = df
        self.method = 'gromos'
        self.y = [central_clusters, n_clusters]
        return central_clusters, n_clusters
    
    def cutoffSearch(self, filepath=None, prefix=None):
        if prefix is not None:
            filename = os.path.join(self.root, '{}.cluster.csv'.format(prefix))
        elif filepath is not None:
            filename = filepath
        else:
            filename = os.path.join(self.root, 'cluster.csv')
        cutoffs = []
        n_clusters = []
        for cutoff in np.arange(0.06, 0.22, 0.02):
            central_clusters, _ = self.gromosCluster(filepath=filename, cutoff=cutoff)
            cutoffs.append(cutoff)
            n_clusters.append(len(central_clusters))
        plt.scatter(cutoffs, n_clusters)
        plt.plot(cutoffs, n_clusters)
        plt.xlabel('RMSD Cutoff', fontweight='bold', fontsize=16)
        plt.ylabel('Number of Clusters', fontweight='bold', fontsize=16)
        plt.xticks(np.arange(min(cutoffs), max(cutoffs), 0.02), fontsize=14)
        # plt.yticks(np.arange(min(n_clusters), max(n_clusters)+1, 2), fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()



    def neighbor_search(self, df, cutoff):
        _max = None
        _max_index = None
        n_neighbors = {}
        neighbors = {}
        df.columns = df.index
        for column in df.columns:
            d = df[column]
           
            n_neighbors[column] = len((d[d <= cutoff]))
            neighbors[column] = d[d<=cutoff]
            if _max is None:
                _max = len((d[d <= cutoff]))
                _max_index = column
            else:
                if _max < len((d[d <= cutoff])):
                    _max = len((d[d <= cutoff]))
                    _max_index = column
                elif _max == len((d[d <= cutoff])):
                    if random.choice([0,1]) == 0:
                        _max = len((d[d <= cutoff]))
                        _max_index = column
                else:
                    pass
        df = df.drop(neighbors[_max_index].index, axis=0)
        df = df.drop(neighbors[_max_index].index, axis=1)
        return df, _max_index, _max

    def getCentroids(self):
        # clf = NearestCentroid()
        # clf.fit(self.X, self.y)
        # self.centroids = clf.centroids_ 
        # return self.centroids
        centroids = {}
        for i in np.unique(self.y):
            x = self.X[self.y == i, 0]
            y = self.X[self.y == i, 1]
            mean_x = x.mean()
            mean_y = y.mean()
            centroids[i] = np.array([mean_x, mean_y])
        self.centroids = centroids
        return centroids

    def findK(self):
        distortions = []
        for i in range(1, 21):
            print('Finding disortions for {}...'.format(i))
            km = KMeans(
                n_clusters=i, init='k-means++',
                n_init=10, max_iter=300,
                tol=1e-04, random_state=0
            )
            km.fit(self.X)
            distortions.append(km.inertia_)
        plt.plot(range(1, 21), distortions, marker='o', color='tab:blue')
        plt.xticks([i for i in range(1,21)])
        plt.title('Distortion', weight='bold', fontsize=22)
        plt.xlabel('Number of clusters', fontsize=14, weight='bold')
        plt.ylabel('Distortion', fontsize=14, weight='bold')
        if not sys.platform == 'linux':
            plt.show()
        else:
            plt.savefig('findk.png')
            print('Since using agg, plotted to findk.png')
        plt.close()
        n_clusters = input('Enter number of clusters: ') 
        n_clusters = int(n_clusters)
        self.n_clusters = n_clusters
        return n_clusters

    def plotClusters(self, title='Cluster', output=None, legend=True):
        # print(self.X)
        # print(self.y)
        color = plt.cm.tab20(np.linspace(0, 1,20))
        matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler('color', color) 
        fig, ax = plt.subplots()
        fig.set_size_inches(8,6)
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
        # for i in range(0,self.n_clusters):
        #     plt.scatter(self.X[self.y == i, 0], self.X[self.y == i, 1],
        #             s=50, edgecolor='black',
        #             label='Cluster {}'.format(i+1),
        #             c=colors[i], alpha=0.3)
        df = pd.DataFrame(self.X)
        df['labels'] = self.y
        for i in np.unique(self.y):
            print('Cluster {}'.format(i))
            d = df[df['labels'] == i]
            ax.scatter(d[0], d[1],
                    s=50, edgecolor='black',
                    label='Cluster {}'.format(i+1), alpha=0.3)
        
        # plot the centroids
        if self.centroids is not None:
            if not isinstance(self.centroids, dict):
                plt.scatter(
                    self.centroids[:, 0], self.centroids[:, 1],
                    s=250, marker='*',
                    c='yellow', edgecolor='black',
                    label='Centroids'
                )
            else:
                for key in self.centroids.keys():
                    plt.scatter(self.centroids[key][0], self.centroids[key][1],
                    s=250, marker='*', c='yellow', edgecolor='black'
                )
        # if legend:
        #     plt.legend(scatterpoints=1,fontsize=10)
        plt.title(title, weight='bold', fontsize=22)
        if not sys.platform == 'linux':
            if output is not None:
                if 'cluster.log' not in os.listdir(self.root):
                    log = open(self.log, 'w')
                else:
                    log = open(self.log, 'a')
                log.write('Plotting {}\n'.format(output))
                log.close()
                plt.savefig(output)
            # plt.show()
        else:
            if hasattr(self, 'prefix'):
                plt.savefig('{}_kmeans.png'.format(prefix))
            else:
                plt.savefig(output)
        self.fig = fig
        self.ax = ax
        # self.parent.fig = fig
        # self.parent.ax = ax
        plt.close()

    def getCentralStructureIndices(self, df):
        indices = []
        sizes = []
        if self.method != 'gromos':
            for i in self.sortByClusterSize():
                d = df[df['labels'] == i]
                centroid = self.centroids[i]
                _min = None
                _centroid = None
                for index in d.index:
                    point = np.array([d.loc[index,0], d.loc[index,1]])
                    dist = np.linalg.norm(point-centroid)
                    if (_min is None) or (dist < _min):
                        _min = dist
                        _centroid = index
                indices.append(_centroid)
                sizes.append(len(d))
            return zip(indices, sizes)
        else:
            print(self.y[0])
            print(self.y[1])
            times = self.y[0]
            df['index'] = df.index
            _times = [int(list(np.where(df['time'] == i))[0]) for i in times]
            return zip(_times, self.y[1])


    def writePDB(self, traj=None, top=None, prefix=None):
        # if self.argmins is None:
        #     self.kmeans()
        if self._traj is None:
            self.loadTrajectory()
            traj = self._traj
        else:
            traj = self._traj
        if 'cluster.log' not in os.listdir(self.root):
            log = open(os.path.join(self.root, 'cluster.log'), 'w')
        else:
            log = open(os.path.join(self.root, 'cluster.log'), 'a')
        if self.method != 'gromos':
            if not isinstance(self.centroids, dict):
                if self.centroids is None:
                    self.getCentroids()
                else:
                    centroids = {}
                    i = 0
                    for entry in self.centroids:
                        index = np.unique(self.y)[i]
                        centroids[index] = entry
                        i += 1
                    self.centroids = centroids
            df = pd.DataFrame(self.X)
            df['labels'] = self.y
        else:
            df = self.df
            _ind = df.index
            df.index = [i for i in range(0, len(df))]
            df.columns =[i for i in range(0, len(df))]
            df['time'] = _ind
        log.write('Cluster data:\n')
        if prefix is not None:
            prefix_log = open(os.path.join(self.root, '{}.cluster.log'.format(prefix)), 'w')
        size_data = []
        k = 0
        # _centroid is index in trajectory, d is # of members in cluster
        for _centroid, d in self.getCentralStructureIndices(df):
            frame = traj._xyz[_centroid]
            chain_index = 0
            chain_id = 'A'
            contents = ['REMARK t={}\n'.format(traj.time[_centroid])]
            contents.append('REMARK cluster size = {}\n'.format(d))
            contents.append('REMARK clust % = {}\n'.format((d/traj.n_frames)*100))
            size_data.append(['{}'.format(str(k).zfill(2)), 
                            traj.time[_centroid],
                            d,
                            (d/traj.n_frames)*100])
            if prefix is None:
                log.write('Cluster {}:\n'.format(str(k).zfill(2)))
                log.write('**********\n')
                log.write('t = {}\n'.format(traj.time[_centroid]))
                log.write('cluster size = {}\n'.format(d))
                log.write('clust % = {}\n'.format((d/traj.n_frames)*100))
            else:
                # print('Cluster {}:'.format(str(i).zfill(2)))
                # print('**********')
                prefix_log.write('Cluster {}:\n'.format(str(k).zfill(2)))
                prefix_log.write('**********\n')
                prefix_log.write('t = {}\n'.format(traj.time[_centroid]))
                # print('t = {}'.format(traj.time[_centroid]))
                prefix_log.write('cluster size = {}\n'.format(d))
                # print('cluster size = {}'.format(len(d)))
                prefix_log.write('clust % = {}\n'.format((d/traj.n_frames)*100))
                # print('clust % = {}\n'.format((len(d)/traj.n_frames)*100))
            for z in range(0, len(frame)):
                atom = traj.topology._atoms[z]
                # print(atom.residue.__dict__.keys())
                # print(atom.residue.chain.index)
                if atom.residue.chain.index > chain_index:
                    chain_index = atom.residue.chain.index
                    chain_id = chr(ord(chain_id) + 1)
                x, y, z = map(self.fixCoordinates, frame[z])
                line = ['ATOM', str(atom.index), atom.name, atom.residue.name, chain_id, str(atom.residue.resSeq), x, y, z, '1.00', '0.00', atom.element.symbol]
                contents.append(line)
            if prefix is None:
                output = os.path.join(self.root, 'clust{}.pdb'.format(str(k).zfill(2)))
            else:
                output = os.path.join(self.root, '{}.clust{}.pdb'.format(prefix, str(k).zfill(2)))
            writePDB(contents, output)
            if prefix is None:
                log.write('Wrote {}\n\n'.format(output))
            else:
                prefix_log.write('Wrote {}\n\n'.format(output))
            if self.verbose:
                print('Wrote {}'.format(output))
            k += 1
        size_df = pd.DataFrame(size_data, columns=['clust', 't', 'nstructs', 'percent'])

        if prefix is None:
            log.write('Clustering job complete. End time {}\n'.format(self.now()))
            size_df.to_csv(os.path.join(self.root, 'size.csv'))
        else:
            prefix_log.write('Clustering job complete. End time {}\n'.format(self.now())) 
            size_df.to_csv(os.path.join(self.root, '{}.size.csv'.format(prefix)))
        if self.verbose:
            print('Clustering job for {} complete\n\n'.format(self.system_name))
        if prefix is not None:
            prefix_log.close()
        log.close()

    def sortByClusterSize(self):
        df = pd.DataFrame(self.X)
        df['labels'] = self.y
        sizes = {}
        sizes_to_sort = []
        for i in np.unique(self.y):
            d = df[df['labels'] == i]
            sizes[len(d)] = i
            sizes_to_sort.append(len(d))
        sizes_sorted = sorted(sizes_to_sort, reverse=True)
        y_sorted = []
        for i in sizes_sorted:
            y_sorted.append(sizes[i])
        return y_sorted

    def getOligomerSize(self, prefix='cluster', cutoff=10):
        df = pd.read_csv(os.path.join(self.root, '{}.size.csv'.format(prefix)), index_col=0, header=0)
        oligomer_types_all = []
        for clust in df['clust']:
            if prefix == 'cluster':
                filename = os.path.join(self.root, 'clust{}.pdb'.format(str(clust).zfill(2)))
            else:
                filename = os.path.join(self.root, '{}.clust{}.pdb'.format(prefix, str(clust).zfill(2)))
            protein = Protein(filename)
            chains = list(protein.chains.keys())
            backbone = protein.backbone
            connectivity = {}
            for i in range(0, len(chains)):
                backbone_a = backbone[chains[i]]
                connectivity[chains[i]] = []
                for k in range(0, len(chains)):
                    if i == k:
                        continue
                    backbone_b = backbone[chains[k]]
                    mindist = protein.mindistAtoms(backbone_a, backbone_b)
                    if mindist <= cutoff:
                        connectivity[chains[i]].append(1)
                    else:
                        connectivity[chains[i]].append(0)
            oligomer_types = []
            for chain, connections in connectivity.items():
                connection_sum = sum(connections) + 1
                oligomer_type = '{}-mer'.format(connection_sum)
                if oligomer_type not in oligomer_types:
                    oligomer_types.append(oligomer_type)
            oligomer_types_all.append('/'.join(oligomer_types))
        df['oligomer_type'] = oligomer_types_all
        return df

                

    def mpCaller(self, data):
        pass

    def fixCoordinates(self, xyz):
        return xyz*10


# run clustering example
# def runClustering():
#     clust = Cluster('cat.pbc.nowat.xtc', 'nowat.top.gro', stride=10)
#     clust.run(nprocs=6, output='cluster.csv')

# if __name__ == '__main__':
#     # run clustering

#     clust = Cluster('cat.pbc.nowat.xtc', 'nowat.top.gro', stride=10)
#     clust.load('cluster.csv')
#     clust.analyze.pca()
#     clust.analyze.kmeans()






