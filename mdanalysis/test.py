import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
from sklearn import preprocessing

from collections.abc import Iterable

class Network:
    
    def __init__(self, nodes=None, edges=None, node_ids=None, labels=None):
        
        self.nodes = [] if nodes is None else nodes
        self.edges = [] if edges is None else edges
        if nodes is not None:
            self.labels = [[node.label for node in self.nodes if node.bipartite == 0],
                           [node.label for node in self.nodes if node.bipartite == 1]]
#         self.node_ids = [] if node_ids is None else node_ids
        self.labels = None

    @classmethod
    def fromXPM(cls, xpm, num_residues, num_peptides, igncaps=True):
        from xpm import XpmParser
        x = XpmParser(xpm, num_residues, num_peptides, igncaps=True)
        d = x.averageAll()
        t = d.reindex(index=d.index[::-1])
        df = Network.normalize(t)
        df = df.apply(lambda x: 1 - x)
    
        return Network.fromDataFrame(df)

    @classmethod
    def fromMultiXPM(cls, xpm, num_residues_1, num_residues_2, num_peptides_1, num_peptides_2, igncaps=True):
        from mixedxpm import MixedParser
        x = MixedParser(xpm, num_residues_1, num_residues_2, num_peptides_1, num_peptides_2, igncaps=igncaps)
        d = x.averageResiduesMixed()
        t = d.reindex(index=d.index[::-1])
        df = Network.normalize(t)
        df = df.apply(lambda x: 1 - x)
        
        return Network.fromDataFrame(df)
        
    @classmethod
    def fromDataFrame(cls, df):
        nodes = []
        edges = []
        node_id = 0
        for column in df.columns:
            node = Node.fromData(label=column, bipartite=0, id_=node_id)
            nodes.append(node)
            node_id += 1
        bp_end = node_id
        for index in df.index:
            node = Node.fromData(label=index, bipartite=1, id_=node_id)
            nodes.append(node)
            node_id += 1
        b0 = nodes[:bp_end]
        b1 = nodes[bp_end:]
        for i in range(0, len(b0)):
            s = df.iloc[:, i]
            k = 0
            for weight in s:
                edge = Edge(node1=b0[i], node2=b1[k], weight=weight)
                if (weight !=0 ) and (weight != 1):
                    edges.append(edge)
                k += 1
#         node_ids = [i for i in range(0, node_id)]
#         return cls(nodes=nodes, edges=edges, node_ids=node_ids)
        return cls(nodes=nodes, edges=edges)
    
    @property
    def edgelist(self):
        edgelist = []
        for edge in self.edges:
            edgelist.append((edge.nodes[0].id, edge.nodes[1].id))
        return edgelist
        
    @property
    def weightlist(self):
        weightlist = []
        for edge in self.edges:
            weightlist.append(edge.weight)
        return sorted(weightlist)

    @property
    def bins(self):
        #df = pd.DataFrame([edge.weight for edge in self.edges], index=[i for i in range(len(self.edges))], columns=['edges'])
        bins = []
        x = pd.cut(self.weightlist, 10)
        alphas = [0.005, 0.01, 0.015, 0.03, 0.06, 0.08, 0.1, 0.15, 0.3, 0.5]
        i = 0
        for item in x.categories:
            bins.append(Interval(item, alpha=alphas[i]))
            i += 1
        return bins

    @property
    def aaProperties(self):
        properties = {
            'G':'nonpolar',
            'A':'nonpolar',
            'V':'nonpolar',
            'L':'nonpolar',
            'I':'nonpolar',
            'M':'nonpolar',
            'F':'aromatic',
            'W':'aromatic',
            'P':'nonpolar',
            'S':'polar',
            'T':'polar',
            'C':'polar',
            'Y':'polar',
            'N':'polar',
            'Q':'polar',
            'D':'acidic',
            'E':'acidic',
            'K':'basic',
            'R':'basic',
            'H':'basic'
        }
        colors = {
            'nonpolar':'#c7dd92',
            'polar':'#ecc7e0',
            'acidic':'#dc967d',
            'basic':'#afd5f8',
            'aromatic':'#f2edb4'
        }
        props = {}
        if self.labels is None:
            raise ValueError('Must have self.labels, or pass labels here')
        for grp in self.labels:
            for label in grp:
                for char in label:
                    if char.isalpha():
                        prop = properties[char]
                        col = colors[prop]
                        props[label] = col
        return props
        
        
    def addNode(self, label, bipartite=None):
        node_id = len(self.nodes)
        data = {
            'label':label,
            'bipartite':bipartite,
            'id':node_id
        }
        node = Node(data)
#         self.nodes.append(node)
        self.nodes(node)
        
    def addEdge(self, nodes=None, ids=None, weight=None):
        if (ids is None) and (nodes is None):
            raise ValueError('ids OR node objects must be specified')
        if (nodes is not None):
            if (len(nodes) != 2):
                raise ValueError('length of nodes must be 2')
            edge = Edge(nodes[0], nodes[1], weight=weight)
#             self.edges(edge)
            self.edges.append(edge)
            i = 0
            for node in self.nodes:
                if nodes[0] is node:
                    self.nodes[i].edges.append(nodes[1])
                if nodes[1] is node:
                    self.nodes[i].edges.append(nodes[0])
                i += 1
                
        if (ids is not None):
            if (len(ids) != 2):
                raise ValueError('length of ids must be 2')
            edge = Edge(self.nodes[ids[0]], self.nodes[ids[1]], weight=weight)
            self.edges.append(edge)
            self.nodes[ids[0]].edges.append(self.nodes[ids[1]])
            self.nodes[ids[1]].edges.append(self.nodes[ids[0]])
            return edge
    
    def relabelNode(self, label, id_=None, node=None):
        if id_ is not None:
            self.nodes[id_].label = label
            return self.nodes[id_]
        if node is not None:
            i = 0
            for n in self.nodes:
                if node is n:
                    self.nodes[i].label = label
                    return self.nodes[i]
                    
    def relabelNodes(self, l0, l1):
        i = 0
        j = 0
        for node in self.nodes:
            if node.bipartite == 0:
                node.label = l0[i]
                i += 1
            else:
                node.label = l1[j]
                j += 1
        
    def getNetworkXGraph(self):
        b0 = []
        b1 = []
        for node in self.nodes:
            if node.bipartite == 0:
                b0.append(node.id)
            else:
                b1.append(node.id)
        G = nx.Graph()
        G.add_nodes_from(b0, bipartite=0)
        G.add_nodes_from(b1, bipartite=1) 
        return G

    def getPos(self, G=None):
        if G is None:
            G = self.getNetworkXGraph()
        b0 = []
        for node in self.nodes:
            if node.bipartite == 0:
                b0.append(node.id)
        pos = nx.bipartite_layout(G, nodes=b0)
        p = {}
        for key, value in pos.items():
            v = np.array([value[0], value[-1]*-1])
            p[key] = v
            
        return p
    
    def drawNodes(self, ax=None, size=500, color=None, alpha=None, edgecolors=None, linewidths=None):
        G = self.getNetworkXGraph()
        nodelist = list(G.nodes())
        pos = self.getPos(G)
        pos = np.asarray([pos[v] for v in nodelist])
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);
            
        node_collection = ax.scatter(
            pos[:,0],
            pos[:,1],
            c=color,
            s=size,
            alpha=alpha,
            edgecolors=edgecolors,
            linewidths=linewidths
        )
        ax.tick_params(
            axis="both",
            which="both",
            bottom=False,
            left=False,
            labelbottom=False,
            labelleft=False,
        )
        ax.axis('off')
        return node_collection
    
    def drawNodesByProperty(self, ax=None, size=500, alpha=None, edgecolors=None, linewidths=None, labels=None):
        G = self.getNetworkXGraph()
        nodelist = list(G.nodes())
        pos = self.getPos(G)
        pos = np.asarray([pos[v] for v in nodelist])
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);
        if self.labels is None:
            if labels is None:
                raise ValueError('pass labels')
            else:
                if (len(labels) == 2) and (isinstance(labels, list)):
                    l0 = labels[0]
                    l1 = labels[1]
                elif (isinstance(labels, list)):
                    l0 = labels
                    l1 = labels
                else:
                    raise ValueError("for 'list', expected a list of format [[], []] or []")
        self.labels = [l0,l1]
        
        node_collection = []
        i = 0
        
        for l in self.labels:
            for key, value in self.aaProperties.items():
                n = ax.scatter(
                pos[i,0],
                pos[i,1],
                c=value,
                s=size,
                alpha=alpha,
                edgecolors=edgecolors,
                linewidths=linewidths
            )
                i += 1
                node_collection.append(n)
        
        ax.tick_params(
            axis="both",
            which="both",
            bottom=False,
            left=False,
            labelbottom=False,
            labelleft=False,
        )
        ax.axis('off')
        return node_collection
    
    def drawEdges(self, ax=None, width=1.0, alpha=None, edge_color='k', style='solid', 
                  edge_cmap=None, edge_vmin=None,edge_vmax=None,):
        
        from matplotlib.collections import LineCollection
        
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);   
        
        G = self.getNetworkXGraph()
        pos = self.getPos(G)
        edgelist = self.edgelist
        edge_pos = np.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])
        
        
        edge_collection = LineCollection(
            edge_pos,
            colors=edge_color,
            linewidths=width,
            antialiaseds=(1,),
            linestyle=style,
            transOffset=ax.transData,
            alpha=alpha,
        )
        
#         edge_collection.set_cmap(edge_cmap)
#         edge_collection.set_clim(edge_vmin, edge_vmax)
        
        edge_collection.set_zorder(0)
        ax.add_collection(edge_collection)
        return edge_collection

    def drawEdgesByWeight(self, ax=None, width=1.0, edge_color='k', style='solid'):

        from matplotlib.collections import LineCollection
        
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);   
        
        G = self.getNetworkXGraph()
        pos = self.getPos(G)
        edge_collection = []
#         low_bins = []
#         med_bins = []
#         top_bins = []
        for edge in self.edges:
            for b in self.bins:
                if (b.left < edge.weight) and (b.right > edge.weight):
                    alpha = b.alpha
                    if alpha > 1:
                        alpha = 1
            
            edge_pos = np.asarray([(pos[edge.nodes[0].id], pos[edge.nodes[1].id])])
#             alpha = edge.weight
            if alpha > 1:
                alpha = 1
            e = LineCollection(
                edge_pos,
                colors=edge_color,
                linewidths=width,
                antialiaseds=(1,),
                linestyle=style,
                transOffset=ax.transData,
                alpha=alpha,
            )

            e.set_zorder(0)
            ax.add_collection(e)
            edge_collection.append(e)
#         print(len(low_bins), len(med_bins), len(top_bins))
        return edge_collection


    def drawLabels(self, ax=None, labels=None, font_size=12, font_family='sans-serif'):
        
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);   
        
        G = self.getNetworkXGraph()
        pos = self.getPos(G)
        
        if labels is not None:
            if (len(labels) == 2) and (isinstance(labels, list)):
                self.relabelNodes(labels[0], labels[1])
            elif (isinstance(labels, list)):
                self.relabelNodes(labels, labels)
            else:
                raise ValueError("for 'list', expected a list of format [[], []] or []")
        
        text_items = {}
        for node in self.nodes:
            (x, y) = pos[node.id]
            label = node.label
            if not isinstance(label, str):
                label = str(label)  # this makes "1" and 1 labeled the same
            
            t = ax.text(x, y, label, size=font_size, family=font_family, horizontalalignment='center', 
                        verticalalignment='center', transform=ax.transData, clip_on=True)
            text_items[node.id] = t

        
        return text_items
    
    @staticmethod
    def normalize(df):
        
        def norm(x, min_, max_):
            return (x - min_) / (max_ - min_)
        
        mi = None
        ma = None
        for column in df.columns:
            for index in df.index:
                val = df.loc[index, column]
                if mi is None:
                    mi = val
                if ma is None:
                    ma = val
                if val > ma:
                    ma = val
                if val < mi:
                    mi = val
        print(mi,ma)
        df = df.apply(lambda x: (x-mi)/(ma-mi))
        return df
        
class Node:
    
    def __init__(self, data):
        self.label = data['label']
        self.bipartite = data['bipartite']
        self.id = data['id']
        self.edges = []
    
    @classmethod
    def fromData(cls, label, bipartite, id_):
        data = {
            'label':label,
            'bipartite':bipartite,
            'id':id_
        }
        return cls(data)
class Edge:
    
    def __init__(self, node1, node2, weight):
        self.nodes = (node1, node2)
        self.weight = weight

class Interval:
    
    def __init__(self, interval, alpha=None):
        self.left = interval.left
        self.right = interval.right
#         self.alpha = (self.left + self.right) / 2
        if alpha is None:
            self.alpha = self.right
        else:
            self.alpha = alpha

n = Network.fromMultiXPM(['D:/Work/grant/MDmat_mean/Mix1_MDmat_mean.xpm', 'D:/Work/grant/MDmat_mean/Mix2_MDmat_mean.xpm', 
                     'D:/Work/grant/MDmat_mean/Mix3_MDmat_mean.xpm', 'D:/Work/grant/MDmat_mean/Mix4_MDmat_mean.xpm'],
                    num_residues_1=7, num_peptides_1=3, num_residues_2=10, num_peptides_2=3, igncaps=True)
fig = plt.figure(figsize=(16,12));
ax = fig.add_subplot(121);
mix_labels = [ab_labels, iapp_labels]
n.drawNodesByProperty(ax=ax, labels=mix_labels, size=900, edgecolors='k')
n.drawEdgesByWeight(ax=ax)
n.drawLabels(ax=ax, labels=mix_labels)
plt.show()