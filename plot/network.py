import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import networkx as nx
import numpy as np
import pandas as pd
from pymd.utilities.library import code_conversions
from pymd.utilities.custom_colormaps import create_colormap
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from pymd.utilities.library import code_conversions
import mdtraj 

class Network:
    
    def __init__(self, nodes=None, edges=None, pdb=None, node_ids=None, labels=None):
        
        self.nodes = [] if nodes is None else nodes
        self.edges = [] if edges is None else edges
        self.pdb = pdb
        self.labels = labels
        if self.pdb is not None:
            self.top = mdtraj.load(self.pdb).topology
        else:
            self.top = None
        if labels is None:
            if self.nodes == []:
                self.labels = [[node.label for node in self.nodes if node.bipartite == 0],
                            [node.label for node in self.nodes if node.bipartite == 1]]
            else:
                self.labels = []
        self.properties = {
            'G':'nonpolar',
            'A':'nonpolar',
            'V':'nonpolar',
            'L':'nonpolar',
            'I':'nonpolar',
            'M':'nonpolar',
            'F':'aromatic',
            'W':'aromatic',
            'P':'nonpolar',
            'Y':'aromatic',
            'S':'polar',
            'T':'polar',
            'C':'polar',
            'N':'polar',
            'Q':'polar',
            'D':'acidic',
            'E':'acidic',
            'K':'basic',
            'R':'basic',
            'H':'polar'
        }
#         self.node_ids = [] if node_ids is None else node_ids
        self.G = nx.Graph()
        self.getNetworkXGraph()
        self.axis = 0

    @classmethod
    def fromDataFrame(cls, df, pdb=None, scale=True):
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
        if scale:
            scaler = MinMaxScaler((0,1))
            scaler.fit(df)
            X = scaler.transform(df)
            d = pd.DataFrame(X)
            d.index = df.index
            d.columns = df.columns
            df = d
        for i in range(0, len(b0)):
            s = df.iloc[:, i]
            k = 0
            for weight in s:
                edge = Edge(node1=b0[i], node2=b1[k], weight=weight)
                if (weight !=0 ):
                    edges.append(edge)
                k += 1
        return cls(nodes=nodes, edges=edges, pdb=pdb)

    @property
    def edgelist(self):
        edgelist = []
        for edge in self.edges:
            edgelist.append((edge.nodes[0], edge.nodes[1]))
        return edgelist
    
    def aaProperties(self, labels):
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
                    'Y':'aromatic',
                    'S':'polar',
                    'T':'polar',
                    'C':'polar',
                    'N':'polar',
                    'Q':'polar',
                    'D':'acidic',
                    'E':'acidic',
                    'K':'basic',
                    'R':'basic',
                    'H':'polar'
                }
        colors = {
            'nonpolar':'#c7dd92',
            'polar':'#ecc7e0',
            'acidic':'#dc967d',
            'basic':'#afd5f8',
            'aromatic':'#f2edb4'
        }
        props = {}
        for grp in labels:
            for label in grp:
                for char in label:
                    if char.isalpha():
                        prop = properties[char]
                        col = colors[prop]
                        props[label] = col
        return props
    
    def unique(self, prop='bipartite'):
        props = []
        if hasattr(self.nodes[0], prop):
            for node in self.nodes:
                if node.__dict__[prop] not in props:
                    props.append(node.__dict__[prop])
        return props
    
    def assignProps(self, labels, props):
        self.properties = {
            'G':'nonpolar',
            'A':'nonpolar',
            'V':'nonpolar',
            'L':'nonpolar',
            'I':'nonpolar',
            'M':'nonpolar',
            'F':'aromatic',
            'W':'aromatic',
            'P':'nonpolar',
            'Y':'aromatic',
            'S':'polar',
            'T':'polar',
            'C':'polar',
            'N':'polar',
            'Q':'polar',
            'D':'acidic',
            'E':'acidic',
            'K':'basic',
            'R':'basic',
            'H':'polar'
        }
        if (len(labels) == 2) and (isinstance(labels, list)):
            labs = labels[0] + labels[1]
        elif (isinstance(labels, list)):
            if len(labels) == len(self.nodes):
                labs = labels
            else:
                labs = labels + labels
        else:
            raise ValueError('wrong data type for labels :(((')
        nodes = []
        for node, l in zip(self.nodes, labs):
            code = l[0]
            color = props[code]
            node.color = color
            node.code = code
            node.label = l
            node.prop = self.properties[l[0]]
            nodes.append(node)
        self.nodes = nodes
        self.updateEdges()
        return self.nodes
    
    def updateEdges(self):
        edges = []
        for edge in self.edges:
            n1 = edge.nodes[0].id
            n2 = edge.nodes[1].id
            nodes = (self.nodes[n1], self.nodes[n2])
            edge.nodes = nodes
            edges.append(edge)
        self.edges = edges
        self.getEdgePos()
        return self.edges

    def getNetworkXGraph(self):
        b0 = []
        b1 = []
        for node in self.nodes:
            if node.bipartite == 0:
                b0.append(node.id)
            else:
                b1.append(node.id)
        self.G = nx.Graph()
        self.G.add_nodes_from(b0, bipartite=0)
        self.G.add_nodes_from(b1, bipartite=1) 
        self.getNodePos()
        self.getEdgePos()
        return self.G

    def getNodePos(self, G=None):
        if G is None:
            G = self.G
        b0 = []
        for node in self.nodes:
            if node.bipartite == 0:
                b0.append(node.id)
        pos = nx.bipartite_layout(G, nodes=b0)
        p = []
        for key, value in pos.items():
            v = np.array([value[0], value[-1]*-1])
            v = [v[0], v[1]]
            p.append(v)
        pos = np.array(p)
        nodes = []
        for (node, xy) in zip(self.nodes, pos):
            node.pos = np.array((node.bipartite, xy[1]))
            nodes.append(node)
        self.nodes = nodes
        self.updateEdges()
        return pos
    
    @property
    def npos(self):
        pos = [node.pos for node in self.nodes]
        return np.array(pos)
    def getEdgePos(self, G=None):
        if G is None:
            G = self.G
        pos = []
        edges = []
        for edge in self.edges:
            xy = np.array((edge.nodes[0].pos, edge.nodes[1].pos))
            pos.append(xy)
            edge.pos = xy
            edges.append(edge)
        self.edges = edges
        return pos
            
    def scaleNodePos(self, axis=1, s=0.2):
        nodes = []
        b = 0
        b0 = [node for node in self.nodes if node.bipartite == 0]
        b1 = [node for node in self.nodes if node.bipartite == 1]
        bs = [b0, b1]
        last = 0
        for b in bs:
            for i, node in enumerate(b):
                if i == 0:
                    nodes.append(node)
                    last = node.pos[axis]
                else:
                    node.pos[axis] = last - s
                    nodes.append(node)
                    last = node.pos[axis]
        self.nodes = nodes
        self.updateEdges()
        return self.nodes

    def transpose(self):
        nodes = []
        for node in self.nodes:
            node.pos = np.array((node.pos[1], node.pos[0]))
            nodes.append(node)
        self.node = nodes
        self.updateEdges()
        if self.axis == 0:
            self.axis = 1
        else:
            self.axis = 0
        return self
    
    def mirror(self):
        nodes = []
        bs = []
        for b in np.unique(self.npos[:,self.axis]):
            pos = self.npos[self.npos[:, self.axis] == b][::-1]
            i = 0
            for node in self.nodes:
                if node.bipartite == b:
                    try:
                        node.pos = pos[i]
                    except:
                        print(i)
                    nodes.append(node)
                    i += 1
        self.nodes = nodes
        self.updateEdges()
        return self



    
    def drawNodes(self, ax=None, size=500, color='tab:blue', alpha=None, edgecolors='black', linewidths=None):
        # self.getNetworkXGraph()
        # pos = self.getNodePos()
        # pos = np.asarray([pos[v] for v in nodelist])
        # print(pos)
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);
            ax.set_facecolor('white')
        pos = np.array([node.pos for node in self.nodes])
        node_collection = ax.scatter(
            pos[:,0],
            pos[:,1],
            c=color,
            s=size,
            alpha=alpha,
            edgecolors=edgecolors,
            linewidths=linewidths
        )

        return ax
    
    def drawNodesByProperty(self, labels, ax=None, size=500, edgecolors='black', linewidths=None, transp=[]):
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);
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
        props = self.aaProperties(labels)
        self.assignProps(labels, props)
        node_collection = []

        # pos = self.getNodePos()
        for (i, node) in enumerate(self.nodes):
            if transp == []:
                alpha = 1
            else:
                if node.prop in transp:
                    alpha = 0.3
                else:
                    alpha = 1
            n = ax.scatter(
                node.pos[0],
                node.pos[1], 
                c=node.color,
                s=size,
                edgecolor=edgecolors,
                linewidth=linewidths,
                alpha=alpha
            )
            node_collection.append(n)
        return ax

    def drawEdges(self, ax=None, width=1.0, alpha='weight', color='dimgray', style='solid', 
                  edge_cmap=None, edge_vmin=None,edge_vmax=None,):
        
        from matplotlib.collections import LineCollection
        
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);   
        
        # G = self.getNetworkXGraph()
        # edgelist = self.edgelist
        # edge_pos = np.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])
    
        for edge in self.edges:
            pos = edge.pos
            if alpha=='weight':
                a = edge.weight
                if a > 1:
                    a = 1
            else:
                a = alpha
            ax.plot(pos[:,0], pos[:,1], color=color, alpha=a, zorder=0, linewidth=width, linestyle=style)
        return ax
    
    def drawLabels(self, labels, ax=None, transp=[], **kwargs):
        if ax is None:
            fig = plt.figure(figsize=(16,12));
            ax = fig.add_subplot(121);
            ax.set_facecolor('white')
        props = self.aaProperties(labels)
        self.assignProps(labels, props)
        for node in self.nodes:
            if transp == []:
                alpha = 1
            else:
                if node.prop in transp:
                    alpha = 0.3
                else:
                    alpha = 1
            ax.text(node.pos[0], node.pos[1], node.label, ha='center', va='center', alpha=alpha, **kwargs)
        return ax

    def setLabels(self, pdb=None, chain_index='all', ignore=[]):
        if (pdb is None) and (self.top is None):
            raise ValueError('No .pdb loaded. Pass a .pdb file in this function or when instantiating Network object.')
        if (self.top is None) and (pdb is not None):
            self.pdb = pdb
            self.top = mdtraj.load(self.pdb).topology
        conversions = code_conversions()
        labels = []
        for residue in self.top.residues:
            if residue.name in ignore:
                continue
            if (chain_index != 'all') and (residue.chain.index != chain_index):
                continue
            l = conversions[residue.name.lower()].upper() + '$_{' + str(residue.resSeq) + '}$'
            labels.append(l)
        self.labels = labels
        return labels

    def draw(self, ax=None, size=400, edge_color='dimgray', linewidth=1, style='solid', transp=[], labels=None, font_dict={}):
        if labels is None:
            if self.labels is not None:
                labels = self.labels
            if self.pdb is not None:
                labels = self.setLabels(self.pdb)
        self = self.mirror()
        ax = self.drawNodesByProperty(labels, ax=ax, size=size, edgecolors='black', transp=transp)
        ax = self.drawEdges(ax=ax, width=linewidth, style=style, color=edge_color)
        ax = self.drawLabels(labels, ax=ax, transp=transp, **font_dict)
        ax.tick_params(
            axis="both",
            which="both",
            bottom=False,
            left=False,
            labelbottom=False,
            labelleft=False,
        )
        ax.set_clip_on(False)
        ax.axis('off')
        return ax


class Node:
    
    def __init__(self, data):
        self.label = data['label']
        self.bipartite = data['bipartite']
        self.id = data['id']
        self.edges = []
        self.pos = np.array([0,0])
        self.code = 'X'
        self.color = 'tab:blue'    
        self.label = 'X'
        self.alpha = 1
        self.prop = None
    @classmethod
    def fromData(cls, label, bipartite, id_):
        data = {
            'label':label,
            'bipartite':bipartite,
            'id':id_
        }
        return cls(data)
    
    @property
    def x(self):
        return self.pos[0]
    
    @property
    def y(self):
        return self.pos[1]

class Edge:
    
    def __init__(self, node1, node2, weight):
        self.nodes = (node1, node2)
        self.weight = weight
        self.pos = np.array([])
class Interval:
    
    def __init__(self, interval, alpha=None):
        self.left = interval[0]
        self.right = interval[1]
#         self.alpha = (self.left + self.right) / 2
        if alpha is None:
            self.alpha = self.right
        else:
            self.alpha = alpha