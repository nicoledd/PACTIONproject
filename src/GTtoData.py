import networkx as nx
import numpy as np


class Node:
    def __init__(self, name, segments, props):
        self.name = name
        self.segments = segments
        self.props = props
    def __eq__(self, other):
        if len(self.segments) != len(other.segments):
            return False
        return np.all([seg1 == seg2 for (seg1,seg2) in zip(self.segments,other.segments)])


class Segment:
    def __init__(self, cna, snvs):
        self.cna = np.array(cna)
        self.snvs = np.array(snvs)
    def __eq__(self, other):
        return np.array_equal(self.cna, other.cna) and np.array_equal(self.snvs, other.snvs)



class GTdata:
    def __init__(self, G):
        self.G = G
        self.k = 0
        self.numclusters = 0
        self.cnatrees = []
        self.genotypetrees = []
        self.clones = []
        self.cnaclones = []
        self.edges = []
        self.cnaedges = []

        self.computek()
        self.computenumclusters()
        self.computecnatrees()
        self.computegenotypetrees()
        self.computeclones()
        self.computecnaclones()
        #self.computeedges()
        #self.computecnaedges()

    def computek(self):
        for node in self.G.nodes:
            rows = self.G.nodes[node]['label'].split('\\n')
            self.k = len(rows[0].split(';'))-1
            return

    def computegenotypetrees(self):
        return

    '''
    def computecnaedges(self):
        nodes = []
        for node in G.nodes:
            rows = G.nodes[node]['label']
            props = None
            segments = []
            for segmentidx in range(self.k):
                cna = rows[0][segmentidx]
                segments.append(Segment(cna, None))
            nodes.append(Node(node, segments, props))
        T = nx.DiGraph()
        T.add_nodes(nodes)
        for edge in G.edges:
            u,v = edge
            rowsU = self.G.nodes[u]['label'].split('\\n')
            rowsV = self.G.nodes[v]['label'].split('\\n')
            T.add_edge(nodeU, nodeV)
        self.cnaedges = list(T.edges)


    def computeedges(self):
        nodes = []
        for node in G.nodes:
            rows = G.nodes[node]['label']
            props = None
            segments = []
            for segmentidx in range(self.k):
                cna = rows[0][segmentidx]
                snvs = rows[segmentidx+1]
                segments.append(Segment(cna, snvs))
            nodes.append(Node(node, segments, props))
        T = nx.DiGraph()
        T.add_nodes(nodes)
        for edge in G.edges:
            u,v = edge
            rowsU = self.G.nodes[u]['label'].split('\\n')
            rowsV = self.G.nodes[v]['label'].split('\\n')
            T.add_edge(nodeU, nodeV)
        self.edges = list(T.edges)
    '''


    def computeclones(self):
        for node in self.G.nodes:
            rows = self.G.nodes[node]['label'].split('\\n')
            cnasofnode = rows[0].split(';')
            cnasofnode = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnode]
            cnasofnode = cnasofnode[:-1]
            cnasofnode = [tuple(map(int, ele.split(','))) for ele in cnasofnode]
            snvsofnode = []
            for segmentIdx in range(len(rows)):
                if rows[segmentIdx] == "":
                    snvsofnode.append([])
                elif rows[segmentIdx][0] == '0' or rows[segmentIdx][0] == '1':
                    currSnvs = rows[segmentIdx].split(';')
                    currSnvs = [ele.strip() for ele in currSnvs]
                    currSnvs = [ele for ele in currSnvs if ele != ""]
                    snvsofnode.append(currSnvs)
            props = None
            segments = []
            for segmentidx in range(self.k):
                cna = cnasofnode[segmentidx]
                snvs = snvsofnode[segmentidx]
                segments.append(Segment(cna, snvs))
            self.clones.append(Node(node, segments, props))

    def computecnaclones(self):
        for node in self.G.nodes:
            rows = self.G.nodes[node]['label'].split('\\n')
            cnasofnode = rows[0].split(';')
            cnasofnode = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnode]
            cnasofnode = cnasofnode[:-1]
            cnasofnode = [tuple(map(int, ele.split(','))) for ele in cnasofnode]
            props = None
            segments = []
            for segmentidx in range(self.k):
                cna = cnasofnode[segmentidx]
                snvs = None
                segments.append(Segment(cna, snvs))
            self.clones.append(Node(node, segments, props))


    def computenumclusters(self):
        for edge in self.G.edges:
            u,v = edge
            rowsU = self.G.nodes[u]['label'].split('\\n')
            rowsV = self.G.nodes[v]['label'].split('\\n')
            snvsU = []
            snvsV = []
            for segmentIdx in range(len(rowsU)):
                if rowsU[segmentIdx] == "":
                    snvsU.append([])
                elif rowsU[segmentIdx][0] == '0' or rowsU[segmentIdx][0] == '1':
                    currSnvs = rowsU[segmentIdx].split(';')
                    currSnvs = [ele.strip() for ele in currSnvs]
                    currSnvs = [ele for ele in currSnvs if ele != ""]
                    snvsU.append(currSnvs)
            for segmentIdx in range(len(rowsV)):
                if rowsV[segmentIdx] == "":
                    snvsV.append([])
                elif rowsV[segmentIdx][0] == '0' or rowsV[segmentIdx][0] == '1':
                    currSnvs = rowsV[segmentIdx].split(';')
                    currSnvs = [ele.strip() for ele in currSnvs]
                    currSnvs = [ele for ele in currSnvs if ele != ""]
                    snvsV.append(currSnvs)
            if snvsU != snvsV:
                self.numclusters += 1
    
    def computecnatrees(self):
        self.cnatrees = [nx.DiGraph() for _ in range(self.k)]
        for edge in self.G.edges:
            u,v = edge
            rowsU = self.G.nodes[u]['label'].split('\\n')
            cnasofnodeu = rowsU[0].split(';')
            cnasofnodeu = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnodeu]
            cnasofnodeu = cnasofnodeu[:-1]
            cnasofnodeu = [tuple(map(int, ele.split(','))) for ele in cnasofnodeu]
            rowsV = self.G.nodes[v]['label'].split('\\n')
            cnasofnodev = rowsV[0].split(';')
            cnasofnodev = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnodev]
            cnasofnodev = cnasofnodev[:-1]
            cnasofnodev = [tuple(map(int, ele.split(','))) for ele in cnasofnodev]
            for segmentidx in range(self.k):
                if cnasofnodeu[segmentidx] != cnasofnodev[segmentidx]:
                    self.cnatrees[segmentidx].add_edge(cnasofnodeu[segmentidx], cnasofnodev[segmentidx])
                else:
                    if cnasofnodeu[segmentidx] not in self.cnatrees[segmentidx].nodes:
                        self.cnatrees[segmentidx].add_node(cnasofnodeu[segmentidx])
             

def getGTdata(filename):
    G = nx.Graph(nx.nx_pydot.read_dot(filename))
    return GTdata(G)