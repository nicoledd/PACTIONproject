import networkx as nx
import numpy as np
import pandas as pd
from fractions import Fraction



'''
class: Tree
attributes:
    T: dict. Holds information for each tree node
    k: int. Number of segments
    m: int. Number of samples
    snvs: pandas dataframe.
    cnas: dict.
    segments: dict. Holds information for each tree segment
methods:
    computeSegments
    computeCNAs
    computeVAFs
    obtainSNVprops
    computem
'''
class Tree:
    def __init__(self, T):

        '''
        self.T.nodes[node]['name'] : String. name of node
        self.T.nodes[node]['props'] : Float list. Proportions for node
        self.T.nodes[node]['segments'][segmentidx] : Dict. Information for each segment
            self.T.nodes[node]['segments'][segmentidx]['cna']: Tuple? the cna state
            self.T.nodes[node]['segments'][segmentidx]['snvs'] : Pandas dataframe? list of snv states
        '''
        self.T = T
        self.k = self.computek()
        self.m = self.computem()
        self.dot = self.dotFile()
        
        # accessed through - 'u', 'v', 'segment', 'snvidx', 'cnaidx', 'snv_val', and samples
        self.snvs = self.computeSNVold()


        '''
        self.segments[kidx] : dict. info on each segment
            self.segments[kidx]['cnas'] : dict. info on CNAs of segment
                self.segments[kidx]['cnas']['df'] : pandas dataframe. CNAs and props
                self.segments[kidx]['cnas']['dot'] : dot string. CNA graph of segment
                self.segments[kidx]['cnas']['nx'] : networkx graph. CNA graph of segment
                    self.segments[kidx]['cnas']['nx'].nodes[node]['props'] : list? gives props of the cnas
            self.segments[kidx]['snvs'] : dict. info on SNVs of segment
                self.segments[kidx]['snvs'][snvidx] : dict. info on SNVs of segment
                    self.segments[kidx]['snvs'][snvidx]['nx'] : networkx graph. genotype graph of snv
                    self.segments[kidx]['snvs'][snvidx]['dot'] : dot string. genotype graph of snv
                    self.segments[kidx]['snvs'][snvidx]['vaf'] : list. gives vafs of snv
                    self.segments[kidx]['snvs'][snvidx]['node'] : the node the snv originated in
                    self.segments[kidx]['snvs'][snvidx]['ref'] : num reference alleles (tot - mut)
                    self.segments[kidx]['snvs'][snvidx]['mut'] : num mutated alleles
            self.segments[kidx]['n']: int. number of snvs in segment
        '''
        self.segments = {}
        for kidx in range(self.k):
            self.segments[kidx] = {}
        self.computeSegments()


    def computeSegments(self):
        self.computeCNAsegments()
        self.computeSNVsegments()


    def computeCNAsegments(self):
        def dotFile():
            for kidx in range(self.k):
                A = "digraph T {\n"
                for edge in self.segments[kidx]['cnas']['nx'].edges:
                    u,v = edge
                    A += "\"" + str(u) + "\" -> \"" + str(v) + "\"\n"
                A += "}"
                self.segments[kidx]['cnas']['dot'] = A
        def computeTrees():
            for kidx in range(self.k):
                self.segments[kidx]['cnas'] = {}
                self.segments[kidx]['cnas']['nx'] = nx.DiGraph()
                for edge in self.T.edges:
                    u,v = edge
                    uState = self.T.nodes[u]['segments'][kidx]['cna']
                    vState = self.T.nodes[v]['segments'][kidx]['cna']
                    if uState == vState:
                        self.segments[kidx]['cnas']['nx'].add_node(str(uState))
                    else:
                        self.segments[kidx]['cnas']['nx'].add_edge(str(uState), str(vState))
        def computeProps():
            for kidx in range(self.k):
                for state in self.segments[kidx]['cnas']['nx'].nodes:
                    self.segments[kidx]['cnas']['nx'].nodes[state]['props'] = [0 for _ in range(self.m)]
            for kidx in range(self.k):
                for node in self.T.nodes:
                    props = self.T.nodes[node]['props']
                    state = str(self.T.nodes[node]['segments'][kidx]['cna'])
                    for midx in range(self.m):
                        self.segments[kidx]['cnas']['nx'].nodes[state]['props'][midx] += props[midx]
        def computeDfs():
            for kidx in range(self.k):
                rows = []
                for cna in self.segments[kidx]['cnas']['nx'].nodes:
                    row = [cna] + self.segments[kidx]['cnas']['nx'].nodes[cna]['props']
                    rows.append(row)
                self.segments[kidx]['cnas']['df'] = pd.DataFrame(rows, columns =['copy_number_state'] + ['sample_' + str(i) for i in range(self.m)])

        computeTrees()
        dotFile()
        computeProps()
        computeDfs()



    def computeSNVsegments(self):
        def dotFile():
            for kidx in range(self.k):
                A = "digraph T {\n"
                for snvidx in range(len(self.segments[kidx]['snvs'])):
                    for edge in self.segments[kidx]['snvs'][snvidx]['nx'].edges:
                        u,v = edge
                        A += "\"" + str(u) + "\" -> \"" + str(v) + "\"\n"
                    A += "}"
                    self.segments[kidx]['snvs'][snvidx]['dot'] = A
        def computeTrees():
            for kidx in range(self.k):
                for node in self.T.nodes:
                    snvs = self.T.nodes[node]['segments'][kidx]['snvs']
                    self.segments[kidx]['n'] = len(snvs)
                    self.segments[kidx]['snvs'] = {}
                    for snvidx in range(len(snvs)):
                        self.segments[kidx]['snvs'][snvidx] = {}
                        self.segments[kidx]['snvs'][snvidx]['nx'] = nx.DiGraph()
                        continue
            for edge in self.T.edges:
                u,v = edge
                for kidx in range(self.k):
                    u1,u2 = self.T.nodes[u]['segments'][kidx]['cna']
                    v1,v2 = self.T.nodes[v]['segments'][kidx]['cna']
                    uSnvs = self.T.nodes[u]['segments'][kidx]['snvs']
                    vSnvs = self.T.nodes[v]['segments'][kidx]['snvs']
                    for snvidx in range(len(uSnvs)):
                        u3,u4 = self.T.nodes[u]['segments'][kidx]['snvs'][snvidx]
                        v3,v4 = self.T.nodes[v]['segments'][kidx]['snvs'][snvidx]
                        if [u1, u2, u3, u4] == [v1, v2, v3, v4]:
                            self.segments[kidx]['snvs'][snvidx]['nx'].add_node(str([u1,u2,u3,u4]))
                        else:
                            self.segments[kidx]['snvs'][snvidx]['node'] = v
                            self.segments[kidx]['snvs'][snvidx]['nx'].add_edge(str([u1,u2,u3,u4]), str([v1,v2,v3,v4]))
        def computeVafs():
            for kidx in range(self.k):
                for snvidx in range(self.segments[kidx]['n']):
                    self.segments[kidx]['snvs'][snvidx]['vaf'] = [0 for _ in range(self.m)]
                    for node in self.T.nodes:
                        for matpat in range(2):
                            if self.T.nodes[node]['segments'][kidx]['snvs'][snvidx][matpat] > 0:
                                for midx in range(self.m):
                                    self.segments[kidx]['snvs'][snvidx]['vaf'][midx] += self.T.nodes[node]['props'][midx]*self.T.nodes[node]['segments'][kidx]['snvs'][snvidx][matpat]
        def compute_refs():
            for kidx in range(self.k):
                for snvidx in range(self.segments[kidx]['n']):
                    self.segments[kidx]['snvs'][snvidx]['ref'] = []
                    self.segments[kidx]['snvs'][snvidx]['mut'] = []
                    for midx in range(self.m):
                        f = Fraction(self.segments[kidx]['snvs'][snvidx]['vaf'][midx]).limit_denominator(1000)
                        self.segments[kidx]['snvs'][snvidx]['ref'].append(f.denominator - f.numerator)
                        self.segments[kidx]['snvs'][snvidx]['mut'].append(f.numerator)

        computeTrees()
        dotFile()
        computeVafs()
        compute_refs()


   

    def dotFile(self):
        A = "digraph T {\n"
        for node in self.T.nodes:
            A += self.T.nodes[node]['name'] + "[label=\"name: " + str(node) + "\\n"
            A += "props: " + str(self.T.nodes[node]['props']) + "\\n"
            for kidx in range(self.k):
                A += str(self.T.nodes[node]['segments'][kidx]['cna']) + " , " + str(self.T.nodes[node]['segments'][kidx]['snvs'])
                A += "\\n"
            A += "\"]\n"
        for edge in self.T.edges:
            u,v = edge
            A += str(u) + " -> " + str(v) + "\n"
        A += "}"
        return A


    def computek(self):
        for node in self.T.nodes:
            return len(list(self.T.nodes[node]['segments']))
    
    def computem(self):
        for node in self.T.nodes:
            return len(list(self.T.nodes[node]['props']))





    def obtainSNVprops(self, kidx, snvidx, cnaidx, val):
        def obtainVafDenom():
            vafdenom = [0 for _ in range(self.m)]
            for node in self.T.nodes:
                propsnode = self.T.nodes[node]['props']
                cnasnode = self.T.nodes[node]['segments'][kidx]['cna']
                for midx in range(self.m):
                    vafdenom[midx] += int(cnasnode[cnaidx]) * propsnode[midx]
            return vafdenom
        propsret = [0 for _ in range(self.m)]
        for node in self.T.nodes:
            propsnode = self.T.nodes[node]['props']
            snvsnode = self.T.nodes[node]['segments'][kidx]['snvs']
            cnasnode = self.T.nodes[node]['segments'][kidx]['cna']
            if snvidx<len(snvsnode) and snvsnode[snvidx][cnaidx]==val:
                for midx in range(self.m):
                    propsret[midx] += int(cnasnode[cnaidx]) * propsnode[midx]
        vafdenom = obtainVafDenom()
        for i in range(self.m):
            propsret[i] /= vafdenom[i]
        return propsret


    def computeSNVold(self):
        SNVs = []
        numsnvs = 0

        for edge in self.T.edges:
            u,v = edge
            uSNVs = []
            vSNVs = []
            for kidx in range(self.k):
                uSNVs.append(self.T.nodes[u]['segments'][kidx]['snvs'])
                vSNVs.append(self.T.nodes[v]['segments'][kidx]['snvs'])
            
            for kidx in range(self.k):
                for snvidx in range(len(uSNVs[kidx])):
                    for cnaidx in [0,1]:
                        if uSNVs[kidx][snvidx][cnaidx] < vSNVs[kidx][snvidx][cnaidx]:
                            val = vSNVs[kidx][snvidx][cnaidx]
                            props = self.obtainSNVprops(kidx, snvidx, cnaidx, val)
                            SNVs.append([numsnvs, u, v, int(kidx)] + props)
                            numsnvs += 1

        SNVs = np.array(SNVs)
        SNVs = pd.DataFrame(SNVs, columns=['snv_mut', 'u', 'v', 'segment'] + ['sample'+str(midx) for midx in range(self.m)])
        return SNVs

