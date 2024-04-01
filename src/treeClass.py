import networkx as nx
import pygraphviz as pgv

from networkx.drawing.nx_agraph import to_agraph


import numpy as np
import pandas as pd




class Tree:
    def __init__(self, T):

        # self.T.nodes[node]['name'] : gives string name of node
        # self.T.nodes[node]['props'] : gives a float array of props for node
        # self.T.nodes[node]['segments'][segmentidx] : gives a 'cna' dict and 'snvs' list
            # self.T.nodes[node]['segments'][segmentidx]['cna'] : the cna state
            # self.T.nodes[node]['segments'][segmentidx]['snvs'] : list of snv states
        self.T = T
        self.k = self.computek()
        self.m = self.computem()
        self.dot = self.dotFile()
        
        # accessed through - 'u', 'v', 'segment', 'snvidx', 'cnaidx', 'snv_val', and samples
        self.snvs = self.computeVAFs()
        self.cnas = self.computeCNAs()

        # self.segments[kidx] : gives a 'cnas' dict and 'snvs' dict and 'n' num of snvs int
            # self.segments[kidx]['cnas'] : information about cnas
                # self.segments[kidx]['cnas']['nx'] : gives cna graph of this segment
                    # self.segments[kidx]['cnas']['nx'][node]['props'] : gives props of the cnas
                # self.segments[kidx]['cnas']['dot'] : gives cna graph of this segment, dot form
            # self.segments[kidx]['snvs'][snvidx] : information about snvs
                # self.segments[kidx]['snvs'][snvidx]['nx'] : gives genotype graph of this snv
                # self.segments[kidx]['snvs'][snvidx]['dot'] : gives genotype graph of this snv, dot form
                # self.segments[kidx]['snvs'][snvidx]['vaf'] : gives vafs of this snv
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
                    for midx in range(self.m):
                        self.segments[kidx]['cnas']['nx'].nodes[state]['sample' + str(midx)] = 0
            for kidx in range(self.k):
                for node in self.T.nodes:
                    props = self.T.nodes[node]['props']
                    state = str(self.T.nodes[node]['segments'][kidx]['cna'])
                    for midx in range(self.m):
                        self.segments[kidx]['cnas']['nx'].nodes[state]['sample'+str(midx)] += props[midx]
        computeTrees()
        dotFile()
        computeProps()



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
                            self.segments[kidx]['snvs'][snvidx]['nx'].add_edge(str([u1,u2,u3,u4]), str([v1,v2,v3,v4]))
        def computeVafs(): # TODO compute vafs of tree
            pass
        computeTrees()
        dotFile()
        computeVafs()


   

    def dotFile(self): # TODO: create dot file of tree, image display

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


    def computeVAFs(self):
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



    def computeCNAs(self):
        self.cnas = {}
        for kidx in range(self.k):
            self.cnas[kidx] = set()
            for node in self.T.nodes:
                currCNA = self.T.nodes[node]['segments'][kidx]['cna']
                self.cnas[kidx].add(currCNA)
        for kidx,cnaSet in self.cnas.items():
            cnas = list(cnaSet)
            self.cnas[kidx] = {}
            for cna in cnas:
                self.cnas[kidx][cna] = {}
                for sampleidx in range(self.m):
                    self.cnas[kidx][cna]['sample' + str(sampleidx)] = 0
        for kidx in range(self.k):
            for node in self.T.nodes:
                props = self.T.nodes[node]['props']
                currCNA = self.T.nodes[node]['segments'][kidx]['cna']
                for sampleidx in range(self.m):
                    self.cnas[kidx][currCNA]['sample'+str(sampleidx)] += props[sampleidx]
        return self.cnas