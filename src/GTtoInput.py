import networkx as nx
import numpy as np
import pandas as pd

class GTinput:
    def __init__(self, G):
        self.G = G
        self.k = None
        self.m = None
        self.snvs = [] # accessed through int id - 'u','v', 'segment', 'sample'
        self.cnas = None
        self.computek()
        self.computem()
        self.computesnvs()
        self.computecnas()

    def computek(self):
        for node in self.G.nodes:
            rows = self.G.nodes[node]['label'].split('\\n')
            self.k = len(rows[0].split(';'))-1
            return
    
    def computem(self):
        for node in self.G.nodes:
            rows = self.G.nodes[node]['label'].split('\\n')
            props = rows[-1].split(' ')
            self.m = len(props)-2
            return

    def computesnvs(self):
        numsnvs = 0
        for edge in self.G.edges:
            u,v = edge
            rowsU = self.G.nodes[u]['label'].split('\\n')
            rowsV = self.G.nodes[v]['label'].split('\\n')
            snvsU = []
            snvsV = []
            for segmentIdxxx in range(len(rowsU)):
                if rowsU[segmentIdxxx] == "":
                    snvsU.append([])
                elif rowsU[segmentIdxxx][0] == '0' or rowsU[segmentIdxxx][0] == '1':
                    currSnvs = rowsU[segmentIdxxx].split(';')
                    currSnvs = [ele.strip() for ele in currSnvs]
                    currSnvs = [ele for ele in currSnvs if ele != ""]
                    snvsU.append(currSnvs)
            for segmentIdxxx in range(len(rowsV)):
                if rowsV[segmentIdxxx] == "":
                    snvsV.append([])
                elif rowsV[segmentIdxxx][0] == '0' or rowsV[segmentIdxxx][0] == '1':
                    currSnvs = rowsV[segmentIdxxx].split(';')
                    currSnvs = [ele.strip() for ele in currSnvs]
                    currSnvs = [ele for ele in currSnvs if ele != ""]
                    snvsV.append(currSnvs)
            
            for segmentIdx in range(len(snvsU)):
                for snvIdx in range(len(snvsU[segmentIdx])):

                    if segmentIdx > len(snvsU)-1 or snvIdx > len(snvsU[segmentIdx])-1 or segmentIdx>len(snvsV)-1 or snvIdx>len(snvsV[segmentIdx])-1:
                        continue
                    
                    if snvsU[segmentIdx][snvIdx] != snvsV[segmentIdx][snvIdx]:
                        currU = u
                        currV = v
                        currSegment = segmentIdx-1
                        
                        if snvsU[segmentIdx][snvIdx][0] != snvsV[segmentIdx][snvIdx][0]:
                            cnaIdx = 0
                            cnaIdx2 = 0
                        else:
                            cnaIdx = 1
                            cnaIdx2 = 2
                        
                        # iterate through every node and add its proportions
                        # if there is snv mutation there
                        props = [0 for _ in range(self.m)]
                        for nodee in self.G.nodes:
                            rowss = self.G.nodes[nodee]['label'].split('\\n')
                            cnasofnode = rowss[0].split(';')
                            cnasofnode = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnode]
                            cnasofnode = cnasofnode[:-1]
                            cnasofnode = [tuple(map(int, ele.split(','))) for ele in cnasofnode]
                            propsofnode = rowss[-1].strip('[').strip(']').strip('\"').lstrip().rstrip().split()[:-1]
                            propsofnode = [float(ele) for ele in propsofnode]
                            snvsofnode = []
                            for segmentIdxx in range(len(rowss)):
                                if rowss[segmentIdxx] == "":
                                    snvsofnode.append([])
                                elif rowss[segmentIdxx][0] == '0' or rowss[segmentIdxx][0] == '1':
                                    currrSnvs = rowss[segmentIdxx].split(';')
                                    currrSnvs = [ele.strip() for ele in currrSnvs]
                                    currrSnvs = [ele for ele in currrSnvs if ele != ""]
                                    snvsofnode.append(currrSnvs)

                            if segmentIdx > len(snvsofnode)-1 or snvIdx > len(snvsofnode[segmentIdx])-1 or segmentIdx>len(snvsV)-1 or snvIdx>len(snvsV[segmentIdx])-1:
                                continue
                            if snvsofnode[segmentIdx][snvIdx] == snvsV[segmentIdx][snvIdx]:
                                for mIdx in range(self.m):
                                    props[mIdx] += int(cnasofnode[segmentIdx][cnaIdx])*propsofnode[mIdx]
                        self.snvs.append([numsnvs, u,v,segmentIdx] + props)
                        numsnvs += 1

        self.snvs = np.array(self.snvs)
        self.snvs = pd.DataFrame(self.snvs, columns=['snv_mut', 'u', 'v', 'segment'] + ['sample'+str(midx) for midx in range(self.m)])
        



    def computecnas(self):
        self.cnas = {}
        for segmentidx in range(self.k):
            self.cnas[segmentidx] = set()
            for node in self.G.nodes:
                rows = self.G.nodes[node]['label'].split('\\n')
                cnasofnode = rows[0].split(';')
                cnasofnode = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnode]
                cnasofnode = cnasofnode[:-1]
                cnasofnode = [tuple(map(int, ele.split(','))) for ele in cnasofnode]
                self.cnas[segmentidx].add(cnasofnode[segmentidx])
        for segmentidx,cnaSet in self.cnas.items():
            cnas = list(cnaSet)
            self.cnas[segmentidx] = {}
            for cna in cnas:
                self.cnas[segmentidx][cna] = {}
                for sampleidx in range(self.m):
                    self.cnas[segmentidx][cna]['sample'+str(sampleidx)] = 0
        for segmentidx in range(self.k):
            for node in self.G.nodes:
                rows = self.G.nodes[node]['label'].split('\\n')
                propsofnode = rows[-1].strip('\"').strip('[').strip(']')[1:-1]
                propsofnode = list(map(float, propsofnode.split(' ')))
                cnasofnode = rows[0].split(';')
                cnasofnode = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnode]
                cnasofnode = cnasofnode[:-1]
                cnasofnode = [tuple(map(int, ele.split(','))) for ele in cnasofnode]
                currNode = cnasofnode[segmentidx]
                for sampleidx in range(self.m):
                    self.cnas[segmentidx][currNode]['sample'+str(sampleidx)] += propsofnode[sampleidx]

            

def getGTinput(filename):
    G = nx.Graph(nx.nx_pydot.read_dot(filename))
    return GTinput(G)