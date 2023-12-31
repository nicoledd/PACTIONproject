import csv
from csv import writer
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib.patches import Rectangle
from networkx.drawing.nx_pydot import graphviz_layout
import clonelib


from linear_programs.paction_solver import PCTIsolver
from preprocess.clean_data import getIndexOfTrueGraphs, processSolutionTree, processTreeFile, ifTrueTreesAreFoundIn, getTreeEdges, readProportionMatrices, getAllPossibleTrees
from postprocess.write import writeText, writeTree, writeClones
from modality.modality_class import Modality
from progressive_caller import solveProgressivePCI, solveProgressivePaction


def runAlgo(GTinput):

    s = snvToIdxDict(GTinput.snvs) # snv name -> idx

    trees = []

    for segmentidx in range(GTinput.k):
        currSnvDf = snvDf.loc[snvDf['segment'] == segmentidx]
        currCnaDf = pd.read_csv(args.c[segmentidx], sep='\t')
        T = segmentTree(currSnvDf, currCnaDf, s, GTinput.m)
        '''
        plt.figure(figsize=(15,15))
        pos = graphviz_layout(T, prog="dot")
        nx.draw(T, with_labels = True, font_size=10, pos=pos)
        plt.savefig('segment_' + str(segmentidx) + '.png')
        '''
        trees.append(T)

    finalT = solveProgressivePaction(trees, GTinput.m)
    for node in finalT.nodes:
        firstSegment = None
        nodeName = ""
        for segment in node:
            if len(segment) == 2:
                if firstSegment == None:
                    firstSegment = (segment,)
                else:
                    firstSegment = firstSegment + (segment,)
            else:
                nodeName += str(segment) + '\n'
        finalT.nodes[node]['realName'] = str(firstSegment) + '\n' + str(nodeName)
            
    #plt.figure(figsize=(100,100))
    #labels = nx.get_node_attributes(finalT, 'realName') 
    #pos = nx.nx_agraph.graphviz_layout(finalT, prog="dot")
    #nx.draw(finalT, labels=labels, font_size=10, pos=pos)
    #plt.savefig('finalTree.png')

    predictedNodes = []
    for nodeName in finalT.nodes:
        nodeName0 = list(nodeName)
        firstSection = [ele for ele in nodeName0 if type(ele[0])==int]
        nodeName1 = [tuple(firstSection)] + [ele for ele in nodeName0 if type(ele[0])!=int]
        segments = []
        for segment in nodeName1:
            copies = [segment[0][0], segment[0][1]]
            mutations = [0 for _ in range(len(segment)-1)]
            for segIdx in range(len(segment)-1):
                mutIdx0 = snvDf.at[segment[segIdx+1][0],'snv_mut']
                mutIdx1 = tuple(mutIdx0.split(','))
                mutIdx2 = int(mutIdx1[2])

                mutations[mutIdx2] = segment[segIdx+1][1]
            s = Segment(copies, mutations)
            segments.append(s)
        props = []
        for sampleIdx in range(GTinput.m):
            props.append(finalT.nodes[nodeName]['sample_' + str(sampleIdx)])
        newNode = Node(nodeName, segments, props)
        predictedNodes.append(newNode)

    realNodes = getRealNodes(args.g)  
    recall = getRecall(realNodes, predictedNodes)
    print(recall)

    return finalT




# processing and comparison funcs

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
    def __init__(self, copies, mutated):
        self.copies = np.array(copies)
        self.mutated = np.array(mutated)
    def __eq__(self, other):
        return np.array_equal(self.copies, other.copies) and np.array_equal(self.mutated, other.mutated)


def getRealNodes(filename):
    G = nx.Graph(nx.nx_pydot.read_dot(filename))
    nodes = []
    for node in list(G.nodes):
        label0 = G.nodes[node]['label']
        label1 = label0.split('\\n')
        label2 = [ele.replace("\"",'').replace('(', '').replace(')','').strip() for ele in label1]
        label25 = [ele.replace(" ","") if ele[-1]!=']' else ele for ele in label2]
        label3 = [ele for ele in label25 if ele[-1]==';' or ele[-1]==']']
        copies = [ele for ele in label3[0].split(';') if ele!='']
        segments = []
        for i in range(len(copies)):
            currCopy0 = copies[i]
            currCopy1 = currCopy0.split(',')
            currCopy3 = [int(currCopy1[0]), int(currCopy1[1])]
            mutated0 = label3[i+1]
            mutated1 = mutated0.split(';')
            mutated2 = [1 if ele.count('1')>0 else 0 for ele in mutated1]
            s = Segment(currCopy3, mutated2)
            segments.append(s)
        
        props0 = label3[-1]
        props1 = props0.replace(']','').replace('[','').strip().split(' ')
        props2 = [float(ele) for ele in props1]
        nnode = Node(node, segments, props2)
        nodes.append(nnode)

    return nodes


def getRecall(realNodes, predictedNodes):
    found = 0
    for node in predictedNodes:
        for realn in realNodes:
            if realn==node:
                print('realn', realn.segments[0].copies)
                print('node', node.segments[0].copies)
                found += 1
                break
    print('num real nodes', len(realNodes))
    print('found', found)
    print('num predicted', len(predictedNodes))
    return found/len(realNodes)



# paction funcs

def snvToIdxDict(snvDf):
    snvToIdx = {}
    for idx, row in snvDf.iterrows():
        snvToIdx[row['snv_mut']] = idx
    return snvToIdx

def segmentTree(currSnvDf, currCnaDf, s, nsamples):
    cnatree, gentrees = processSegment(currSnvDf, currCnaDf)
    currTrees = []
    C = getCnaTree(currCnaDf, cnatree)
    currTrees.append(C)

    G = {} # snv -> gentree
    currSnvIndices = [] # tells us which indices we are processing
    if gentrees==None:
        gentrees = {}
        for snv in list(currSnvDf['snv_mut']):
            gentrees[snv] = [((1,1,0,0),(1,1,0,1))]
    for snv,gentree in gentrees.items():
        G[snv] = findSnvTree(currCnaDf, currSnvDf, snv, gentree)
        currSnvIndices.append(s[snv])
        currTrees.append(G[snv])

    tmpT = solveProgressivePCI(currTrees, nsamples, currCnaDf)
    T = postprocessCombinedTree(tmpT, currSnvIndices, nsamples)
    return T


def postprocessCombinedTree(tmpT, currSnvIndices, nsamples):
    T = nx.DiGraph()
    tmpd = {}
    for node in tmpT.nodes:
        Tnode = (node[0],)
        for i in range(1, len(list(node))):
            Tnode = Tnode + ((currSnvIndices[i-1], node[i][2]+node[i][3]),)
        tmpd[node] = Tnode
    for node in tmpT.nodes:
        T.add_node(tmpd[node])
        for i in range(nsamples):
            T.nodes[tmpd[node]]['sample_' + str(i)] = tmpT.nodes[node]['sample_' + str(i)]
    for edge in tmpT.edges:
        T.add_edge(tmpd[edge[0]], tmpd[edge[1]])
    return T



def getVAF(cnaDf, sample):
    numer = sum((cnaDf['xybarsum'])*cnaDf[sample])
    denom = sum((cnaDf['xysum'])*cnaDf[sample])
    return numer/denom


def findSnvTree(cnaDf, snvDf, snv, gentree):
    snvRow = snvDf.set_index('snv_mut').loc[snv]
    nsamples = cnaDf.columns.str.startswith('sample_').sum()
    cnaDf['cns_tuple'] = cnaDf['copy_number_state'].map(eval)
    cnaDf['xysum'] = cnaDf['cns_tuple'].map(sum)
    cnaDf = cnaDf.set_index('cns_tuple')
    T = nx.DiGraph()
    T.add_edges_from(gentree)
    xy_to_bar = {}
    expanded_xy = None
    expanded_bars = None
    for node in T.nodes:
        xy = node[:2]
        bar = node[2:]
        if xy in xy_to_bar:
            assert expanded_xy is None, "why are there 2 expanded nodes"
            expanded_xy = xy
            expanded_bars = [xy_to_bar[xy], bar]
            xy_to_bar[xy] = (np.nan,np.nan)
        else:
            xy_to_bar[xy] = bar
    assert expanded_xy is not None, "why is there no expanded node"
    cnaDf['xybarsum'] = [sum(xy_to_bar[xy]) for xy in cnaDf.index]
    cnaDf0 = cnaDf.copy(deep=True)
    cnaDf1 = cnaDf.copy(deep=True)
    cnaDf0['xybarsum'].fillna(sum(expanded_bars[0]), inplace=True)
    cnaDf1['xybarsum'].fillna(sum(expanded_bars[1]), inplace=True)    
    for s in range(nsamples):
        sample = f'sample_{s}'
        sample_no_underscore = f'sample{s}'
        vaf = snvRow[sample_no_underscore]
        vaf0 = getVAF(cnaDf0, sample)  # vaf if expanded node is 100% #0
        vaf1 = getVAF(cnaDf1, sample)  # vaf if expanded node is 100% #1
        if vaf0==vaf1 and vaf0==0:
            prop0 = 0
            prop1 = 0
        else:        
            assert vaf0 != vaf1, "cant linear interpolate if it's the same"
            expanded_prop = cnaDf[sample][expanded_xy]
            assert min(vaf0,vaf1) <= vaf and vaf <= max(vaf0,vaf1), "we r gonna have negative proportions"
            prop0 = (vaf1-vaf)/(vaf1-vaf0) * expanded_prop
            prop1 = (vaf-vaf0)/(vaf1-vaf0) * expanded_prop
        for node in T.nodes:
            xy = node[:2]
            if xy == expanded_xy:
                bar = node[2:]
                if bar == expanded_bars[0]:
                    T.nodes[node][sample] = prop0
                elif bar == expanded_bars[1]:
                    T.nodes[node][sample] = prop1
                else:
                    raise Exception("what")
            else:
                T.nodes[node][sample] = cnaDf[sample][xy]
    return T


def getCnaTree(currCnaDf, cnatree):
    C = nx.DiGraph()
    if cnatree==None:
        C.add_nodes_from([(1,1)])
    else:
        C.add_edges_from(cnatree)
    for j,r in currCnaDf.iterrows():
        node = eval(r.copy_number_state)
        for k,v in r.items():
            if 'sample_' in k:
                C.nodes[node][k] = v
    return C

def getCnaTrees(currCnaDf):
    cnaStates = []
    for item in list(currCnaDf['copy_number_state']):
        item = item.replace('(','').replace(')','')
        cnaStates.append(tuple(map(int, item.split(', '))))
    cnaTrees = clonelib.get_cna_trees(set(cnaStates), 1, 1)
    return cnaTrees


def filterByNode(allNodes, currNode):
    ans = []
    for node in allNodes:
        x0,y0,_,_ = node
        if x0==currNode[0] and y0==currNode[1]:
            ans.append(node)
    return ans
        

def drawCnaGraphs(cnaTrees, cnaDf):
    for i in range(len(cnaTrees)):
        cnatree = cnaTrees[i]
        C = nx.DiGraph()
        C.add_edges_from(cnatree)
        labels = {}
        for node in list(C.nodes):
            tmp = cnaDf.loc[cnaDf['copy_number_state']==str(node)]
            props = str(tmp.iloc[0]["copy_number_state"]) + "\n" + str(round(tmp.iloc[0]["sample_0"],3)) + "\n" + str(round(tmp.iloc[0]["sample_1"],3))
            labels[node] = props
        pos = graphviz_layout(C, prog="dot")
        if i == 1:
            nx.draw(C, with_labels=True, font_size=9, node_size=1400, node_color="yellow", labels=labels, pos=pos)
        else:
            nx.draw(C, with_labels=True, font_size=9, node_size=1400, node_color="#ADD8E6", labels=labels, pos=pos)
        plt.savefig("cnaTree" + str(i) + ".png")
        plt.clf()

def drawGenGraph(genTrees, cnaIdx, cnaDf, minVafs, maxVafs, snvDf):

    for j in range(len(genTrees)):
        gentree = genTrees[j]
        G = nx.DiGraph()
        G.add_edges_from(gentree)
        labels = {}
        for longNode in list(G.nodes):
            node = (longNode[0], longNode[1])
            tmp = cnaDf.loc[cnaDf['copy_number_state']==str(node)]
            props = str(tmp.iloc[0]["copy_number_state"]) + "\n" + str(round(tmp.iloc[0]["sample_0"],3)) + "\n" + str(round(tmp.iloc[0]["sample_1"],3))
            labels[longNode] = str(props) + "\n" + str(longNode)
        pos = graphviz_layout(G, prog="dot")
        x = np.array(minVafs[j])
        y = np.array(maxVafs[j])
        x = np.round(x, 3)
        y = np.round(y, 3)
        fits = getFits(minVafs[j], maxVafs[j], snvDf)
        txt = "minvaf " + str(x) + "\n" + "maxvaf " + str(y) + "\n" + "fits " + str(fits) + "\n" + "rows " + str([ele[-8] for ele in list(snvDf['snv_mut'])])
        plt.figtext(0.5, 0.01, txt)

        if cnaIdx == 1:
            nx.draw(G, with_labels=True, font_size=9, node_size=1400, node_color="yellow", labels=labels, pos=pos)
        else:
            nx.draw(G, with_labels=True, font_size=9, node_size=1400, node_color="#ADD8E6", labels=labels, pos=pos)
        plt.savefig("cnaTree" + str(cnaIdx) + "_genTree" + str(j) + ".png")
        plt.clf()


def getFits(minVaf, maxVaf, snvDf):
    n = len(snvDf.columns) - 2
    fits = []
    for snv in list(snvDf['snv_mut']):
        tmp = True
        for i in range(n):
            props = snvDf.loc[snvDf['snv_mut']==snv]['sample'+str(i)]
            vaf = float(props.iloc[0])
            if vaf < minVaf[i] or vaf > maxVaf[i]:
                tmp = False
                break
        fits.append(tmp)
    return fits
            

def processSegment(snvDf, cnaDf):

    cnaTreeAns = None
    genTreesAns = None
    minError = 100

    cnaTrees = getCnaTrees(cnaDf)

    #drawCnaGraphs(cnaTrees, cnaDf)

    for i in range(len(cnaTrees)):
        cnaTree = cnaTrees[i]
        C = nx.DiGraph()
        C.add_edges_from(cnaTree)
        gentrees = clonelib.get_genotype_trees(cnaTree)
        minVafs = []
        maxVafs = []
        currError = 0
        for gentree in gentrees:
            G = nx.DiGraph()
            G.add_edges_from(gentree)
            minVaf, maxVaf = getVafRange(C, G, cnaDf)
            minVafs.append(minVaf)
            maxVafs.append(maxVaf)

        snvgentrees = {}
        for snv in list(snvDf['snv_mut']):
            error, snvgentrees[snv] = assignSnv(snv, snvDf, gentrees, minVafs, maxVafs)
            currError += error
        
        if currError < minError:
            minError = currError
            cnaTreeAns = cnaTree
            genTreesAns = snvgentrees

        #drawGenGraph(gentrees, i, cnaDf, minVafs, maxVafs, snvDf)

    return cnaTreeAns, genTreesAns


def assignSnv(snv, snvDf, gentrees, minVafs, maxVafs):
    minError = 1000
    minTree = None
    for i in range(len(gentrees)):
        error = snvError(snv, snvDf, minVafs[i], maxVafs[i])
        if error < minError:
            minError = error
            minTree = gentrees[i]
    return minError, minTree


def snvError(snv, snvDf, minVaf, maxVaf):
    n = len(snvDf.columns) - 2
    error = 0
    for i in range(n):
        props = snvDf.loc[snvDf['snv_mut']==snv]['sample'+str(i)]
        vaf = float(props.iloc[0])
        if vaf < minVaf[i]:
            error += abs(minVaf[i] - vaf)
        elif vaf > maxVaf[i]:
            error += abs(maxVaf[i] - vaf)
    return error


def calcVaf(nodes, cnaDf):
    n = len(cnaDf.columns) - 2
    numer = np.zeros(n)
    denom = np.zeros(n)
    for node in nodes:
        x0,y0,x1,y1 = node
        propsRow = cnaDf.loc[cnaDf['copy_number_state'] == '(' + str(x0) + ', ' + str(y0) + ')']
        props = [list(propsRow['sample_' + str(m)])[0] for m in range(n)]
        numer += np.array([(x1+y1)*ele for ele in props])
        denom += np.array([(x0+y0)*ele for ele in props])
    return numer/denom


def getVafRange(C, G, cnaDf):
    n = len(cnaDf.columns) - 2
    sectionedNodes = []
    for node in list(C.nodes):
        sectionedNodes.append(filterByNode(list(G.nodes), node))
    combs = [p for p in itertools.product(*sectionedNodes)]
    vafMin = [10 for _ in range(n)]
    vafMax = [-1 for _ in range(n)]
    for comb in combs:
        vaf = calcVaf(list(comb), cnaDf)
        vafMin = [min(vaf[i], vafMin[i]) for i in range(n)]
        vafMax = [max(vaf[i], vafMax[i]) for i in range(n)]
    return vafMin, vafMax


