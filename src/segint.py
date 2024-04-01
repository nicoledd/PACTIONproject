import pandas as pd
import networkx as nx
import numpy as np
import itertools
import clonelib


from linear_programs.paction_solver import PCTIsolver
from progressive_caller import solveProgressivePCI, solveProgressivePaction


'''
name: segInt
purpose:
params:
    GTinput: Tree object.
return:
    todo
'''
def segInt(GTinput):

    snvDf = GTinput.snvs
    cnaDf = GTinput.cnas

    s = snvToIdxDict(GTinput.snvs) # snv name -> idx

    trees = []
    cnatrees = []

    for segmentidx in range(GTinput.k):
        currSnvDf = snvDf.loc[snvDf['segment'] == str(segmentidx)]
        currCnaDf0 = cnaDf[segmentidx]
        tmpp = []
        for cna,samples in currCnaDf0.items():
            tmp = [str(cna)]
            for sampleidx, prop in samples.items():
                tmp.append(prop)
            tmpp.append(tmp)
        currCnaDf = pd.DataFrame(tmpp, columns =['copy_number_state'] + ['sample_' + str(ti) for ti in range(GTinput.m)])
        T, cnatree = segmentTree(currSnvDf, currCnaDf, s, GTinput.m)
        trees.append(T)
        cnatrees.append(cnatree)

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

    return finalT, cnatrees



# paction funcs

def snvToIdxDict(snvDf):
    snvToIdx = {}
    for idx, row in snvDf.iterrows():
        snvToIdx[row['snv_mut']] = idx
    return snvToIdx

def segmentTree(currSnvDf, currCnaDf, s, nsamples):
    cnatree, gentrees = processSegment(currSnvDf, currCnaDf, nsamples)
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
    return T, cnatree


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
        if vaf0 == vaf1:
            prop0 = 0
            prop1 = 0
        else:
            expanded_prop = cnaDf[sample][expanded_xy]
            vaf = float(vaf)
            vaf1 = float(vaf1)
            vaf0 = float(vaf0)
            if vaf < 0:
                prop0 = 0
                prop1 = 0
            else:
                prop0 = (vaf1-vaf)/(vaf1-vaf0) * expanded_prop
                prop1 = (vaf-vaf0)/(vaf1-vaf0) * expanded_prop
        for node in T.nodes:
            xy = node[:2]
            if xy == expanded_xy:
                bar = node[2:]
                if bar == expanded_bars[0]:
                    prop0 = max(prop0, 0)
                    T.nodes[node][sample] = prop0
                elif bar == expanded_bars[1]:
                    prop1 = max(prop0, 0)
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


def processSegment(snvDf, cnaDf, n):

    cnaTreeAns = None
    genTreesAns = None
    minError = 100

    cnaTrees = getCnaTrees(cnaDf)

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
            minVaf, maxVaf = getVafRange(C, G, cnaDf, n)
            minVafs.append(minVaf)
            maxVafs.append(maxVaf)

        snvgentrees = {}
        for snv in list(snvDf['snv_mut']):
            error, snvgentrees[snv] = assignSnv(snv, snvDf, gentrees, minVafs, maxVafs, n)
            currError += error
        
        if currError < minError:
            minError = currError
            cnaTreeAns = cnaTree
            genTreesAns = snvgentrees

    return cnaTreeAns, genTreesAns


def assignSnv(snv, snvDf, gentrees, minVafs, maxVafs, n):
    minError = 1000
    minTree = None
    for i in range(len(gentrees)):
        error = snvError(snv, snvDf, minVafs[i], maxVafs[i], n)
        if error < minError:
            minError = error
            minTree = gentrees[i]
    return minError, minTree


def snvError(snv, snvDf, minVaf, maxVaf, n):
    error = 0
    for i in range(n):
        props = snvDf.loc[snvDf['snv_mut']==snv]['sample'+str(i)]
        vaf = float(props.iloc[0])
        if vaf < minVaf[i]:
            error += abs(minVaf[i] - vaf)
        elif vaf > maxVaf[i]:
            error += abs(maxVaf[i] - vaf)
    return error


def calcVaf(nodes, cnaDf, n):
    numer = np.zeros(n)
    denom = np.zeros(n)
    for node in nodes:
        x0,y0,x1,y1 = node
        propsRow = cnaDf.loc[cnaDf['copy_number_state'] == '(' + str(x0) + ', ' + str(y0) + ')']
        props = [list(propsRow['sample_' + str(m)])[0] for m in range(n)]
        numer += np.array([(x1+y1)*ele for ele in props])
        denom += np.array([(x0+y0)*ele for ele in props])
    return numer/denom


def getVafRange(C, G, cnaDf, n):
    sectionedNodes = []
    for node in list(C.nodes):
        sectionedNodes.append(filterByNode(list(G.nodes), node))
    combs = [p for p in itertools.product(*sectionedNodes)]
    vafMin = [10 for _ in range(n)]
    vafMax = [-1 for _ in range(n)]
    for comb in combs:
        vaf = calcVaf(list(comb), cnaDf, n)
        vafMin = [min(vaf[i], vafMin[i]) for i in range(n)]
        vafMax = [max(vaf[i], vafMax[i]) for i in range(n)]
    return vafMin, vafMax


