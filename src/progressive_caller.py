import csv
from csv import writer
import pandas as pd
import networkx as nx
import numpy as np
from collections import deque

from linear_programs.paction_solver import PCIsolver, PCTIsolver
from preprocess.clean_data import getIndexOfTrueGraphs, processSolutionTree, processTreeFile, ifTrueTreesAreFoundIn, getTreeEdges, readProportionMatrices, getAllPossibleTrees
from postprocess.write import writeText, writeTree, writeClones
from modality.modality_class import Modality


def solveProgressivePCI(trees, nsamples, cnaDf):
    vafs = {}
    nsamples = cnaDf.columns.str.startswith('sample').sum()
    vafs[0] = [-1 for _ in range(nsamples)]
    for i in range(1, len(trees)):
        vaf = obtainVaf(trees[i], cnaDf, nsamples)
        vafs[i] = vaf

    if trees == []:
        return nx.DiGraph()
    treeA = trees[0]

    for i in range(1, len(trees)):
        treeB = trees[i]
        tupleToIdA, idToTupleA = getProcessedClones(list(nx.topological_sort(treeA)))
        tupleToIdB, idToTupleB = getProcessedClones(list(nx.topological_sort(treeB)))
        treeAFake = translateTo(treeA, tupleToIdA, nsamples, False)
        treeBFake = translateTo(treeB, tupleToIdB, nsamples, True)
        currSol = solvePCI(treeAFake, treeBFake, treeA, treeB, nsamples, vafs[i], cnaDf, tupleToIdA, tupleToIdB)
        treeA = translateFrom(currSol.G, idToTupleA, idToTupleB, nsamples)
        treeB = None
   
    return treeA


def obtainVaf(tree, cnaDf, nsamples):
    numer = [0 for _ in range(nsamples)]
    denom = [0 for _ in range(nsamples)]
    for node in tree.nodes:
        for i in range(nsamples):
            prop = cnaDf.loc[cnaDf['copy_number_state']==str((node[0],node[1]))]['sample_'+str(i)].iloc[0]
            numer[i] += prop*(node[2]+node[3])
            denom[i] += prop*(node[0]+node[1])
    vaf = np.array(numer)/np.array(denom)
    return vaf


def getProcessedClones(trueClones):
    tupleToName = {}
    nameToTuple = {}
    for i in range(len(trueClones)):
        tupleToName[trueClones[i]] = i
        nameToTuple[i] = trueClones[i]
    return tupleToName,nameToTuple

def translateTo(tree, d, nsamples, treeBOrNot):
    T = nx.DiGraph()
    for u,v in list(tree.edges):
        T.add_edge(d[u], d[v])
    for u in list(tree.nodes):
        T.add_node(d[u])
        for i in range(nsamples):
            T.nodes[d[u]]['sample_'+str(i)] = tree.nodes[u]['sample_'+str(i)]
        if treeBOrNot == False:
            if type(u[0]) is int:
                T.nodes[d[u]]['x'] = u[0]
                T.nodes[d[u]]['y'] = u[1]
            else:
                T.nodes[d[u]]['x'] = u[0][0]
                T.nodes[d[u]]['y'] = u[0][1]
        else:
            T.nodes[d[u]]['x'] = u[0]
            T.nodes[d[u]]['y'] = u[1]
            T.nodes[d[u]]['xbar'] = u[2]
            T.nodes[d[u]]['ybar'] = u[3]
    return T


def translateFrom(tree, e, e2, nsamples):
    T = nx.DiGraph()
    for uv,wx in list(tree.edges):
        u,v = uv
        w,x = wx
        edge1 = (e[u],e2[v])
        edge2 = (e[w],e2[x])
        if type(edge1[0][0]) is not int:
            edge1 = edge1[0] + edge1[1:]
            edge2 = edge2[0] + edge2[1:]
        T.add_edge(edge1, edge2)
    for u,v in list(tree.nodes):
        node = (e[u],e2[v])
        if type(node[0][0]) is not int:
            node = node[0] + node[1:]
        T.add_node(node)
        for i in range(nsamples):
            T.nodes[node]['sample_'+str(i)] = tree.nodes[(u,v)]['sample_'+str(i)]
    return T



def solvePCI(Afake, Bfake, A, B, nsamples, vafs, cnaDf, tupleToIdA, tupleToIdB):
    T = refine(A,B)
    G = nx.DiGraph()
    for edge in T.edges:
        u,v = edge
        au = u[:-1]
        av = v[:-1]
        bu = u[-1]
        bv = v[-1]
        if len(au)==1:
            G.add_edge((tupleToIdA[au[0]],) + (tupleToIdB[bu],), (tupleToIdA[av[0]],) + (tupleToIdB[bv],))
        else:
            G.add_edge((tupleToIdA[au],) + (tupleToIdB[bu],), (tupleToIdA[av],) + (tupleToIdB[bv],))
    solver = PCIsolver(Afake, Bfake, set(G.nodes), list(G.edges), nsamples, vafs, cnaDf)
    solver.solve()
    return solver


def refine(A,B):
    if isinstance(list(A.nodes)[0][0], int):
        G = nx.relabel_nodes(A, {n: (n,) for n in A.nodes})
    else:
        G = A.copy()

    splitedge = None
    splitcna = None
    for be in B.edges:
        if be[0][:2] == be[1][:2]:
            assert splitedge is None
            splitedge = be
            splitcna = be[0][:2]

    cna_to_new_node = {splitcna: splitedge[1]}
    for bn in B.nodes:
        if bn[:2] != splitcna:
            cna_to_new_node[bn[:2]] = bn

    splitnode = None
    for anode in nx.topological_sort(G):
        if tuple(anode[0]) == tuple(splitcna):
            splitnode = anode
            break

    splitnode0 = splitnode + (splitedge[0],)
    if G.in_degree(splitnode) > 0:
        assert G.in_degree(splitnode) == 1
        par = list(G.in_edges(splitnode))[0][0]
        G.remove_edge(par, splitnode)
        G.add_edge(par, splitnode0)
    G.add_edge(splitnode0, splitnode)

    nx.relabel_nodes(G, {
        n: n if n == splitnode0 else n + (cna_to_new_node[n[0][:2]],)
    for n in G.nodes}, copy=False)

    return G


def findSpanningTree(A,B):
    G = nx.DiGraph()
    T = nx.DiGraph()
    for na in A.nodes:
        for nb in B.nodes:
            if A.nodes[na]['x']==B.nodes[nb]['x'] and A.nodes[na]['y']==B.nodes[nb]['y']:
                T.add_node((na,nb))
    for i in range(len(T.nodes)):
        n1 = list(T.nodes)[i]
        a,b = n1
        for j in range(i+1,len(T.nodes)):
            n2 = list(T.nodes)[j]
            c,d = n2
            if ((a,c) in A.edges or a==c) and ((b,d) in B.edges or b==d):
                T.add_edge(n1,n2)
    # get minimum spanning tree
    newT = nx.minimum_spanning_tree(T.to_undirected())
    E = set(newT.edges())
    newE = [e for e in T.edges() if e in E or reversed(e) in E]
    G.add_edges_from(newE)
    return G


def spanningTrees(self,G):
    def build_tree(H, edges):
        if nx.is_connected(H):
            yield H
        else:
            for i in range(len(edges)):
                if edges[i][1] not in nx.algorithms.descendants(H, edges[i][0]):
                    H1 = nx.Graph(H)
                    H1.add_edge(*edges[i])
                    for H2 in build_tree(H1, edges[i+1:]):
                        yield H2
    E = nx.Graph()
    E.add_nodes_from(G)
    return build_tree(E, [e for e in G.edges])



def solveProgressivePaction(trees, nsamples):
    tupleToId = {}
    idToTuple = {}
    if trees == []:
        return nx.DiGraph()
    treeA = trees[0]
    for i in range(1, len(trees)):
        treeB = trees[i]
        tupleToId, idToTuple = getProcessedClones(list(nx.topological_sort(treeA)))
        tupleToIdB, idToTupleB = getProcessedClones(list(nx.topological_sort(treeB)))
        treeAFake = translateToFinal(treeA, tupleToId, nsamples)
        treeBFake = translateToFinal(treeB, tupleToIdB, nsamples)
        currSol = solvePaction(treeAFake, treeBFake, nsamples)
        treeA = translateFrom(currSol.G, idToTuple, idToTupleB, nsamples)
        treeB = None
    return treeA

def translateToFinal(tree, d, nsamples):
    T = nx.DiGraph()
    for u,v in list(tree.edges):
        T.add_edge(d[u], d[v])
    for u in list(tree.nodes):
        T.add_node(d[u])
        for i in range(nsamples):
            T.nodes[d[u]]['sample_'+str(i)] = tree.nodes[u]['sample_'+str(i)]
    return T


def solvePaction(treeA, treeB, nsamples):
    solver = PCTIsolver(treeA, treeB, nsamples)
    solver.solve()
    return solver