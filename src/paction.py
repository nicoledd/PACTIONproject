import sys
import argparse
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt

from algo import runAlgo
from GTtoData import getGTdata
from GTtoInput import getGTinput
from compare import runCompare

def get_relations(directed_tree, snv_introductions):
    isancestor = nx.all_pairs_node_connectivity(directed_tree)
    ancestries = []
    clusters = []
    unrelateds = []
    for n0 in directed_tree.nodes:
        clusters += [
            (snv1, snv2)
            for snv1 in snv_introductions[n0]
            for snv2 in snv_introductions[n0]
            if snv1 < snv2
        ]
    for n1 in directed_tree.nodes:
        for n2 in directed_tree.nodes:
            if n1 >= n2:  # i hope inequality operators are well defined on these nodes...
                continue
            if isancestor[n1][n2]:
                ancestries += [
                    (snv1,snv2)
                    for snv1 in snv_introductions[n1] 
                    for snv2 in snv_introductions[n2]
                ]
            elif isancestor[n2][n1]:
                ancestries += [
                    (snv2,snv1)
                    for snv1 in snv_introductions[n1] 
                    for snv2 in snv_introductions[n2]
                ]
            else:
                unrelateds += [
                    (min(snv1, snv2), max(snv1, snv2))
                    for snv1 in snv_introductions[n1] 
                    for snv2 in snv_introductions[n2]
                ]
    return unrelateds, clusters, ancestries

def get_relations_predicted(Tdata):
    def get_snvs(node):
        ans = []
        segments = [[]]
        for x in node:
            if isinstance(x[0], int):
                segments[0].append(x)
            else:
                segments.append(x)
        for seg in segments:
            for tup in seg[1:]:
                if tup[1] > 0:
                    ans.append(tup[0])
        return ans
    snv_per_node = {node: set(get_snvs(node)) for node in Tdata.nodes}
    snv_introductions = {node: set() for node in Tdata.nodes}
    for n1,n2 in Tdata.edges:
        assert len(snv_introductions[n2]) == 0
        snv_introductions[n2] = snv_per_node[n2]-snv_per_node[n1]
    return get_relations(Tdata, snv_introductions)

def get_relations_true(GTinput):
    snv_introductions = defaultdict(list)
    for v, muts in GTinput.snvs.groupby('v')['snv_mut']:
        snv_introductions[v] = list(map(int,muts))

    GTInputGDirected = nx.DiGraph()
    GTInputGDirected.add_edges_from(GTinput.G.edges)
    return get_relations(GTInputGDirected, snv_introductions)

def get_recall(true_values, pred_values):
    common = set(true_values) & set(pred_values)
    return len(common) / len(true_values) if len(true_values) else 0

def get_tree_recalls(GTinput, Tdata):
    true_unrelateds, true_clusters, true_ancestries = get_relations_true(GTinput)
    pred_unrelateds, pred_clusters, pred_ancestries = get_relations_predicted(Tdata)
    ipr = get_recall(true_unrelateds, pred_unrelateds)
    cpr = get_recall(true_clusters, pred_clusters)
    apr = get_recall(true_ancestries, pred_ancestries)
    return ipr, cpr, apr

def main():

    print('k, n, s, predicted, true')
    iprs = []
    cprs = []
    aprs = []
    for k in ['10', '100']: # 1. for-loops are just parameters for the type of tree
        for n in ['10', '50', '100']:
            for s in ['1', '2', '3']:
                g = "../data/k" + k + "_n" + n + "_s" + s + "/T.dot" # 2. g: name of tree file
                o = "../results/"
                #GTdata = getGTdata(g)
                GTinput = getGTinput(g) # 3. GTinput: object that represents the ground truth tree
                Tdata, usedIndices = runAlgo(GTinput) # 4. Tdata: either an object or nx graph that represents the estimated tree, ignore usedIndices

                # compare: number of nodes (x)
                # compare: cna states of nodes
                
                # 5. below here: code that computes the recall
                num = 0
                trueCnas = []

                for node in GTinput.G.nodes:
                    rows = GTinput.G.nodes[node]['label'].split('\\n')
                    cnasofnode = rows[0].split(';')
                    cnasofnode = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in cnasofnode]
                    cnasofnode = cnasofnode[:-1]
                    cnasofnode = [tuple(map(int, ele.split(','))) for ele in cnasofnode]
                    trueCnas.append(cnasofnode)

                for node in Tdata.nodes:
                    currCnas = []
                    lines = Tdata.nodes[node]['realName'].split('\n')
                    currIdx = 0

                    for segmentidx in range(GTinput.k):
                        if segmentidx not in usedIndices or currIdx > len(lines)-1:
                            currCnas.append((1,1))
                        else:
                            if lines[currIdx] == "":
                                currCnas.append((1,1))
                            else:
                                currCnas.append((int(lines[currIdx][2]), int(lines[currIdx][5])))
                            currIdx += 1

                    # if there's any node in GTinput.cnas that look like this, plus one
                    if currCnas in trueCnas:
                        num += 1

                # 6. todo: write code that computes the ancestral recall for true tree vs estimated tree


                # GTinput, Tdata
                ipr, cpr, apr = get_tree_recalls(GTinput, Tdata)
                print("k,n,s", k,n,s)
                print(f"{ipr=} {cpr=} {apr=}")
                iprs.append(ipr)
                cprs.append(cpr)
                aprs.append(apr)
    plt.figure()
    plt.hist(iprs)
    plt.xlabel('ipr (in distinct paths, recall)')
    plt.ylabel('count')
    plt.savefig('ipr.jpg')
    plt.close()

    plt.figure()
    plt.hist(cprs)
    plt.xlabel('cpr (clustered, recall)')
    plt.ylabel('count')
    plt.savefig('cpr.jpg')
    plt.close()

    plt.figure()
    plt.hist(aprs)
    plt.xlabel('apr (ancestral, recall)')
    plt.ylabel('count')
    plt.savefig('apr.jpg')
    plt.close()

    

                


main()



