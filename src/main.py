import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import itertools
import logging

from segint import segInt
from treeClass import Tree


# recall functions

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
    GTInputGDirected.add_edges_from(GTinput.T.edges)
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




# include in other function
def processClonesim(f):
    def computeCNstates(rows):
        CNstates = rows[0].split(';')
        CNstates = [ele.strip(' ').strip('\"').strip('(').strip(')') for ele in CNstates]
        CNstates = CNstates[:-1]
        CNstates = [tuple(map(int, ele.split(','))) for ele in CNstates]
        return CNstates
    def computeSNVs(rows):
        SNVs = []
        for segmentIdx in range(len(rows)):
            if rows[segmentIdx] == "":
                SNVs.append([])
            elif rows[segmentIdx][0] in [str(i) for i in range(10)]:
                currSnvs = rows[segmentIdx].split(';')
                currSnvs = [ele.strip() for ele in currSnvs]
                currSnvs = [ele for ele in currSnvs if ele != ""]
                currSnvs = [ele.split('|') for ele in currSnvs]
                for i in range(len(currSnvs)):
                    for j in range(len(currSnvs[i])):
                        currSnvs[i][j] = int(currSnvs[i][j])
                SNVs.append(currSnvs)
        return SNVs
    def computek(T):
        for node in T.nodes:
            rows = T.nodes[node]['label'].split('\\n')
            k = len(rows[0].split(';'))-1
            return k
    def computeNodes(G, T):
        k = computek(G)
        for node in G.nodes:
            rows = G.nodes[node]['label'].split('\\n')
            CNstates = computeCNstates(rows)
            SNVs = computeSNVs(rows)
            T.nodes[node]['props'] = [float(ele) for ele in rows[-1].split(" ") if ele[0]=='0']
            for kidx in range(k):
                if 'segments' not in T.nodes[node]:
                    T.nodes[node]['segments'] = {}
                T.nodes[node]['segments'][kidx] = {}
                T.nodes[node]['segments'][kidx]['cna'] = CNstates[kidx]
                T.nodes[node]['segments'][kidx]['snvs'] = SNVs[kidx]
        return T
    R = nx.DiGraph(nx.nx_pydot.read_dot(f))
    T = nx.DiGraph()
    T.add_nodes_from(R.nodes)
    T.add_edges_from(R.edges)
    for node in T.nodes:
        T.nodes[node]['name'] = node
    T = computeNodes(R, T)
    return T





'''
class: Params
attributes:
    maxCNstates: list. max numbers of copy-number states allowed in simulation
    segments: list. number of segments in simulation
    SNVs: list. number of SNVs in simulation
    seeds: list. seed values in simulation
methods:
    simulate
    estimate
'''
class Params:

    def __init__(self, maxCNstates, segments, SNVs, seeds):
        self.maxCNstates = maxCNstates
        self.segments = segments
        self.SNVs = SNVs
        self.seeds = seeds

    '''
    name: simulate
    purpose: run a series of terminal commands that generate simulated-tree-files in a director
    params:
        None
    return:
        None
    '''
    def simulate(self):
        def delete_files_in_directory(directory_path):
            try:
                files = os.listdir(directory_path)
                for file in files:
                    file_path = os.path.join(directory_path, file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
            except OSError:
                print("Error occurred while deleting files.")
        delete_files_in_directory('../data/')
        for i in self.maxCNstates:
            os.system('../clonesim/build/generatecnatrees -maxCN ' + str(i) + ' > ../data/cnatrees' + str(i) + '.txt')
        for k,n,s,c in itertools.product(self.segments, self.SNVs, self.seeds, self.maxCNstates):
            cmd = " ".join(['../clonesim/build/simulate -r -s', str(s), '-k', str(k), '-n', str(n), '-S ../data/cnatrees' + str(i) + '.txt -dot ../data/c'+str(c)+'_k'+str(k)+'_n'+str(n)+'_s'+str(s)+'_T.dot > ../data/o.txt'])
            os.system(cmd)
    

    def decifer(self, G):
        decifer_input_file = 'decifer_input.tsv'
        file = open(decifer_input_file, 'w')
        file.write('0 #characters\n0 #samples\n#sample_index	sample_label	character_index	character_label	ref	var\n')
        file.close()

        sample = []
        ref = []
        var = []
        cna_major = []
        cna_minor = []
        cna_props = []

        # go through graph by node
        for node in G.T.nodes:
            for kidx in range(G.k):
                curr_cna_major, curr_cna_minor = G.T.nodes[node]['segments'][kidx]['cna']
                for snvidx in range(G.segments[kidx]['n']):
                    for sample_idx in range(G.m):
                        curr_var = G.segments[kidx]['snvs'][snvidx]['mut'][sample_idx]
                        curr_ref = G.segments[kidx]['snvs'][snvidx]['ref'][sample_idx]
                        sample.append(sample_idx)
                        cna_major.append(curr_cna_major)
                        cna_minor.append(curr_cna_minor)
                        ref.append(curr_ref)
                        var.append(curr_var)
                        node = G.segments[kidx]['snvs'][snvidx]['node']
                        cna_state_of_node = G.T.nodes[node]['segments'][kidx]['cna']
                        df = G.segments[kidx]['cnas']['df']
                        df2 = df[df['copy_number_state'] == str(cna_state_of_node)]
                        val = df2.iloc[0]['sample_' + str(sample_idx)]
                        cna_props.append(val)

        mut = [i for i in range(len(sample))]

        file = open(decifer_input_file, 'a')
        for i in range(len(sample)):
            file.write('\t'.join(str(ele) for ele in [sample[i],sample[i],mut[i],mut[i],ref[i],var[i],cna_major[i],cna_minor[i],cna_props[i]]))
            file.write('\n')
        file.close()


    '''
    name: estimate
    purpose:
        simulated tree - algorithm -> estimated tree
        create metrics plots
    params:
        None
    return:
        None
    '''
    def estimate(self):

        def cnaTreeRecall(G, cnatrees, k):
            pred = []
            true = []
            for tree in cnatrees:
                T = nx.DiGraph()
                if tree == None:
                    continue
                for edge in tree:
                    u,v = edge
                    T.add_edge(str(u), str(v))
                pred.append(T)
            for kidx in range(k):
                true.append(G.segments[kidx]['cnas']['nx'])
            trueIndices = [0 for _ in range(len(true))]
            for treeP in pred:
                trueIdx = 0
                for treeT in true:
                    if set(treeP.edges) == set(treeT.edges):
                        trueIndices[trueIdx] = 1
                    trueIdx += 1
            return sum(trueIndices)/k
    
        def plotHist(vals, xlabel, ylabel, ofile):
            plt.figure()
            plt.hist(vals)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.savefig(ofile)
            plt.close()

        iprs = []
        cprs = []
        aprs = []

        numNodesG = []
        numNodesT = []
            
        for k,n,s,c in itertools.product(self.segments, self.SNVs, self.seeds, self.maxCNstates):
            
            print('k,n,s,c', k,n,s,c)

            ifile = "".join(["../data/c", str(c), "_k", str(k), "_n", str(n), "_s", str(s), "_T.dot"])
            G = Tree(processClonesim(ifile))
            self.decifer(G)
            Tdata = segInt(G)
            ipr, cpr, apr = get_tree_recalls(G, Tdata)
            numnodest = len(Tdata)
            numnodesg = len(G.T.nodes)
            numNodesG.append(numnodesg)
            numNodesT.append(numnodest)
            iprs.append(ipr)
            cprs.append(cpr)
            aprs.append(apr)
        return

 

'''
name: main
purpose:
    Create new simulations if bool is set.
    Run experiment (simulated tree - algorithm -> predicted tree) if bool is set.
params:
    maxCNstates: list. max numbers of copy-number states allowed in simulation
    segments: list. number of segments in simulation
    SNVs: list. number of SNVs in simulation
    seeds: list. seed values in simulation
    createNewSimulations: bool. whether to create new simulations or not
    estimateTree: bool. whether to run experiment or not
return:
    None
'''
def main(maxCNstates, segments, SNVs, seeds, createNewSimulations, estimateTree):
    # 1. get decifer input
    # 2. run decifer
    # 3. obtain decifer output
    
    def createSimulationsCaller(createNewSimulations, pObj):
        if createNewSimulations:
            pObj.simulate()

    def estimateTreeCaller(estimateTree, pObj):
        if estimateTree:
            pObj.estimate()

    pObj = Params(maxCNstates, segments, SNVs, seeds)
    createSimulationsCaller(createNewSimulations, pObj)
    estimateTreeCaller(estimateTree, pObj)



if __name__ == "__main__":
    maxCNstates = [2]
    segments = [5]
    SNVs = [20]
    seeds = [1] #[0,1] + [i for i in range(4,12)] + [i for i in range(13,23)] + [i for i in range(24,33)]
    createNewSimulations = True
    estimateTree = True
    logging.basicConfig(level=logging.DEBUG)
    main(maxCNstates, segments, SNVs, seeds, createNewSimulations, estimateTree)

