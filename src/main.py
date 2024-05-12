import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import itertools
import logging

from segint import segInt
from treeClass import Tree
from enumerate import Enumerate


# coordinate descent
def optimize_cluster_centers(self, dcfs, clust_group_map): 
        new_dcfs = np.zeros(len(dcfs))

        for q in range(len(dcfs)): #looping over clusters
            # new_dcfs[q] = minimize_scalar(scalar_obj_val, args=(clust_group_map, q, self.data), 
            #                               method='bounded', bounds=[0,1]).x
            result = minimize(objective_function, x0=[dcfs[q]], args=(clust_group_map, q, self.data), bounds=[(0.025,1)])
            if result ==0:
                new_dcfs[q] = self.rng.uniform()
            else:
                new_dcfs[q]= result.x[0]

           
        return new_dcfs

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



def enum_input_mut_cluster_tree(T):
    old_node_names = []
    new_name_to_old_name = {}
    old_name_to_new_name = {}
    for node in T.nodes:
        old_node_names.append(int(node))
    old_node_names = sorted(old_node_names)
    for i in range(len(old_node_names)):
        new_name_to_old_name[i] = old_node_names[i]
        old_name_to_new_name[old_node_names[i]] = i
    newT = nx.DiGraph()
    for node in T.nodes:
        newT.add_node(old_name_to_new_name[int(node)])
    for edge in T.edges:
        u,v = edge
        newT.add_edge(old_name_to_new_name[int(u)], old_name_to_new_name[int(v)])
    return newT


def enum_input_copy_number_state_tree(T):
    newT = nx.DiGraph()
    for node in T.nodes:
        newT.add_node((int(node[1]),int(node[-2])))
    for edge in T.edges:
        u = (int(edge[0][1]),int(edge[0][-2]))
        v = (int(edge[1][1]),int(edge[1][-2]))
        newT.add_edge(u,v)
    return newT




class Params:

    def __init__(self, maxCNstates:list, segments:list, SNVs:list, seeds:list):
        self.maxCNstates = maxCNstates
        self.segments = segments
        self.SNVs = SNVs
        self.seeds = seeds


    def simulate(self) -> None:
        def delete_files_in_directory(directory_path:str) -> None:
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


    def estimate(self) -> None:

        for k,n,s,c in itertools.product(self.segments, self.SNVs, self.seeds, self.maxCNstates):
            
            print('k,n,s,c', k,n,s,c)

            ifile = "".join(["../data/c", str(c), "_k", str(k), "_n", str(n), "_s", str(s), "_T.dot"])
            G = Tree(processClonesim(ifile))
            mut_cluster_tree = enum_input_mut_cluster_tree(G.mut_cluster_tree)

            for kidx in range(G.k):
                copy_number_state_tree = enum_input_copy_number_state_tree(G.segments[kidx]['cnas']['nx'])
                #
                refinements = Enumerate(mut_cluster_tree, copy_number_state_tree).solve()
                # 0. how to get the optimal tree out of all refinements?
                # 1. for each refined tree, for snvs which are in ambiguous positions, use vaf to figure out place
                # 2. to get proportions, use a linear program to minimize correction difference

            # Tdata = segInt(G) ilp is no longer needed

        return

 

def main(maxCNstates:list, segments:list, SNVs:list, seeds:list, createNewSimulations:bool, estimateTree:bool) -> None:

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
    SNVs = [15]
    seeds = [1] #[0,1] + [i for i in range(4,12)] + [i for i in range(13,23)] + [i for i in range(24,33)]
    createNewSimulations = True
    estimateTree = True
    logging.basicConfig(level=logging.DEBUG)
    main(maxCNstates, segments, SNVs, seeds, createNewSimulations, estimateTree)


# input: (hard-coded)
# 1. copy-number states and proportions
# 2. snvs and vafs
# 3. refinement of a) snv mutation cluster tree and b) copy number state tree
# output:
# 1. proportions of tree
# metric:
# 1. obtain the final result of the error function


def foo(C, S, T):
    
    return
