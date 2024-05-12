import networkx as nx
import numpy as np
from scipy.optimize import minimize, LinearConstraint



def compute_cluster_vaf(S: np.ndarray, T: nx.DiGraph, props: list[float]):
    node_to_prop = dict(zip(sorted(T.nodes), props))
    n_clusters = max(node[1]+1 for node in T.nodes)

    # the node that introduced each cluster
    cluster_start_node = [None]*n_clusters
    for node1, node2 in T.edges:
        if node1[1] != node2[1]:
            cluster_start_node[node2[1]] = node2
    
    apc = nx.all_pairs_node_connectivity(T)
    # list of descendant nodes for each cluster
    cluster_descendants = [
        [node for node in T.nodes if c_start == node or apc[c_start][node]]
        for c_start in cluster_start_node
    ]
    
    cluster_vaf = [[None,None] for _ in range(n_clusters)]
    denom = [
        sum(node[0][p] * node_to_prop[node] for node in T.nodes)
        for p in range(2)
    ]
    for c in range(n_clusters):
        for p in range(2):
            numer = sum(node[0][p] * node_to_prop[node] for node in cluster_descendants[c])
            cluster_vaf[c][p] = numer/denom[p]
    return cluster_vaf



def compute_snv_assignment(S, cluster_vaf, T):
    n_clusters = max(node[1]+1 for node in T.nodes)
    snv_to_cluster = []
    total_err = 0
    for snv_idx, snv_vaf in enumerate(S):
        def d_cluster(cluster_idx):
            return min(np.square(S[snv_idx]-cv) for cv in cluster_vaf[cluster_idx])
        cluster_errors = [d_cluster(c) for c in range(n_clusters)]
        best_cluster = min(range(n_clusters), key=lambda c: cluster_errors[c])
        total_err += cluster_errors[best_cluster]
        snv_to_cluster.append(best_cluster)
    return snv_to_cluster, total_err


def compute_proportions(S,T,snv_to_cluster,old_props):
    n_nodes = len(T.nodes)
    nodes = sorted(T.nodes)
    def err_func(test_props):
        cluster_vafs = compute_cluster_vaf(S, T, test_props)
        err = 0
        for snv_idx, cluster_idx in enumerate(snv_to_cluster):
            err += min(np.square(S[snv_idx]-cv) for cv in cluster_vafs[cluster_idx])
        return err
    res = minimize(
        err_func,
        x0=old_props,
        bounds=[(0,1)]*n_nodes,
        constraints=LinearConstraint(A=np.ones(n_nodes), lb=1, ub=1),
        method='trust-constr',
    )
    assert res.success, "we failed :("
    return list(res.x), err_func(res.x)

def coordinate_descent(S: np.ndarray, T: nx.DiGraph):

    n_nodes = len(T.nodes)
    props = [1/n_nodes] * n_nodes

    err = -1
    max_iter = 100
    num_iter = 0
    while (err == -1 or err >= 1e-10) and num_iter < max_iter:
        n_clusters = max(node[1]+1 for node in T.nodes)
        cluster_vaf = compute_cluster_vaf(S, T, props)
        snv_to_cluster, err = compute_snv_assignment(S,cluster_vaf,T)
        props, err = compute_proportions(S,T,snv_to_cluster, props)
        num_iter += 1

    return err






def main():
    
    T = nx.DiGraph([
        [((1,1),-1), ((1,1),0)],
        [((1,1),0),((1,2),0)],
        [((1,1),0),((1,1),1)],
        [((1,1),1),((1,0),1)],
        [((1,0),1),((1,0),2)],
        [((1,1),1),((1,1),3)],
    ])
    for n_snv in [4,10,100]:
        for i in range(10):
            S = np.random.rand(n_snv)
            err = coordinate_descent(S,T)
            print(n_snv, err, i)

    return

main()