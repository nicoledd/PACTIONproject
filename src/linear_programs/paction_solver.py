import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools

class PCIsolver:
    def __init__(self, treeA, treeB, validNodes, validEdges, nsamples, vafs, cnaDf):
        self.n = int(nsamples)
        self.A = treeA
        self.B = treeB
        self.G = nx.DiGraph()
        self.G.add_nodes_from(validNodes)
        self.G.add_edges_from(validEdges)
        self.validNodes = validNodes

        self.correction = None
        self.numSol = None
        self.vafs = vafs
        self.cnaDf = cnaDf

    def solve(self):

        # 1) initialize variables
        nsamples = self.n
        nsnv = len(list(self.A.nodes))
        ncna = len(list(self.B.nodes))

        # 2) create linear programming model and set flags
        model = gp.Model('PCTIsolver')
        model.setParam('OutputFlag',0)
        model.setParam('SolutionLimit',10)

        # 3) initialize linear programming variables
        # x: clones in the final set
        # w, y: the proportion of each clone
        # d: correction variables
        x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
        w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
        d_snv = model.addVars(nsamples, nsnv, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_snv')
        d_cna = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_cna')
        f = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'vaf')


        # 4) Add constraints

        # 4.1) Clone Existance Constraints: encode for the predicted clones
        # encode product w[i,j,k] = y[i,j,k] * x[j,k]
        for i in range(nsamples):
            for j in range(nsnv):
                for k in range(ncna):
                    model.addConstr(w[i,j,k] <= y[i,j,k])
                    model.addConstr(w[i,j,k] <= x[j,k])
                    model.addConstr(w[i,j,k] >= x[j,k] + y[i,j,k] - 1)
        for k in range(ncna):
            for j in range(nsnv):
                if (j,k) in self.validNodes:
                    model.addConstr(x[j,k]==1)
                else:
                    model.addConstr(x[j,k]==0)

        # Correction Constraints: VAF
        # this is the computed vaf
        for i in range(nsamples):
            numerSumm = 0
            denomSumm = 0
            for k in range(ncna):
                for j in range(nsnv):
                    numerSumm += x[j,k]*(self.B.nodes[k]['xbar'] + self.B.nodes[k]['ybar'])*self.B.nodes[k]['sample_'+str(i)]
                    denomSumm += x[j,k]*(self.B.nodes[k]['x'] + self.B.nodes[k]['y'])*self.B.nodes[k]['sample_'+str(i)]
                model.addConstr(denomSumm*self.vafs[i] - numerSumm <= denomSumm*f[i,k])
                model.addConstr(numerSumm - denomSumm*self.vafs[i] <= denomSumm*f[i,k])

        # 4.2) Correction Constraints: difference between input and output matrices
        for i in range(nsamples): # constraint for snv
            for j in range(nsnv):
                summ = gp.LinExpr()
                for k in range(ncna):
                    summ += w[i,j,k]
                model.addConstr(self.A.nodes[j]['sample_'+str(i)] - summ <= d_snv[i,j])
                model.addConstr(summ - self.A.nodes[j]['sample_'+str(i)] <= d_snv[i,j])
        for i in range(nsamples): # constraint for cna
            for k in range(ncna):
                sum = gp.LinExpr()
                for j in range(nsnv):
                    sum += w[i,j,k]
                model.addConstr(self.B.nodes[k]['sample_'+str(i)] - sum <= d_cna[i,k])
                model.addConstr(sum - self.B.nodes[k]['sample_'+str(i)] <= d_cna[i,k])

        # 4.4) Abundance Constraint: Sum of each proportion matrix row must be 1
        for i in range(nsamples):
            summ = gp.LinExpr()
            for j in range(nsnv):
                for k in range(ncna):
                    summ += w[i,j,k]
            model.addConstr(summ == 1)

        # 5) Objective Function: Minimize the correction
        obj_sum = gp.LinExpr()
        for i in range(nsamples):
            for j in range(nsnv):
                obj_sum += d_snv[i,j]
            for k in range(ncna):
                obj_sum += d_cna[i,k] + f[i,k]
        model.setObjective(obj_sum, gp.GRB.MINIMIZE)

        # 6) Run the Model and evaluate

        # 6.2) Run the model
        model.optimize()

        # 6.3) If infeasible, return
        if model.status != gp.GRB.OPTIMAL:
            return
        
        # 6.4) If feasible, store solutions
        else:
            self.processSolution(model, x, w)


    def processSolution(self, model, x, w):
        self.correction = model.objVal
        sol_x = model.getAttr('x', x)
        clones = [key for key, val in sol_x.items() if val >= 0.5]
        proportions = model.getAttr('x', w)
        for clone in clones:
            self.G.add_node(clone)
            for i in range(self.n):
                prop = proportions[i,clone[0],clone[1]]
                if prop < 0.001:
                    prop = 0
                self.G.nodes[clone]['sample_'+str(i)] = prop
        self.numSol = model.SolCount




 


class PCTIsolver:
    def __init__(self, treeA, treeB, nsamples):
        self.A = treeA
        self.B = treeB
        self.n = int(nsamples)
        self.G = nx.DiGraph()
        self.correction = None
        self.numSol = None

    # create dictionaries which store parent of each node
    def createParentDict(self, G):
        d = {}
        for edge in G.edges:
            child = edge[1]
            parent = edge[0]
            if child not in d.keys():
                d[child] = [parent]
            else:
                d[child].append(parent)
        return d

    def solve(self):

        # 1) initialize variables
        nsamples = self.n
        nsnv = len(list(self.A.nodes))
        ncna = len(list(self.B.nodes))

        # 2) create linear programming model and set flags
        model = gp.Model('PCTIsolver')
        model.setParam('OutputFlag',0)
        model.setParam('SolutionLimit',10)

        # 3) initialize linear programming variables
        # x: clones in the final set
        # w, y: the proportion of each clone
        # z: the PREDICTED edges in the cna and snv trees
        # d: correction variables
        x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
        w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
        z_snv = model.addVars(nsnv-1, ncna, vtype = gp.GRB.BINARY, lb = 0, ub = 1, name = 'z_snv')
        z_cna = model.addVars(nsnv, ncna-1, vtype = gp.GRB.BINARY, lb = 0, ub = 1, name = 'z_cna')
        d_snv = model.addVars(nsamples, nsnv, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_snv')
        d_cna = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_cna')

        # 4) Add constraints

        # 4.1) Clone Existance Constraints: encode for the predicted clones
        # encode product w[i,j,k] = y[i,j,k] * x[j,k]
        for i in range(nsamples):
            for j in range(nsnv):
                for k in range(ncna):
                    model.addConstr(w[i,j,k] <= y[i,j,k])
                    model.addConstr(w[i,j,k] <= x[j,k])
                    model.addConstr(w[i,j,k] >= x[j,k] + y[i,j,k] - 1)

        # 4.2) Correction Constraints: difference between input and output matrices
        for i in range(nsamples): # constraint for snv
            for j in range(nsnv):
                summ = gp.LinExpr()
                for k in range(ncna):
                    summ += w[i,j,k]
                model.addConstr(self.A.nodes[j]['sample_'+str(i)] - summ <= d_snv[i,j])
                model.addConstr(summ - self.A.nodes[j]['sample_'+str(i)] <= d_snv[i,j])
        for i in range(nsamples): # constraint for cna
            for k in range(ncna):
                sum = gp.LinExpr()
                for j in range(nsnv):
                    sum += w[i,j,k]
                model.addConstr(self.B.nodes[k]['sample_'+str(i)] - sum <= d_cna[i,k])
                model.addConstr(sum - self.B.nodes[k]['sample_'+str(i)] <= d_cna[i,k])

        # 4.3) Edge Existance Constraints: edge exists <=> endpoints exist
        # encode z_snv[j,k] = x[parent(j), k] * x[j,k]
        # encode z_cna[j,k] = x[j, parent(k)] * x[j,k]
        for edge_idx, edge in enumerate(list(self.A.edges)):
            parent = edge[0] or edge[1]
            child = edge[1] or edge[0]
            for k in range(ncna):
                model.addConstr(z_snv[edge_idx, k] <= x[parent, k])
                model.addConstr(z_snv[edge_idx, k] <= x[child, k])
                model.addConstr(z_snv[edge_idx, k] >= x[parent, k] + x[child, k] - 1)
        for edge_idx, edge in enumerate(list(self.B.edges)):
            parent = edge[0] or edge[1]
            child = edge[1] or edge[0]
            for j in range(nsnv):
                model.addConstr(z_cna[j, edge_idx] <= x[j, parent])
                model.addConstr(z_cna[j, edge_idx] <= x[j, child])
                model.addConstr(z_cna[j, edge_idx] >= x[j, parent] + x[j, child] - 1)

        # 4.4) Abundance Constraint: Sum of each proportion matrix row must be 1
        for i in range(nsamples):
            summ = gp.LinExpr()
            for j in range(nsnv):
                for k in range(ncna):
                    summ += w[i,j,k]
            model.addConstr(summ == 1)

        # 4.5) summation of edges constraints: There exists only one edge of each type
        # encode sum_{k} z_snv[j, k] == 1
        # encode sum_{j} z_cna[j, k] == 1
        for uv in range(nsnv-1):
            summ = gp.LinExpr()
            for k in range(ncna):
                summ += z_snv[uv, k]
            model.addConstr(summ == 1)
        for uv in range(ncna - 1):
            summ = gp.LinExpr()
            for j in range(nsnv):
                summ += z_cna[j, uv]
            model.addConstr(summ == 1)

        # 4.6) No-Orphans Constraints: orphan clones are not allowed
        # encode x[j,k] <= x[parent(j), k] + x[j, parent(k)]
        snv_parent_dict = self.createParentDict(self.A)
        cna_parent_dict = self.createParentDict(self.B)
        for j in range(nsnv):
            for k in range(ncna):
                if j != 0 or k != 0:
                    if j in snv_parent_dict.keys() or k in cna_parent_dict.keys():
                        summ = gp.LinExpr()
                        if j in snv_parent_dict.keys():
                            summ += x[snv_parent_dict[j][0], k]
                        if k in cna_parent_dict.keys():
                            summ += x[j, cna_parent_dict[k][0]]
                        model.addConstr(x[j,k] <= summ)

        # 5) Objective Function: Minimize the correction
        obj_sum = gp.LinExpr()
        for i in range(nsamples):
            for j in range(nsnv):
                obj_sum += d_snv[i,j]
            for k in range(ncna):
                obj_sum += d_cna[i,k]
        model.setObjective(obj_sum, gp.GRB.MINIMIZE)

        # 6) Run the Model and evaluate

        # 6.2) Run the model
        model.optimize()

        # 6.3) If infeasible, return
        if model.status != gp.GRB.OPTIMAL:
            return
        
        # 6.4) If feasible, store solutions
        else:
            self.processSolution(model, x, w, z_snv, z_cna)
            self.checkValidity()


    def processSolution(self, model, x, w, z_snv, z_cna):
        self.correction = model.objVal
        sol_x = model.getAttr('x', x)
        sol_z_snv = model.getAttr('x', z_snv)
        sol_z_cna = model.getAttr('x', z_cna)
        raw_snv_edges = [key for key,val in sol_z_snv.items() if val > 0.5]
        raw_cna_edges = [key for key,val in sol_z_cna.items() if val > 0.5]

        clones = [key for key, val in sol_x.items() if val >= 0.5]        
        proportions = model.getAttr('x', w)
        for clone in clones:
            self.G.add_node(clone)
            for i in range(self.n):
                self.G.nodes[clone]['sample_'+str(i)] = proportions[i,clone[0],clone[1]]

        for uv,j in raw_snv_edges:
            u,v = list(self.A.edges)[uv]
            if (u,j) in self.G.nodes and (v,j) in self.G.nodes:
                self.G.add_edge((u,j),(v,j))
        for i,uv in raw_cna_edges:
            u,v = list(self.B.edges)[uv]
            if (i,u) in self.G.nodes and (i,v) in self.G.nodes:
                self.G.add_edge((i,u),(i,v))

        self.numSol = model.SolCount


    # perform checks on tree digraph, to ensure it is valid
    def checkValidity(self):
        nsnv = len(list(self.A.nodes))
        ncna = len(list(self.B.nodes))

        assert self.G.number_of_nodes()==nsnv+ncna-1, 'num nodes not correct'
        assert self.G.number_of_edges()==nsnv+ncna-2, 'num edges not correct'
        assert self.G.has_node((0,0)), 'root node (0,0) does not exist'
        assert len(self.G.in_edges((0,0)))==0, 'root node (0,0) has parent(s)'
        visited = set()
        def dfs(visited, graph, node):
            if node not in visited:
                visited.add(node)
                for neighbour in graph[node]:
                    dfs(visited, graph, neighbour)
        dfs(visited, self.G, (0,0))
        assert len(visited)==len(self.G), 'graph is not connected'








'''
class PCTIsolver:

    def __init__(self, snv_mat, cna_mat, snv_edges, cna_edges):
        
        # 1) store parameters
        self.snv_mat = snv_mat
        self.cna_mat = cna_mat
        self.snv_edges = snv_edges
        self.cna_edges = cna_edges
        self.S = nx.DiGraph()
        self.S.add_edges_from(self.snv_edges)
        self.C = nx.DiGraph()
        self.C.add_edges_from(self.cna_edges)

        # 2) store all final solutions
        self.correction = None
        self.clones = []
        self.proportions = None # dict[(sample, snv_clone, cna_clone)] = proportion
        self.propDf = None
        self.edges = [] # predicted tree: edge list
        self.G = None # predicted tree: networkx digraph
        self.numSol = 0 # number of optimal solutions


    # create dictionaries which store parent of each node
    def createParentDicts(self):
        # snv_dag: the snv tree in DiGraph format
        # cna_dag: the cna tree in DiGraph format
        G = nx.DiGraph()
        G.add_edges_from(self.snv_edges)
        snv_dag = nx.algorithms.transitive_closure_dag(G)
        G.clear()
        G.add_edges_from(self.cna_edges)
        cna_dag = nx.algorithms.transitive_closure_dag(G)

        # parent dictionaries for snv and cna trees
        # dict[child] = [parent A, parent B, ...]
        snv_parent_dict = {}
        for edge in self.snv_edges:
            child = edge[1]
            parent = edge[0]
            if child not in snv_parent_dict.keys():
                snv_parent_dict[child] = [parent]
            else:
                snv_parent_dict[child].append(parent)
        cna_parent_dict = {}
        for edge in self.cna_edges:
            child = edge[1]
            parent = edge[0]
            if child not in cna_parent_dict.keys():
                cna_parent_dict[child] = [parent]
            else:
                cna_parent_dict[child].append(parent)

        return snv_parent_dict, cna_parent_dict



    def solve(self):

        # 1) initialize variables
        nsamples = self.snv_mat.shape[1]
        nsnv = self.snv_mat.shape[0]
        ncna = self.cna_mat.shape[0]

        # 2) create linear programming model and set flags
        model = gp.Model('PCTIsolver')
        model.setParam('OutputFlag',0)
        model.setParam('SolutionLimit',10)

        # 3) initialize linear programming variables
        # x: clones in the final set
        # w, y: the proportion of each clone
        # z: the PREDICTED edges in the cna and snv trees
        # d: correction variables
        x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
        w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
        z_snv = model.addVars(nsnv-1, ncna, vtype = gp.GRB.BINARY, lb = 0, ub = 1, name = 'z_snv')
        z_cna = model.addVars(nsnv, ncna-1, vtype = gp.GRB.BINARY, lb = 0, ub = 1, name = 'z_cna')
        d_snv = model.addVars(nsamples, nsnv, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_snv')
        d_cna = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_cna')

        # 4) Add constraints

        # 4.1) Clone Existance Constraints: encode for the predicted clones
        # encode product w[i,j,k] = y[i,j,k] * x[j,k]
        for i in range(nsamples):
            for j in range(nsnv):
                for k in range(ncna):
                    model.addConstr(w[i,j,k] <= y[i,j,k])
                    model.addConstr(w[i,j,k] <= x[j,k])
                    model.addConstr(w[i,j,k] >= x[j,k] + y[i,j,k] - 1)

        # 4.2) Correction Constraints: difference between input and output matrices
        for i in range(nsamples): # constraint for snv
            for j in range(nsnv):
                summ = gp.LinExpr()
                for k in range(ncna):
                    summ += w[i,j,k]
                model.addConstr(self.snv_mat[j,i] - summ <= d_snv[i,j])
                model.addConstr(summ - self.snv_mat[j,i] <= d_snv[i,j])
        for i in range(nsamples): # constraint for cna
            for k in range(ncna):
                sum = gp.LinExpr()
                for j in range(nsnv):
                    sum += w[i,j,k]
                model.addConstr(self.cna_mat[k,i] - sum <= d_cna[i,k])
                model.addConstr(sum - self.cna_mat[k,i] <= d_cna[i,k])

        # 4.3) Edge Existance Constraints: edge exists <=> endpoints exist
        # encode z_snv[j,k] = x[parent(j), k] * x[j,k]
        # encode z_cna[j,k] = x[j, parent(k)] * x[j,k]
        for edge_idx, edge in enumerate(self.snv_edges):
            parent = edge[0]
            child = edge[1]
            for k in range(ncna):
                model.addConstr(z_snv[edge_idx, k] <= x[parent, k])
                model.addConstr(z_snv[edge_idx, k] <= x[child, k])
                model.addConstr(z_snv[edge_idx, k] >= x[parent, k] + x[child, k] - 1)
        for edge_idx, edge in enumerate(self.cna_edges):
            parent = edge[0]
            child = edge[1]
            for j in range(nsnv):
                model.addConstr(z_cna[j, edge_idx] <= x[j, parent])
                model.addConstr(z_cna[j, edge_idx] <= x[j, child])
                model.addConstr(z_cna[j, edge_idx] >= x[j, parent] + x[j, child] - 1)

        # 4.4) Abundance Constraint: Sum of each proportion matrix row must be 1
        for i in range(nsamples):
            summ = gp.LinExpr()
            for j in range(nsnv):
                for k in range(ncna):
                    summ += w[i,j,k]
            model.addConstr(summ == 1)

        # 4.5) summation of edges constraints: There exists only one edge of each type
        # encode sum_{k} z_snv[j, k] == 1
        # encode sum_{j} z_cna[j, k] == 1
        for uv in range(nsnv-1):
            summ = gp.LinExpr()
            for k in range(ncna):
                summ += z_snv[uv, k]
            model.addConstr(summ == 1)
        for uv in range(ncna - 1):
            summ = gp.LinExpr()
            for j in range(nsnv):
                summ += z_cna[j, uv]
            model.addConstr(summ == 1)
        
        # 4.6) No-Orphans Constraints: orphan clones are not allowed
        # encode x[j,k] <= x[parent(j), k] + x[j, parent(k)]
        snv_parent_dict, cna_parent_dict = self.createParentDicts()
        for j in range(nsnv):
            for k in range(ncna):
                if j != 0 or k != 0:
                    if j in snv_parent_dict.keys() or k in cna_parent_dict.keys():
                        sum = gp.LinExpr()
                        if j in snv_parent_dict.keys():
                            sum += x[snv_parent_dict[j][0], k]
                        if k in cna_parent_dict.keys():
                            sum += x[j, cna_parent_dict[k][0]]
                        model.addConstr(x[j,k] <= sum)

        # 5) Objective Function: Minimize the correction
        obj_sum = gp.LinExpr()
        for i in range(nsamples):
            for j in range(nsnv):
                obj_sum += d_snv[i,j]
            for k in range(ncna):
                obj_sum += d_cna[i,k]
        model.setObjective(obj_sum, gp.GRB.MINIMIZE)

        # 6) Run the Model and evaluate

        # 6.2) Run the model
        model.optimize()

        # 6.3) If infeasible, return
        if model.status != gp.GRB.OPTIMAL:
            return
        
        # 6.4) If feasible, store solutions
        else:
            self.processSolution(model, x, w, z_snv, z_cna)
            self.checkValidity()


    def processSolution(self, model, x, w, z_snv, z_cna):
        self.correction = model.objVal
        self.numSol = model.SolCount
        sol_x = model.getAttr('x', x)
        sol_w = model.getAttr('x', w)
        sol_z_snv = model.getAttr('x', z_snv)
        sol_z_cna = model.getAttr('x', z_cna)
        raw_snv_edges = [key for key,val in sol_z_snv.items() if val > 0.5]
        raw_cna_edges = [key for key,val in sol_z_cna.items() if val > 0.5]
        for uv,j in raw_snv_edges:
            u,v = self.snv_edges[uv]
            self.edges.append(((u,j),(v,j)))
        for i,uv in raw_cna_edges:
            u,v = self.cna_edges[uv]
            self.edges.append(((i,u),(i,v)))
        self.clones = [key for key, val in sol_x.items() if val >= 0.5]
        self.proportions = model.getAttr('x', w)
        self.correction = model.objVal
        G = nx.DiGraph()
        G.add_edges_from(self.edges)
        self.G = G

        nsamples = self.snv_mat.shape[1]
        data = []
        for clone in self.clones:
            snv,cna = clone
            data.append([clone, snv, cna] + [self.proportions[sample, clone[0], clone[1]] for sample in range(nsamples)])
        df = pd.DataFrame(data, columns=['clone', 'snv', 'cna'] + [f'sample_{idx}' for idx in range(nsamples)])
        self.propDf = df


    # perform checks on tree digraph, to ensure it is valid
    def checkValidity(self):
        nsamples = self.snv_mat.shape[1]
        nsnv = self.snv_mat.shape[0]
        ncna = self.cna_mat.shape[0]

        assert self.G.number_of_nodes()==nsnv+ncna-1, 'num nodes not correct'
        assert self.G.number_of_edges()==nsnv+ncna-2, 'num edges not correct'
        assert self.G.has_node((0,0)), 'root node (0,0) does not exist'
        assert len(self.G.in_edges((0,0)))==0, 'root node (0,0) has parent(s)'

        visited = set()
        def dfs(visited, graph, node):
            if node not in visited:
                visited.add(node)
                for neighbour in graph[node]:
                    dfs(visited, graph, neighbour)
        dfs(visited, self.G, (0,0))
        assert len(visited)==len(self.G), 'graph is not connected'
'''

