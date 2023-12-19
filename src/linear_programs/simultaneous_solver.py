import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools
from modality.modality_class import Modality


class SIMULTsolver:
    def __init__(self, props, trees):
        self.modalities = []
        for i in range(len(props)):
            self.modalities.append(Modality(props[i], [trees[i]]))
            self.modalities[-1].U = self.modalities[-1].U.drop(columns=['genotypes'])
        self.nsamples = self.modalities[0].U.shape[1]

        # 'genotypes' stores the cross product between genotypes of all modalities
        # gives us all possible clones in final clone set
        # ex. [(modalA, modalB,...), (modalA, modalB, ...), ...]
        self.genotypesPool = self.getGenotypesPool()
        self.edgesPool = self.getEdgesPool()


        self.correction = 100
        self.numSol = 0
        self.clones = []
        self.G = nx.DiGraph()
        self.propDf = None

    

    def getNumGenotypes(self):
        ngenotypes = 0
        for modal in self.modalities:
            ngenotypes += modal.U.shape[0]
        return ngenotypes - len(self.modalities) + 1

    def getGenotypesPool(self):
        lens = []
        for modality in self.modalities:
            lens.append([i for i in range(modality.U.shape[0])])
        allPossibleGenotypes = itertools.product(*lens)
        return list(allPossibleGenotypes)


    def solve(self):
        model = gp.Model('PCTIsolver')
        model.setParam('OutputFlag',0)
        model.setParam('SolutionLimit',10)

        # each item is in the format (modalA, modalB,...)
        x = model.addVars(self.genotypesPool, vtype=gp.GRB.BINARY, name='x')

        # each item is in the format (sample#, modalA, modalB,...)
        w = model.addVars(self.nsamples, self.genotypesPool, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        
        # each item is in the format (sample#, modalA, modalB,...)
        y = model.addVars(self.nsamples, self.genotypesPool, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
        
        # each item is in the format (sample#, modal#, genotype#)
        genotypesIterList = self.getGenotypesIterList()
        d = model.addVars(self.nsamples, genotypesIterList, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta')

        # predicted edges - (edgeModal#, edge#, cloneModal#, clone#)
        z = model.addVars(self.edgesPool, vtype = gp.GRB.BINARY, lb = 0, ub = 1, name = 'z')

        # Clone Existance Constraints: encode for the predicted clones
        # encode product w[i,j] = y[i,j] * x[j]
        model = self.addCloneExistanceConstraints(model, w, y, x)

        # Correction Constraints: encode difference between input and output matrices
        model = self.addCorrectionConstraints(model, w, d)

        # Abundance Constraint: enforce that sum of each sample must be 1
        model = self.addAbundanceConstraints(model, w)

        # Summation of edges constraints: enforce that there exists only one edge of each type
        # encode sum_{k} z_snv[j, k] == 1
        # encode sum_{j} z_cna[j, k] == 1  
        model = self.addEdgeSummationConstraints(model, z)

        # Edge Existance Constraints: edge exists <=> endpoints exist
        # encode z_snv[j,k] = x[parent(j), k] * x[j,k]
        # encode z_cna[j,k] = x[j, parent(k)] * x[j,k]
        model = self.addEdgeExistanceConstraints(model, z, x)

        # No-Orphans Constraints: orphan clones are not allowed
        # encode x[j,k] <= x[parent(j), k] + x[j, parent(k)]
        model = self.enforceNoOrphanConstraints(model, x)

        # Objective Function: Minimize the correction
        obj_sum = gp.LinExpr()
        for i in range(self.nsamples):
            for genotype in genotypesIterList:
                tmp = (i,) + genotype
                obj_sum += d[tmp]
        model.setObjective(obj_sum, gp.GRB.MINIMIZE)

        model.optimize()
        if model.status != gp.GRB.OPTIMAL: # if infeasible
            #print('infeasible')
            return
        self.processSolution(model, x, w, z)
        self.checkValidity()



    # TODO: may have issures, double check if needed
    def enforceNoOrphanConstraints(self, model, x):
        parentDicts = self.createParentDicts()
        for genotype in self.genotypesPool:
            if genotype != (0,)*len(self.modalities):
                summ = gp.LinExpr()
                for modalIdx in range(len(self.modalities)):
                    if genotype[modalIdx] != 0:
                        tmp = genotype[:modalIdx] + (parentDicts[(modalIdx, genotype[modalIdx])],) + genotype[modalIdx+1:]
                        summ += x[tmp]
                model.addConstr(x[genotype] <= summ)
        return model


    def createParentDicts(self):
        d = {}
        for modalIdx in range(len(self.modalities)):
            currTree = list(self.modalities[modalIdx].Gs[0].edges)
            for edgeIdx in range(len(currTree)):
                currParent = currTree[edgeIdx][0]
                currChild = currTree[edgeIdx][1]
                d[(modalIdx, currChild)] = currParent
        return d


    def addEdgeExistanceConstraints(self, model, z, x):
        for i in range(len(self.edgesPool)):
            currEdgeRepresent = self.edgesPool[i]
            for j in range(len(self.edgesPool[i])):
                if type(self.edgesPool[i][j]) == tuple:
                    currEdge = self.edgesPool[i][j]
                    currModal = j
                    parentGenotype = self.edgesPool[i][:currModal] + (currEdge[0],) + self.edgesPool[i][currModal+1:]
                    childGenotype = self.edgesPool[i][:currModal] + (currEdge[1],) + self.edgesPool[i][currModal+1:]
                    model.addConstr(z[currEdgeRepresent] <= x[parentGenotype])
                    model.addConstr(z[currEdgeRepresent] <= x[childGenotype])
                    model.addConstr(z[currEdgeRepresent] >= x[parentGenotype] + x[childGenotype] - 1)
        return model


    def addEdgeSummationConstraints(self, model, z):
        d = self.getEdgeIterListDict(self.edgesPool)
        for modalIdx in range(len(self.modalities)):
            currTree = list(self.modalities[modalIdx].Gs[0].edges)
            for edgeIdx in range(len(currTree)):
                summ = gp.LinExpr()
                currEdge = (currTree[edgeIdx][0], currTree[edgeIdx][1])
                for index in d[(modalIdx, currEdge)]:
                    summ += z[self.edgesPool[index]]
                model.addConstr(summ == 1)
        return model


    def getEdgeIterListDict(self, iterList):
        d = {}
        for i in range(len(iterList)):
            for j in range(len(iterList[i])):
                if type(iterList[i][j]) == tuple:
                    if (j,iterList[i][j]) in d:
                        d[(j,iterList[i][j])].append(i)
                    else:
                        d[(j,iterList[i][j])] = [i]
        return d


    # list of potential edges
    # ex. element (1, (0,2), 0) represents the edge (1,0,0) -> (1,2,0)
    # ex. element (0, 0, (0,3)) represents the edge (0,0,0) -> (0,0,3)
    def getEdgesPool(self):
        edgeIterList = set()
        for modalIdx in range(len(self.modalities)):
            currTree = list(self.modalities[modalIdx].Gs[0].edges)
            for edgeIdx in range(len(currTree)):
                for genotype in self.genotypesPool:
                    currTuple = genotype[:modalIdx] + ((currTree[edgeIdx][0], currTree[edgeIdx][1]),) + genotype[modalIdx+1:]
                    edgeIterList.add(currTuple)
        return list(edgeIterList)



    def addCloneExistanceConstraints(self, model, w, y, x):
        for i in range(self.nsamples):
            for j in range(len(self.genotypesPool)):
                    tmp = (i,) + self.genotypesPool[j]
                    model.addConstr(w[tmp] <= y[tmp])
                    model.addConstr(w[tmp] <= x[self.genotypesPool[j]])
                    model.addConstr(w[tmp] >= x[self.genotypesPool[j]] + y[tmp] - 1)
        return model


    def addCorrectionConstraints(self, model, w, d):
        for modalIdx in range(len(self.modalities)): # for each modality
            currModalityProps = self.modalities[modalIdx].U
            for sampleIdx in range(self.nsamples):
                for genotypeIdx in range(len(self.modalities[modalIdx].C)): # for each genotype in current modality
                    summ = gp.LinExpr()
                    currGenotype = self.modalities[modalIdx].C[genotypeIdx]
                    for genotype in self.genotypesPool:  # find every genotype for which the modalIdx-th element is currGenotype, add it
                        if genotype[modalIdx] == currGenotype:
                            tmp = (sampleIdx,) + genotype
                            summ += w[tmp]
                    tmp = (sampleIdx, modalIdx, currGenotype)
                    model.addConstr(currModalityProps.iloc[currGenotype, sampleIdx] - summ <= d[tmp])
                    model.addConstr(summ - currModalityProps.iloc[currGenotype, sampleIdx] <= d[tmp])
        return model

    def addAbundanceConstraints(self, model, w):
        for i in range(self.nsamples):
            summ = gp.LinExpr()
            for genotype in self.genotypesPool:
                tmp = (i,) + genotype
                summ += w[tmp]
            model.addConstr(summ == 1)
        return model


    def getGenotypesIterList(self):
        lens = []
        for i in range(len(self.modalities)):
            modality = self.modalities[i]
            for j in range(modality.U.shape[0]):
                lens.append((i,j))
        return lens
    


    def processSolution(self, model, x, w, z):
        self.correction = model.objVal
        self.numSol = model.SolCount

        sol_x = model.getAttr('x', x)
        self.clones = [key for key, val in sol_x.items() if val >= 0.5]

        sol_z = model.getAttr('x', z)
        raw_edges = [key for key,val in sol_z.items() if val > 0.5]
        for edge in raw_edges:
            for i in range(len(edge)):
                if type(edge[i]) == tuple:
                    u = edge[:i] + (edge[i][0],) + edge[i+1:]
                    v = edge[:i] + (edge[i][1],) + edge[i+1:]
                    self.G.add_edge(u, v)

        props = model.getAttr('x', w)
        data = []
        for clone in self.clones:
            data.append([clone] + [props[(sample,) + clone] for sample in range(self.nsamples)])
        self.propDf = pd.DataFrame(data, columns=['clone'] + [f'sample_{idx}' for idx in range(self.nsamples)])



    # perform checks on tree digraph, to ensure it is valid
    def checkValidity(self):
        ngenotypes = self.getNumGenotypes()

        assert self.G.number_of_nodes()==ngenotypes, 'num nodes not correct'
        assert self.G.number_of_edges()==ngenotypes-1, 'num edges not correct'
        assert self.G.has_node((0,)*len(self.modalities)), 'root node (0,...) does not exist'
        assert len(self.G.in_edges((0,)*len(self.modalities)))==0, 'root node (0,...) has parent(s)'

        visited = set()
        def dfs(visited, graph, node):
            if node not in visited:
                visited.add(node)
                for neighbour in graph[node]:
                    dfs(visited, graph, neighbour)
        dfs(visited, self.G, (0,)*len(self.modalities))
        assert len(visited)==len(self.G), 'graph is not connected'


