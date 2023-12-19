from modality.modality_class import Modality
from preprocess.clean_data import readProportionMatrices, getAllPossibleTrees, getTreeEdges
from linear_programs.simultaneous_solver import SIMULTsolver
import itertools
import networkx as nx


def solveSimultaneousPaction(args):
    modalities = initModalities(args)
    props = [modal.U for modal in modalities]
    allTrees = [modal.Gs for modal in modalities]
    allTreeCombinations = itertools.product(*allTrees)
    minCorrection = 100
    minSolutions = []

    for trees in allTreeCombinations:
        solver = SIMULTsolver(props, trees)
        solver.solve()
        if solver.correction != 100:
            if solver.correction < minCorrection:
                minCorrection = solver.correction
                minSolutions = [solver]
            elif solver.correction == minCorrection:
                minSolutions.append(solver)
    
    f = open("/scratch/data/nsdong2/projectPACTION/newpaction/simult_num_solutions.txt", "a")
    f.write(str(len(minSolutions)) + '\n')
    f.close()




def initModalities(args):
    modalities = []
    for i in range(len(args.p)):
        props = readProportionMatrices(args.p[i])
        clones = list(props['genotypes'])
        trees = []
        if args.t[i] == 'None':
            trees = getAllPossibleTrees(props)
        else:
            trees = [getTreeEdges(args.t[i], clones)]
        modalities.append(Modality(props, trees))
    return modalities

