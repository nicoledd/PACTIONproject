import clonelib


L = [(1, 1), (1, 2), (0,3)]
root_x = 1  
root_y = 1  

# Call the function
result = clonelib.get_cna_trees(set(L), root_x, root_y)

print('L', L)
print('res:')
for res in result:
    print(res)

# print(f"Number of candidate CNA trees: {len(result)}")
# for i in range(len(result)):
#     trees = cnatrees.get_genotype_trees(result[i])
#     print(f"# of genotype trees: {len(trees)}")
#     for j,tree in enumerate(trees):
#         print(f"\nCNA tree: {i} Genotype tree {j}")
#         print(tree)
