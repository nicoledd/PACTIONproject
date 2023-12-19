import networkx as nx

class Modality:
  def __init__(self, props, trees, isCna):
    self.C = list(self.U['genotypes'])

    if isCna == True:
      self.U = props
      self.V = None
      self.Vranges = None
      self.Gs = None
      for i in range(len(self.C)):
        self.C[i] = (int(self.C[i][1]), int(self.C[i][3]))
      self.getCopyStateGraphs()
      self.getVranges()

    else:
      self.U = props['genotypes'] + props.beginswith('sample_') # everything that begins with 'sample'
      self.V = props['genotypes'] + props.beginswith('vaf_') # everything that begins with 'vaf'
      self.Gs = trees
      self.Vranges = None
    

  def getCopyStateGraphs(self):
    self.Gs = []

    assert self.C[0]==(1,1), 'no root (1,1) in copy-number states, or root in incorrect position'

    fullG = nx.DiGraph()
    for i in range(len(self.C)):
      for j in range(i+1, len(self.C)):
        a,b = self.C[i]
        x,y = self.C[j]
        if (a==x and (b==y+1 or b==y-1)) or (b==y and (a==x+1 or a==x-1)):
          fullG.add_edge(self.C[i], self.C[j])
          fullG.add_edge(self.C[j], self.C[i])
    
    self.Gs = [nx.bfs_tree(fullG, source=(1,1))]    

  def getVranges(self):
    self.Vranges = None
    # turn it into a dataframe with a min and max value for each edge in the tree?