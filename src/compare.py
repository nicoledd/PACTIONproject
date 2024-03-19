




def runCompare(GTdata, Tdata, ofile):

    numClusters(GTdata.numclusters, Tdata.numclusters, ofile)
    cnaTreeRecall(GTdata.cnatrees, Tdata.cnatrees, ofile)
    cloneRecall(GTdata.clones, Tdata.clones, ofile)
    cnaCloneRecall(GTdata.cnaclones, Tdata.cnaclones, ofile)
    ancestryRecall(GTdata.G, Tdata.G)


def numClusters(GTnumclusters, Tnumclusters, ofile):
    pass 

def cnaTreeRecall(GTtrees, Ttrees, ofile):
    pass 

def cloneRecall(GTclones, Tclones, ofile):
    pass 

def cnaCloneRecall(GTcnaclones, Tcnaclones, ofile):
    pass 

def ancestryRecall(GTdata, Tdata, ofile):
    pass 

