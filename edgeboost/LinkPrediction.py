import igraph as ig
import numpy as np
import itertools
import random
import scipy.sparse as sp
from edgeboost import Utilities
import time


"""
This moduale contains all of the link prediction functions: AA, Common Neighbors, Jaccard and Common Neighbor (uniform score). Each scorer function is vectorized code that depends on SciPy sparse arrays, this code WILL NOT scale to large graphs that are dense


version 0.1 Only supports undirected networks, will fix this in future update

"""

def common_neighbors_scorer(G):
    """Computes the CN (Common Neighbor) score for all missing edges in G.

    Returns two lists, one containing missing edges and the other containing the common neighbor scores

    returns list, list """

    #get the sprase adjacency matrix
    A = Utilities.igraph_to_sparse_matrix(G)
    A = A.asformat("csr")

    #common neighbor score matrix
    A_common_neighbors =  sp.triu(A*A,k=1)
    
    #read data back into numpy arrays
    x,y = A_common_neighbors.nonzero()
    z  = A_common_neighbors.data
    
    return x,y,z


def local_path_scorer(G,epsilon = 0.1):
    
    
    #get the sprase adjacency matrix
    A = Utilities.igraph_to_sparse_matrix(G)
    A = A.asformat("csr")

    #common neighbor score matrix
    A_common_neighbors =  A*A
    
    #local path equation
    A_local_path = A_common_neighbors + epsilon*(A_common_neighbors*A)    
    
    A_local_path = sp.triu(A_local_path,k=1)

    #read data back into numpy arrays
    x,y = A_local_path.nonzero()
    z  = A_local_path.data
    
    return x,y,z





def adamic_adar_scorer(G):
    """Computes the CN (Common Neighbor) score for all missing edges in G.

    Returns two lists, one containing missing edges and the other containing the common neighbor scores

    returns list, list """

    #get the sprase adjacency matrix
    A = Utilities.igraph_to_sparse_matrix(G)
    A = A.asformat("csr")

    #common neighbor score matrix
    degrees = A.sum(axis=0)
    weights = sp.csr_matrix(np.divide(1,np.log2(degrees)))
    
    A_weighted =  A.multiply(weights)
    A_adamic_adar =  sp.triu(A_weighted*A,k=1)
    
    #read data back into numpy arrays
    x,y = A_adamic_adar.nonzero()
    z  = A_adamic_adar.data
    
    return x,y,z


def jaccard_scorer(G):
    """Computes the CN (Common Neighbor) score for all missing edges in G.

    Returns two lists, one containing missing edges and the other containing the common neighbor scores

    returns list, list """

    #get the sprase adjacency matrix
    A = Utilities.igraph_to_sparse_matrix(G)
    A = A.asformat("csr")

    #common neighbor score matrix
    degArray = np.asarray(A.sum(axis=0))[0]
    
    #common neighbor score matrix
    A_common_neighbors =  sp.triu(A*A,k=1)
    
    #read data back into numpy arrays
    x,y = A_common_neighbors.nonzero()
    z  = A_common_neighbors.data

    #normalize z according to jaccard equation
    z = z/(degArray[x]+degArray[y]-z)
    
    return x,y,z



def random_edge_scorer(G):
    edgeList = []
    scores = []
    x = [random.randint(0, G.vcount()-1) for i in range(10000)]
    y = [random.randint(0, G.vcount()-1) for i in range(10000)]
    edgeList = [z for z in zip(x,y) if G.get_eid(z[0],z[1],error=False) < 0]
    scores = [1./len(edgeList) for i in range(len(edgeList))]
    return edgeList,scores


def link_imputation(G,linkPredictorFunc = common_neighbors_scorer, numImputations = 10):
    """ randomly samples the missing edges in G according to the provided link prediction function.
    Each iteration randomly samples k edges (which is drawn as a random number between 1 and E). 
    The randomly sampled edges are then added as an attribute to the network G['imputation_batches'], which is a list of lists"""  
    
    #compute missing edges and scores
    rowIndices,columnIndices,scores = linkPredictorFunc(G)
    
    linkIndices = np.arange(len(rowIndices))


    #normalize scores to create probability distribution 
    sumScores = np.sum(scores)
    scores = np.divide(scores,sumScores)

    #set the upper bound for the number of edges to be sampled,
    #needed if G has fewer open triangles then current edges
    maxImputed =  min(len(scores),G.ecount())
    imputationBatches = []
    for i in range(numImputations):

        k = np.random.randint(1,maxImputed) #sample num imputed edges
        numsample = min(k+maxImputed,G.ecount())
        linkSample = np.random.choice(linkIndices,numsample,False,scores)
        linkSample = np.array([rowIndices[linkSample],columnIndices[linkSample]]).T
        edgeIndices = np.array(G.get_eids(linkSample,error = False))
        linkSample = linkSample[edgeIndices < 0][0:k]
        
        imputationBatches.append(linkSample) #append to batches
    
    G['imputation_batches'] = imputationBatches
     
    return G



