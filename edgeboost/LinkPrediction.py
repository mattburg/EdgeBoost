import igraph as ig
import numpy as np
import itertools
import random
import scipy.sparse as sp
from edgeboost import Utilities


"""
This moduale contains all of the link prediction functions: AA, Common Neighbors, Jaccard and Common Neighbor (uniform score). Each scorer function is vectorized code that depends on SciPy sparse arrays, this code WILL NOT scale to large graphs that are dense


version 0.1 Only supports undirected networks, will fix this in future update

"""





def aa_score(G,index_1,index_2):
    score = 0.0
    neighbors_1 = set(G.neighbors(index_1))
    neighbors_2 = set(G.neighbors(index_2))
    common_neighbors = neighbors_1.intersection(neighbors_2)
    
    if len(common_neighbors) == 0:
        return score
    
    n_degrees = G.degree(common_neighbors)
    score = 1./np.log2(n_degrees)
    score = np.sum(score)
    return score


def common_neighbors_score(G,index_1,index_2):
    score = 0.0
    neighbors_1 = set(G.neighbors(index_1))
    neighbors_2 = set(G.neighbors(index_2))

    common_neighbors = neighbors_1.intersection(neighbors_2)
    
    return float(len(common_neighbors))

def jaccard_score(G,index_1,index_2):
    score = 0.0
    neighbors_1 = set(G.neighbors(index_1))
    neighbors_2 = set(G.neighbors(index_2))

    common_neighbors = neighbors_1.intersection(neighbors_2)
    union_neighbors = neighbors_1.union(neighbors_2)
    

    return float(len(common_neighbors))/float(len(union_neighbors))

def score_missing_edges(G,score_func = aa_score):
    nodes = [x.index for x in G.vs]
    #node_pairs = itertools.combinations(nodes,2)
    node_pairs = optimized_node_pairs(G)
    score_dict = {}
    scores = []
    counter = 0
    for node_1,node_2 in node_pairs:
        
        counter+=1
        score = score_func(G,node_1,node_2)
        if score > 0.0 and G.get_eid(node_1,node_2,error=False)==-1 :
            scores.append(((G.vs[node_1]['name'],G.vs[node_2]['name']),score))
    return scores

def optimized_node_pairs(G):
    node_pairs = set()
    for node in G.vs:
        node_pairs.update([(x[0],x[1]) for x in itertools.combinations(G.neighbors(node),2)])
    
    return list(node_pairs)


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
    degrees = A.sum(axis=0)
    print degrees
    exit()
    weights = sp.csr_matrix(np.divide(1,np.log2(degrees)))
    
    A_weighted =  A.multiply(weights)
    A_adamic_adar =  sp.triu(A_weighted*A,k=1)
    
    #read data back into numpy arrays
    x,y = A_adamic_adar.nonzero()
    z  = A_adamic_adar.data
    
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
        k = np.random.randint(1,maxImputed) #sampe num imputed edges
        numsample = min(k+maxImputed,G.ecount())
        linkSample = np.random.choice(linkIndices,numsample,False,scores)
        linkSample = np.array([rowIndices[linkSample],columnIndices[linkSample]]).T
        
        edgeIndices = np.array(G.get_eids(linkSample,error = False))
        linkSample = linkSample[edgeIndices < 0][0:k]
        
        imputationBatches.append(linkSample) #append to batches
    
    G['imputation_batches'] = imputationBatches
     
    return G




def main():
    #devolpment and debug testing
    from sys import getsizeof
    import time
    '''N = 10000
    l = []
    for i in range(N):
        l.append(i)
    exit()
    print G.vcount()
    print G.ecount()
    raw_input()'''
    
    #G = ig.Graph.Erdos_Renyi(10000,0.01)
    G = ig.Graph.Barabasi(100,5)
    print [g.degree() for g in G.vs]
    G.vs['name'] = range(G.vcount())
    print G.ecount()
        
    
    G = ig.Graph()
    G.add_vertices(4)
    G.add_edges([(0,1),(1,2),(2,3),(0,3)])
    G.vs['name'] = ["0","1","2","3"]
    
    #DEBUGGING TEST
    t = time.time()
    baseline_scores = score_missing_edges(G,score_func = aa_score)
    print "time old:   ",time.time()-t
    t = time.time()
    x,y,z = adamic_adar_scorer(G)   
    print "new time:    ",time.time() -t
    new_scores = [((i,j),s) for i,j,s in zip(x,y,z)]
    new_scores = [e for e in new_scores if G.get_eid(e[0][0],e[0][1],error = False) < 0]
    new_scores = [((G.vs['name'][index[0]],G.vs['name'][index[1]]),s) for index,s in new_scores]
    
    new_scores = sorted(new_scores)
    baseline_scores = sorted(baseline_scores)

    for s1,s2 in zip(baseline_scores,new_scores):
        if s1[0] == s2[0] and np.abs(s1[1]-s2[1])<0.00001:
            pass
        else:
            print "fuck"

    exit()
    

if __name__ == "__main__":
    main()



