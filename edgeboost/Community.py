import igraph as ig
import numpy as np
from collections import deque
import itertools
from copy import deepcopy
import scipy.sparse as sp
from edgeboost import Utilities
from edgeboost import LinkPrediction
import logging


"""This module contains all code for ensemble community detection and aggregation"""


#minimum size of core communities
MIN_CORE_COMM_SIZE = 2.0

def compute_community_partitions(G_imputed,communityFunc,partitionList = None):
    
    #number of times to run community detection
    numIterations = len(G_imputed['imputation_batches'])
    
    #intialize arrays that store node and community assignments
    nodeCommArray = []
    
    #for each pre-computed batch of imputed links, add links to network and cluster
    #storing each clustering in the two tables
    testcommunities = []
    commCounterId = 1
    for i in xrange(numIterations):
        
        imputedIndices = [j+G_imputed.ecount() for j in range(len(G_imputed['imputation_batches'][i]))]
        G_imputed.add_edges(G_imputed['imputation_batches'][i])
        
        #cluster imputed network
        if partitionList != None:
            #this is for debugging only, will be deleted once devolpment on this code is finished
            communityIterator = partitionList[i]
        else:
            communityIterator = communityFunc(G_imputed)
        nodeCommArray.extend([(n,j+commCounterId) for j,comm in enumerate(communityIterator) for n in comm])
        
        #increment counter that assigns ids to each community
        commCounterId+=len(communityIterator)
        G_imputed.delete_edges(imputedIndices)
    
    #create dense array of node --> community array where each row is a node
    #and each column is a vector representing all of the communities for that node
    

    x = np.array([index[0] for index in nodeCommArray])
    y = np.array([index[1] for index in nodeCommArray])
    nodeCommArray = sp.coo_matrix((np.ones(x.shape),(x,y)),dtype = np.uint32).tocsr()
    return nodeCommArray



def compute_rewire_community_partitions(G,communityFunc,numIterations,rewireProb):
    
    #number of times to run community detection
    
    #intialize arrays that store node and community assignments
    nodeCommArray = []
    
    #for each pre-computed batch of imputed links, add links to network and cluster
    #storing each clustering in the two tables
    testcommunities = []
    commCounterId = 1
    for i in xrange(numIterations):
        
        G_rewired = G.copy()
        G_rewired.rewire_edges(rewireProb)
        communityIterator = communityFunc(G_rewired)
        
        nodeCommArray.extend([(n,j+commCounterId) for j,comm in enumerate(communityIterator) for n in comm])
        
        #increment counter that assigns ids to each community
        commCounterId+=len(communityIterator)
    
    #create dense array of node --> community array where each row is a node
    #and each column is a vector representing all of the communities for that node
    

    x = np.array([index[0] for index in nodeCommArray])
    y = np.array([index[1] for index in nodeCommArray])
    nodeCommArray = sp.coo_matrix((np.ones(x.shape),(x,y)),dtype = np.uint8).tocsr()
    return nodeCommArray




def aggregate_partitions(G,nodeCommArray,N,tau = None,connectStrayNodes=True):
    
    #generate core communities
    #N = len(G['imputation_batches'])
   
    
    neighborsMatrix = nodeCommArray*nodeCommArray.T

    #create copy of neighbors matrix for evaluating thresholds
    scoringMatrix = neighborsMatrix.copy()
    neighborsMatrix = sp.triu(neighborsMatrix,format= "csr")
    
    
    if tau == None:
        
        #compute optimal tau
        tau,connectedComponents = compute_optimal_threshold(neighborsMatrix,scoringMatrix,N)
        coreCommunities = [x for x in connectedComponents if len(x) > MIN_CORE_COMM_SIZE]
    
    elif tau != None:
        
        neighborsMatrix.data[neighborsMatrix.data < tau] = 0
        neighborsMatrix.eliminate_zeros()
        connectedComponents = sp.csgraph.connected_components(neighborsMatrix,directed = False)[1]
        connectedComponents = ig.Clustering(connectedComponents)
        coreCommunities = [x for x in connectedComponents if len(x) > MIN_CORE_COMM_SIZE]
    
    #if rare case (usually degenerate) of no core communities, output CC's as final answer
    if len(coreCommunities) < 1:
        return connectedComponents
    
    #output core communities only, this will result in some of the nodes missing from final partition
    if connectStrayNodes == False:
        return coreCommunities

    #merge stray nodes with core communities    
    coreNodes = reduce(lambda x,y:x+y,coreCommunities)
    strayNodes = [v.index for v in G.vs if v.index not in coreNodes]
    finalCommunities = deepcopy(coreCommunities)
    
    #intialize array to store distances
    commDistanceMatrix = np.zeros([len(coreCommunities),G.vcount()])
    
    #compute distances of stray nodes to each core community 
    for commIndex,comm in enumerate(coreCommunities):
        commMatrix = scoringMatrix[comm]
        commMatrix = commMatrix.astype(np.float16)
        commMatrix = commMatrix.mean(axis=0)
        commDistanceMatrix[commIndex,:] = commMatrix
    maxCommIds = np.argmax(commDistanceMatrix,axis=0)
    
    #add stray nodes to the "closest" core community
    for strayNode in strayNodes:
        finalCommunities[maxCommIds[strayNode]].append(strayNode)
    
    
    return finalCommunities





def compute_stability_score(connectedComponents,numComs,neighborsMatrix,N):
    
    numNodes = neighborsMatrix.shape[0]
    stabilityScores = []
    for commId in range(numComs):
        comm = np.where(connectedComponents == commId)[0]
        commSize = float(len(comm))
        
        if commSize <= 1.0:
            intraCommScore = 0
        else:
            intraCommMatrix = neighborsMatrix[comm]
            intraCommMatrix = intraCommMatrix[:,comm]
            intraCommMatrix = sp.triu(intraCommMatrix,k=1)
            intraCommScore = np.divide( intraCommMatrix.sum()/N,(commSize*(commSize-1))/2. )
        
        weightedScore = (commSize/numNodes)*(intraCommScore) 
        stabilityScores.append(weightedScore)
    
    score = np.sum(stabilityScores) 
    return score

def compute_optimal_threshold(neighborsMatrix,scoringMatrix,N):
    numNodes = neighborsMatrix.shape[0]
    thresholds = np.arange(2,N,1)
    optThreshold = -1
    optScore = np.finfo(np.float64).min
    for threshold in thresholds:
        logging.debug("threshold:   {0}".format(threshold))
        neighborsMatrix.data[neighborsMatrix.data < threshold] = 0
        neighborsMatrix.eliminate_zeros()
        numComs,connectedComponents = sp.csgraph.connected_components(neighborsMatrix,directed = False)
        if numComs==1:
            continue
        score = compute_stability_score(connectedComponents,numComs,scoringMatrix,N)
        if score > optScore:
            optScore = score
            optThreshold = threshold
            optConnectedComponents = ig.Clustering(connectedComponents)
    
    try:
        return optThreshold,optConnectedComponents

    except UnboundLocalError:
        optConnectedComponents = ig.Clustering([0 for n in range(neighborsMatrix.shape[0])])
        optThreshold = 2
        return optThreshold,optConnectedComponents


