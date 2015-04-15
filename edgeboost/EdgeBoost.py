from edgeboost import LinkPrediction 
from edgeboost import Community
import igraph as ig

"""main module that contains different instatiations of the edgeboost framework, currently only contains the CommunityEdgeBoost object that boosts community detection"""


class CommunityEdgeBoost():
    """ Community EdgeBoost Class, takes community detection algorithm C(G) and
    link prediction function L(G) as input.
    
    paramaters:
        
        communityDetector:   community detection function

        linkPredictor:   link prediction function

        numIterations: number of 'boosting' iterations (default = 10)

        
        """
    def __init__(self,communityDetector,linkPredictor,numIterations = 10,threshold = None,mode = "fast"):
        #assign the community function, must be of form f(G_igraph) --> VertexClustering
        
        self.communityDetector = communityDetector
        
        if linkPredictor == "common_neighbors":
            self.linkPredictor = LinkPrediction.common_neighbors_scorer
        elif linkPredictor == "adamic_adar":
            self.linkPredictor = LinkPrediction.adamic_adar_scorer
        
        self.numIterations = numIterations
        if threshold != None:
            self.threshold = int(threshold*numIterations)
        else:
            self.threshold = None
        
    def detect_communities(self,G):
        #create imputated networks
        print "running link imputation"
        G = LinkPrediction.link_imputation(G,self.linkPredictor,self.numIterations)
        
        print "clustering imputed networks"
        nodeTable = Community.compute_community_partitions(G,self.communityDetector)
        
        print "aggregating final partition"
        communities = Community.get_aggregate_partition(G,nodeTable,
                tau = self.threshold)
            
        return communities
    

def read_network(networkPath,communityPath):
    network_file = open(networkPath)
    G = ig.read(network_file,format= "ncol",directed=False)
    G.simplify() 
    
    community_mapping = [line.rstrip().split("\t") for line in open(communityPath)]
    community_mapping = dict( (x[0],int(x[1])) for x in community_mapping )
    community_array = []
    
    for node in G.vs:
        community_array.append(community_mapping[G.vs[node.index]['name']])

    G.vs['community_assignment'] = community_array
    
    return G


def main():
    import cProfile
    import scipy.sparse as sp
    import sys
    import numpy as np
    import time
    sys.path.append("/Users/matthewburgess/Documents/Research/noisynets/src/network_generation/")
    sys.path.append("/Users/matthewburgess/Documents/Research/noisynets/src/utils/")
    import test_networks
    G = read_network("/Users/matthewburgess/Downloads/binary_networks/smalltest_network.dat",
            "/Users/matthewburgess/Downloads/binary_networks/smalltest_community.dat")

    print G
    cd = CommunityEdgeBoost(lambda x:x.community_multilevel(),"common_neighbors",numIterations = 10)
    cd.detect_communities(G)


if __name__ == "__main__":
    main()
















