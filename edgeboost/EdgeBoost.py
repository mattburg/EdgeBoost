from edgeboost import LinkPrediction 
from edgeboost import Community
import igraph as ig
import logging
import sys


"""main module that contains different instatiations of the edgeboost framework, currently only contains the CommunityEdgeBoost object that boosts community detection"""

class CommunityEdgeBoost():
    """ Community EdgeBoost Class, takes community detection algorithm C(G) and
    link prediction function L(G) as input.
    
    paramaters:
        
        communityDetector:   community detection function

        linkPredictor:   link prediction function

        numIterations: number of 'boosting' iterations (default = 10)

        
        """
    
    def __init__(self,communityDetector,linkPredictor,numIterations = 10,threshold = None,connectStrayNodes = True,imputationPercentage = None,verbose = False):
        #assign the community function, must be of form f(G_igraph) --> VertexClustering
        
        if verbose == True:
            logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
        else:
            logging.basicConfig(stream=sys.stdout, level=logging.INFO)


        self.communityDetector = communityDetector
        
        if linkPredictor == "common_neighbors":
            self.linkPredictor = LinkPrediction.common_neighbors_scorer
        elif linkPredictor == "adamic_adar":
            self.linkPredictor = LinkPrediction.adamic_adar_scorer
        elif linkPredictor == "jaccard":
            self.linkPredictor = LinkPrediction.jaccard_scorer
        elif linkPredictor == "local_path":
            self.linkPredictor = lambda x: LinkPrediction.local_path_scorer(x,0.1)
        
        self.numIterations = numIterations
        if threshold != None:
            self.threshold = int(threshold*numIterations)
        else:
            self.threshold = None
        
        self.connectStrayNodes = connectStrayNodes 
        self.linkImputationPercent  = imputationPercentage

    def detect_communities(self,G):
        #create imputated networks
        logging.debug("running link prediction and imputation")
        G = LinkPrediction.link_imputation(G,self.linkPredictor,self.numIterations,self.linkImputationPercent)
        
        logging.debug("clustering imputed networks")
        nodeTable = Community.compute_community_partitions(G,self.communityDetector)
        
        logging.debug("aggregating clusters")
        communities = Community.aggregate_partitions(G,nodeTable,N = self.numIterations,
                tau = self.threshold,connectStrayNodes = self.connectStrayNodes)
            
        return communities
    


class RewireEdgeBoost():
    """ Rewire EdgeBoost Class, takes community detection algorithm C(G) and re-wire probability as input
    
    paramaters:
        
        communityDetector:   community detection function

        linkPredictor:   link prediction function

        numIterations: number of 'boosting' iterations (default = 10)

        
        """
    def __init__(self,communityDetectors,numIterations = 10,rewiringProb = 0.1,threshold = None):
        #assign the community function, must be of form f(G_igraph) --> VertexClustering
        
        self.communityDetectors = communityDetectors
        
        self.rewiringProb = rewiringProb
        
        self.numIterations = numIterations
        if threshold != None:
            self.threshold = int(threshold*numIterations)
        else:
            self.threshold = None
    

    def detect_communities(self,G):
        #create imputated networks
        
        logging.debug("clustering randomly re-wired networks")
        nodeTable = Community.compute_rewire_community_partitions(G,self.communityDetectors,self.numIterations,self.rewiringProb)
        
        logging.debug("aggregating final partition")
        communities = Community.aggregate_partitions(G,nodeTable,N = self.numIterations,
                tau = self.threshold)
            
        return communities



class ThresholdEdgeBoost():
    """ protype threshold booster, cluster network over various threshold settings of a weighted network
    
    paramaters:
        
        communityDetector:   community detection function

        
        """
    def __init__(self,communityDetector,thresholds):
        #assign the community function, must be of form f(G_igraph) --> VertexClustering
        
        self.communityDetector = communityDetector
        
        self.thresholds =  thresholds       
    

    def detect_communities(self,G):
        #create imputated networks
        
        logging.debug("clustering networks")
        nodeTable = Community.compute_threshold_community_partitions(G,self.communityDetector,self.thresholds)
        
        logging.debug("aggregating final partition")
        communities = Community.aggregate_partitions(G,nodeTable,N = len(self.thresholds),
                tau = None)
            
        return communities





