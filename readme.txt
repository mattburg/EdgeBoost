Installation:

sudo python setup.py install 

-------------------------------------------------------------------------------------------------------------------

standard usage:

#create edgeBoost instance
edgeBooster = CommunityEdgeBoost(communityFunc,linkAlgorithm,numIterations)

#find communities for graph G
edgeBooster.detect_communities(G)


-------------------------------------------------------------------------------------------------------------------

communityFunc -- community function on igraph object, must take igraph.Graph as input argument and output
                igraph.VertexClustering object, see example below on how to use existing igraph community functions 
 
linkAlgorithm -- link prediction algorithm: "common_neighbors" , "jaccard" or "adamic_adar"

N  = number of boosting iterations, for small to medium networks use 50 or more

-------------------------------------------------------------------------------------------------------------------

example code:

import igraph as ig
from edgeboost.EdgeBoost import CommunityEdgeBoost

G = ig.Graph.Erdos_Renyi(200,0.1)

#creates EdgeBoost object 
edgeBooster = CommunityEdgeBoost(lambda x:x.community_multilevel(),"common_neighbors",numIterations = 10)

#detect communities
communities = edgeBooster.detect_communities(G)

print communities

