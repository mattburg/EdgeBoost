Installation:

sudo python setup.py install 


usage:

G -- an igraph Graph object

comm_func -- community function on igraph object, use lambda
( e.g  x:x.community_multilevel() for louvain)
 
link_algorithm -- link prediction algorithm, either "common_neighbors" , "jaccard" or
"adamic_adar"

N  = number of boosting iterations, for small networks use 100 or more


example:

import igraph as ig
from edgeboost.EdgeBoost import CommunityEdgeBoost

G = ig.Graph.Erdos_Renyi(200,0.1)
 
edgeBooster = CommunityEdgeBoost(lambda x:x.community_multilevel(),"common_neighbors",numIterations = 10)

communities = edgeBooster.detect_communities(G)

print communities

