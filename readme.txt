hack instillation:


sys.path.append("/edgeboost_path")
import edgeboost.EdgeBoost


usage:

G -- an igraph Graph object
comm_func -- community function on igraph object, use lambda
( e.g  x:x.community_mulilevel() for louvain)
 
link_algorithm -- link prediction algorithm, either "common_neighbors" or
"adamic_adar"

N  = number of boosting iterations, for small networks use 100 or more
 

edgeBooster = CommunityEdgeBoost(comm_func,link_algorithm,numIterations = N)

communities = edgeBooster.detect_communities(G)

(e.g [[0,1,2],[4,5,6]] is a partition of the node vertices into two
partitions [0,1,2] and [4,5,6] on a graph with vertex ids 0-5) 

