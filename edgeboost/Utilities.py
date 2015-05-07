import scipy.sparse as sp
import numpy as np

"""Module for Utility Functions"""


def igraph_to_sparse_matrix(G):
    xs, ys = map(np.array, zip(*G.get_edgelist()))
    if not G.is_directed():
        xs, ys = np.hstack((xs, ys)), np.hstack((ys, xs))
    else:
        xs,ys = xs.T,ys.T
    return sp.coo_matrix((np.ones(xs.shape), (xs, ys)))
