from setuptools import setup
from pkg_resources import parse_version


def is_scipy_installed():
    try:
        import scipy

    except ImportError:
        return False
    
    if parse_version(scipy.__version__) >=  parse_version("0.14.0"):
        return True
    else:
        return False


def is_numpy_installed():
    try:
        import numpy

    except ImportError:
        return False

    if parse_version(numpy.__version__) >=  parse_version("1.7.0"):
        return True
    else:
        return False

def is_igraph_installed():
    try:
        import igraph

    except ImportError:
        return False

    if parse_version(igraph.__version__) >=  parse_version("0.6.0"):
        return True
    else:
        return False


if is_igraph_installed() == False:
    print "Install igraph (version >= 0.6.0) before installing EdgeBoost"
    exit()


setup(name='EdgeBoost',
      version='0.1',
      description='EdgeBoost is a Conensus Clustering framework for complex networks',
      author='Matthew Burgess',
      author_email='mattburg@umich.edu',
      url='http://web.eecs.umich.edu/~mattburg/',
      packages=['edgeboost'],
      install_requires = ['numpy >= 1.7','scipy >= 0.14'])


