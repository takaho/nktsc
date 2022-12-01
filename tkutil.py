import scipy.sparse
import pickle
import gzip
import pandas as pd
import numpy as np
import scipy.io
import os, sys, re
import logging
import argparse
import io

# __loggers = {}

def check_file(filename, date_stat=None):
    """Examine processed file. If the file exists newer than other file, return True
"""
    flag_exist =  os.path.isfile(filename) and os.path.getsize(filename) > 1024
    if not flag_exist:
        return False
    if date_stat is None:
        return True
    if isinstance(date_stat, str) and os.path.isfile(date_stat):
        date_stat = int(os.stat(date_stat).st_mtime)
    mtime = int(os.stat(filename).st_mtime)
    if isinstance(date_stat, int):
        return mtime >= date_stat
    if hasattr(date_stat, 'st_mtime'):
        return mtime >= int(date_stat.st_mtime)
    return True

def get_logger(name=None, stdout=True, logfile=None):
    if name is None:
        name = sys._getframe().f_code.co_name
        pass
    
    logger = logging.getLogger(name)
    # set logging
    for h in logger.handlers:
        h.removeHandler(h)
    def _set_log_handler(logger, handler):#, verbose):
        handler.setFormatter(logging.Formatter('%(asctime)s %(name)s:%(lineno)s %(funcName)s [%(levelname)s]: %(message)s'))
        logger.addHandler(handler)
        return logger
    if logfile is not None:
        _set_log_handler(logger, logging.FileHandler(logfile))
    else:
        stdout = True
    if stdout:
        _set_log_handler(logger, logging.StreamHandler())
    # logger.setLevel(logging.ERROR)
    logger.propagate = False
    return logger

def instanciate_standard_logger(name=None):
    if name is None: name = __name__
    return get_logger(name, stdout=True, logfile='.run.log')

def read_chromosomes(sambamfile, logger=None):
    import collections
    cmd = 'samtools', 'view', '-H', sambamfile
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    chromosomes = collections.OrderedDict()
    with io.TextIOWrapper(proc.stdout, encoding='ascii') as fi:
        for line in fi:
            m = re.match('@SQ\\s+SN:(\\S+)\\s+LN:(\\d+)', line)
            if m:
                chromosomes[m.group(1)] = int(m.group(2))
                if logger:
                    logger.info('{}\t{}bp'.format(m.group(1), m.group(2)))
    proc.stdout.close()
    proc.wait()
    return chromosomes

def phi_correlation(c00, c01=None, c10=None, c11=None):
    if c01 is None:
        c00, c01, c10, c11 = c00[0:4]
    n0_ = c00 + c01
    n1_ = c10 + c11
    n_0 = c00 + c10
    n_1 = c01 + c11
    den = n0_ * n1_ * n_0 * n_1
    if den == 0:
        phi = 0
    else:
        try:
            phi = (c00 * c11 - c01 * c10) / np.sqrt(den)
        except:
            print(c00, c01, c10, c11, n0_, n1_, n_0, n_1)
            raise
    return phi

def get_ordered_leaves(matrix, metrics=None):
    if matrix.shape[0] != matrix.shape[1]:
        return get_ordered_leaves_nondiagonal(matrix, metrics)
    import scipy.cluster.hierarchy
    # print(matrix)
    n = matrix.shape[0]
    X = []
    for i in range(n):
        for j in range(i + 1, n):
            X.append(matrix[i,j])
    Z = scipy.cluster.hierarchy.ward(X)
    ll = scipy.cluster.hierarchy.leaves_list(scipy.cluster.hierarchy.optimal_leaf_ordering(Z, matrix))
    return ll

def get_ordered_leaves_nondiagonal(matrix, metrics=None): # non-diagonal
    import scipy.cluster.hierarchy
    # print(matrix)
    n = matrix.shape[0]
    X = []
    dmat = np.zeros((n, n))
    if metrics is None:
        metrics = lambda x, y: np.dot(x-y, x-y)
    for i in range(n):
        vi = matrix[i,:]
        dmat[i,i] = 0
        for j in range(i + 1, n):
            vj = matrix[j,:]
            # print(vi, vj)
            dmat[i,j] = dmat[j,i] = metrics(vi, vj)
    return get_ordered_leaves(dmat)
    # Z = scipy.cluster.hierarchy.ward(X)
    # ll = scipy.cluster.hierarchy.leaves_list(scipy.cluster.hierarchy.optimal_leaf_ordering(Z, matrix))
    # return ll

