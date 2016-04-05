'''
Combining p-values from dependent tests

A Python port of R code given in below publication

Dai, H., Leeder, J. S., & Cui, Y. (2013). A modified generalized Fisher
method for combining probabilities from dependent tests. Frontiers in
genetics, 5, 32-32.
'''
import numpy as np
from scipy.stats import chi2
from scipy.special import chdtrc as chi2_cdf


class Error(Exception):
    '''
    Base class for exceptions in this module
    '''
    pass


def __handle(vals):
    '''
    Helper function to get rid of boundaries 0 and 1

    Parameters:
        vals: a numpy array of p-values
    '''
    vals[vals >= 0.999999] = 0.999999
    vals[vals < 0.000001] = 0.000001
    return vals


def self_contained(pvals, pmat=None, weights=None):
    '''
    Parameters:
        pvals: a list of p-values to be combined
        pmat: a list of lists of p-values randomly obtained from the data
        weights: a list of weights for p-values in pvals
    '''
    pvals = np.array(pvals, dtype=np.float64)
    if np.isnan(pvals).any():
        raise Error('List pvals contains None value!')
    if pvals.min() < 0:
        raise Error('List pvals contains negative value!')
    if pvals.max() > 1:
        raise Error('List pvals contains > 1 value!')
    pvals = __handle(pvals)

    if pmat is not None:
        pmat = np.matrix(pmat, dtype=np.float64)
        if len(pvals) != len(pmat):
            raise Error('Dimensions of p-values and p-matrix don\'t match!')
        if np.isnan(pmat).any():
            raise Error('List pmat contains None value!')
        if pmat.min() < 0:
            raise Error('List pmat contains negative value!')
        if pmat.max() > 1:
            raise Error('List pmat contains > 1 value!')
        pmat = __handle(pmat)

    if weights is None:
        weights = np.repeat(2, pvals.size)
    else:
        if pvals.size != len(weights):
            raise Error('Dimensions of p-values and weights don\'t match!')
        weights = np.array(weights, dtype=np.float64)
        if np.isnan(weights).any():
            raise Error('List pvals contains None value!')
        if weights.min() < 0:
            raise Error('List weights contains negative value!')

    # e: E(T)
    e = weights.sum()

    # var: Var(T)
    if pmat is None:
        var = 2 * weights.sum()
    else:
        var = np.cov(
                np.apply_along_axis(
                    lambda q: chi2.ppf(1 - q, weights),
                    0,
                    pmat)).sum()

    # v: new degrees of freedom based on Satterthwaite's approximation
    v = 2 * (e**2 / var)
    # T: Lancaster test statistic
    T = chi2.ppf(1 - pvals, weights).sum()
    c = var / 2 / e
    pval = chi2_cdf(v, T / c)

    if pmat is None:
        cor = np.diag(np.repeat(1, pvals.size))
    else:
        cor = np.corrcoef(pmat)

    return {
        'pval': pval,
        'cor': cor
    }


def competitive(pvals, pmat=None, weights=None, n=100000):
    '''
    Parameters:
        pvals: a list of p-values to be combined
        pmat: a list of lists of p-values randomly obtained from the data
        weights: a list of weights for p-values in pvals
        n: number of iterations to compute random pvals when pmat is not given
    '''
    if pmat is None:
        pmat = np.random.randn(len(pvals), n)
    else:
        if len(pvals) != len(pmat):
            raise Error('Dimensions of p-values and p-matrix don\'t match!')
        pmat = np.matrix(pmat, dtype=np.float64)

    random_pvals = []
    for i in xrange(pmat.shape[1]):
        t = self_contained([pmat[0, i], pmat[1, i]], pmat=None, weights=weights)
        random_pvals.append(t['pval'])

    t = self_contained(pvals, pmat=None, weights=weights)

    return {
        'pval': (t['pval'] > random_pvals).mean()
    }
