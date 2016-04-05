# Combining p-values from dependent tests

A Python port of R code given in below publication

Dai, H., Leeder, J. S., & Cui, Y. (2013). A modified generalized Fisher
method for combining probabilities from dependent tests. Frontiers in
genetics, 5, 32-32.

## Installation

    pip install git+https://github.com/gungorbudak/cdpvals.git

## Usage

    from cdpvals import self_contained
    from cdpvals import competitive


    pvals = [0.06, 0.15]
    pmat = [ [0.02,0.06,0.07,0.01,0.02,0.09,0.01], [0.01,0.10,0.12,0.14,0.07,0.09,0.10] ]
    print self_contained(pvals, pmat)
    print competitive(pvals, pmat)

## Documentation

This package is using Lancaster procedure (a generalized Fisher's method) with weight functions including Satterthwaite's approximation to model correlations among p-values.

### self_contained

Description

#### Parameters

pvals: a list of p-values to be combined
pmat: a list of lists of p-values randomly obtained from the data
weights: a list of weights for p-values in pvals

### competitive

Description

#### Parameters

pvals: a list of p-values to be combined
pmat: a list of lists of p-values randomly obtained from the data
weights: a list of weights for p-values in pvals
n: number of iterations to compute random pvals when pmat is not given
