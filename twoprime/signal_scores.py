#! /usr/bin/env python

from __future__ import absolute_import, print_function

'''signal_scores.py

scores based on RiboMeth-seq supplementary info described in 

    Profiling of Ribose Methylations in RNA by High-Throughput
    Sequencing (2015) Angewandte Chemie. Birkedal *et al.*

'''

__author__ = 'Jay Hesselberth <jay.hesselberth@gmail.com>'
__license__ = 'MIT'

from numpy import nanmean, nanstd, nansum

# default from Birkedal et al. paper
FLANK_SIZE = 6

# weights used by scoreB and scoreC, listed in order from left to right
WEIGHTS = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.0)

SUM_SCALES = sum(WEIGHTS[abs(i)]
                 for i in range(-FLANK_SIZE, FLANK_SIZE + 1))

__all__ = ['scoreA', 'scoreB', 'scoreC']


def scoreA(gdarchive, chrom, pos, trackname, verbose = False):
    '''
    "Score A is used for detection of ribose methylated positions and is
    based on both the average and standard deviation of the neighboring
    positions. It scales between 0 and 1 with relatively low signal to
    noise ratio. A Matthews's  correlation coefficient[10] is used to
    determine the optimal cut-off value."

    .. math::

        S_i = max \begin{cases}
                  1 - \frac{2n_i + 1}  
                           {\frac{1}{2} |\mu_l - \sigma_l| + n_i 
                            + \frac{1}{2} |\mu_r - \sigma_r| + 1}
                  \\ 
                  0
                  \end{cases}

    Args:
        gdarchive (genomedata.Genome): genomedata archive
        chrom (str): chromosome in genomedata
        pos (int): queried postition
        trackname (str): name of signal track
        verbose (Optional[bool]): maximum verbosity

    Returns:
        float: calculated score or 0.0, whichever is greater
    '''
    n_i = gdarchive[chrom][pos]

    l_data = gdarchive[chrom][pos-FLANK_SIZE:pos, trackname]
    r_data = gdarchive[chrom][pos:pos+FLANK_SIZE, trackname]

    mean_l = nanmean(l_data)
    std_l = nanstd(l_data)
    mean_r = nanmean(r_data)
    std_r = nanstd(r_data)

    score_numer = 2.0 * n_i + 1.0
    score_denom = 0.5 * abs(mean_l - std_l) + \
                  n_i + \
                  0.5 * abs(mean_r - std_r) + \
                  1.0

    score = 1.0 - (score_numer / score_denom)

    return max(score, 0.0)

def _calc_flanks(gdarchive, chrom, pos, trackname):
    '''
    calculate scaled flanks for scoreB and scoreC.

    Args:
        gdarchive (genomedata.Genome): genomedata archive
        chrom (str): chromosome in genomedata
        pos (int): queried postition
        trackname (str): name of signal track

    Returns:
        tuple of floats: scaled left and right flanks
    '''

    l_data = gdarchive[chrom][pos-FLANK_SIZE:pos, trackname]
    r_data = gdarchive[chrom][pos:pos+FLANK_SIZE, trackname]

    scaled_l_flank = nansum(x * y for x, y in zip(l_data, weights))
    scaled_r_flank = nansum(x * y for x, y in zip(l_data, reversed(weights)))

    return (scaled_l_flank, scaled_r_flank)


def scoreB(gdarchive, chrom, pos, trackname, verbose = False):
    '''
    "Score B is based on the weighted average of neighboring positions. It
    has a very high signal to noise ratio and is used for inspection of
    the data in genome browser format. This score is very sensitive to
    small differences at high-level of methylation and is not suitable for
    comparisons of these positions."

    .. math::

        S_i = \frac{|n_i - \frac{1}{2}
                           \left(\frac{\sum_{j=i-\delta}^{i-1} \omega_j n_j}
                                      {\sum_{j=i-\delta}^{i-1} \omega_j} +
                                 \frac{\sum_{j=i+1}^{i+\delta} \omega_j n_j}
                                      {\sum_{j=i+1}^{i+\delta} \omega_j}
                           \right)|
                   }
                   {n_i + 1}

    Args:
        gdarchive (genomedata.Genome): genomedata archive
        chrom (str): chromosome in genomedata
        pos (int): queried postition
        trackname (str): name of signal track
        verbose (Optional[bool]): maximum verbosity

    Returns:
        float: calculated score
    '''

    scaled_l_flank, scaled_r_flank = _calc_flanks(gdarchive, chrom, pos, trackname)

    n_i = gdarchive[chrom][pos, trackname]

    score_numer = abs(n_i - 0.5 * (scaled_l_flank / SUM_SCALES) + \
                                  (scaled_r_flank / SUM_SCALES))
    score_denom = n_i + 1.0
   
    score = score_numer / score_denom

    return score


def scoreC(gdarchive, chrom, pos, trackname, verbose = False):
    '''
    "Score C is a normalized version of score B and expresses the percent
    methylation at a given position using the flanking positions to estimate
    the value corresponding to 100% of molecules methylated at the position.
    In contrast to the other scores, score C is only calculated for the
    selected positions."

    .. math::

        S_i = max \begin{cases}
                  1 - 
                  \frac{n_i}{
                  \frac{1}{2}
                  \left(\frac{\sum_{j=i-\delta}^{i-1} \omega_j n_j}
                             {\sum_{j=i-\delta}^{i-1} \omega_j} +
                        \frac{\sum_{j=i+1}^{i+\delta} \omega_j n_j}
                                      {\sum_{j=i+1}^{i+\delta} \omega_j}
                  \right)}
                  \\ 
                  0
                  \end{cases}

    Args:
        gdarchive (genomedata.Genome): genomedata archive
        chrom (str): chromosome in genomedata
        pos (int): queried postition
        trackname (str): name of signal track
        verbose (Optional[bool]): maximum verbosity

    Returns:
        float: calculated score or 0.0, whichever is greater
    '''

    scaled_l_flank, scaled_r_flank = _calc_flanks(genome, chrom, pos)

    n_i = gdarchive[chrom][pos, trackname]

    score_deom = 0.5 * ((scaled_l_flank / SUM_SCALES) + \
                        (scaled_r_flank / SUM_SCALES))

    score = 1.0 - (n_i / score_denom)

    return max(score, 0.0)


