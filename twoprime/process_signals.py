#! /usr/bin/env python

from __future__ import absolute_import, print_function

'''process_signals.py

Process signals from twoprime-seq using different scoring schemes.

Signal must first be loaded in genomedata format
'''

__author__ = 'Jay Hesselberth <jay.hesselberth@gmail.com>'


import sys

from genomedata import Genome
from .signal_scores import scoreA, scoreB, scoreC

__all__ = ['process_signals']


def process_signals(gdarchive, trackname, score_type, verbose = False):
    '''
    Process signals from a genomedata archive.

    Args:
        genome (genomedata.Genome): genomedata archive
        trackname (str): name of track in genomedata archive
        score_type (str): scoring function
        verbose (Optional[bool]): maximum verbosity

    Returns:
        None. Writes signals to stdout or specified file.
    '''

    score_func = _score_func(score_type)

    with Genome(gdarchive) as genome:

        for chrom in genome:
            for contig, continuous in chrom.itercontinuous():
               
                import pdb; pdb.set_trace()

                score = score_func(genome, chrom, pos, trackname, verbose)
                fields = (chrom, pos, pos + 1, score)
                print(map(str, fields), sep='\t')

def _score_func(score_type):
    ''' determine scoring function'''
    
    if score_type == 'A':
        return scoreA
    elif score_type == 'B':
        return scoreB
    elif score_type == 'C':
        return scoreC
    else:
        raise ValueError, "Unknown score type: %s" % score_type


def parse_args(args):

    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('gdarchive', help='genomedata archive')

    parser.add_argument('--trackname', help="track name")

    parser.add_argument('--score-type', choices=['A','B','C'], default='A',
            help="scoring function")

    parser.add_argument('--verbose', action='store_true',
            help="maximum verbosity")

    args = parser.parse_args()

    return args

def main(args=sys.argv[1:]):

    args = parse_args(args)

    kwargs = {'trackname':args.trackname,
              'score_type':args.score_type,
              'verbose':args.verbose}

    return process_signals(args.gdarchive, **kwargs)

if  __name__ == '__main__':
    sys.exit(main())

