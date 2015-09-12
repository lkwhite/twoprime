#! /usr/bin/env python

from __future__ import absolute_import, print_function

'''end_coverage.py

Takes a BAM file on input and writes bedGraph files for `pos` and `neg`
strands containing combined coverage of `5p` and `3p` ends.

Makes sure acquired `3p` ends are acquired within the read, i.e. those
reads are soft-clipped.
'''

__author__ = 'Jay Hesselberth <jay.hesselberth@gmail.com>'
__license__ = 'MIT'

import sys

from gzip import open as gzip_open
from functools import partial
from collections import defaultdict, Counter

from path import path
from pysam import Samfile

__all__ = ['end_coverage']


def end_coverage(bamfilename, verbose = False):
    '''
    Generate coverage of 5p and 3p ends in bedGraph format from a
    specified BAM file. 

    Args:
        bamfilename (str): specific BAM file
        verbose (Optional[bool]): maximum verbosity

    Returns:
        None. Writes signals to stdout or specified file.
    '''

    dictcounter = partial(defaultdict, Counter)
    coverage = defaultdict(dictcounter)  # c['pos']['chr'][pos] = #

    with Samfile(bamfilename, 'rb') as bamfile:

        for record in bamfile:

            # identify reads that have bona fide 3p end mappings, which
            # have additional, non-genomic sequence at their 3p ends.
            #
            # These have to be distinguished from reads with 3p ends that
            # map precisely to the genome, which presumably represent
            # fragments that were not completely covered by the
            # sequencing read.

            chrom = bamfile.getrname(record.rname)

            if not record.is_reverse:
                
                strand = 'pos'

                pos_5p = record.query_alignment_start
                pos_3p = record.query_alignment_end

                coverage[strand][chrom][pos_5p] += 1

                # test whether the aligned 3p end is less than the
                # reference end
                if pos_3p + 1 < record.reference_end:
                    coverage[strand][chrom][pos_3p] += 1
            
            elif record.is_reverse:
                
                strand = 'neg'

                pos_5p = record.query_alignment_end
                pos_3p = record.query_alignment_start

                coverage[strand][chrom][pos_5p] += 1

                # test whether the aligned 3p end is less than the
                # reference end
                if pos_3p + 1 > record.reference_start:
                    coverage[strand][chrom][pos_3p] += 1
     
    # write bedgraph files
    strands = ('pos', 'neg')
    for strand in strands:

        # XXX: make fname for each strand
        fname = 'coverage.{strand}.bg.gz'.format(strand=strand)

        with gzip_open(fname, 'wb') as outfile:
            for chrom, counts in sorted(coverage[strand].items()):
                for pos, count in sorted(counts.items()):

                    fields = (chrom, pos, pos + 1, count)
                    print(*map(str, fields), sep='\t', file=outfile)


def parse_args(args):

    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('bamfile', help='BAM filename')

    parser.add_argument('--verbose', action='store_true',
            help="maximum verbosity")

    args = parser.parse_args()

    return args


def main(args=sys.argv[1:]):

    args = parse_args(args)

    kwargs = {'verbose':args.verbose}

    return end_coverage(args.bamfile, **kwargs)


if  __name__ == '__main__':
    sys.exit(main())

