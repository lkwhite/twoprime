#! /usr/bin/env python

from __future__ import absolute_import, print_function

'''end_coverage.py

Generate coverage of 5p and 3p ends from a BAM file.
'''

__author__ = 'Jay Hesselberth <jay.hesselberth@gmail.com>'


import sys

from pysam import Samfile
from pybedtools import BedTool

__all__ = ['end_coverage']


def end_coverage(bamfilename, strand, end, verbose = False):
    '''
    Generate coverage of 5p and 3p ends in bedGraph format from a
    specified BAM file

    Args:
        bamfilename (str): specific BAM file
        strand (str): pos(`+`) or neg(`-`)
        end (str): `5` or `3`
        verbose (Optional[bool]): maximum verbosity

    Returns:
        None. Writes signals to stdout or specified file.
    '''

    # XXX determine bam outfile name

    with Samfile(bamfilename, 'rb') as samfile,
         Samfile(outfilename, 'wb') as outbamfile:

        for record in samfile:

            # XXX: determine records to process based on end and strand, save
            # to bamfile

            # the major issue is to identify reads that have bona fide 3p
            # end mappings, which have additional, non-genomic sequence at
            # their 3p ends. These have to be distinguished from reads with
            # 3p ends that map precisely to the genome, which presumably
            # represent sequences that were not completely sequenced by
            # the sequencing read length.
            
    # setup bedtool arguments for end and strand
    bedtool_args = {'bg':True, 'strand':strand}

    if end == '5':
        bedtools_args['5':True]
    elif end == '3':
        bedtools_args['3':True]

    coverage_tool = BedTool(outfile).genome_coverage(**bedtool_args)

    # XXX determine bedGraph outfile name

    # XXX write bedGraph data, ideally compressed 

def parse_args(args):

    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('bamfile', help='BAM filename')

    parser.add_argument('--strand', choices=['+', '-'], default=None,
            help="scoring function")

    parser.add_argument('--end', choices=['5', '3'], default=None,
            help="scoring function")

    parser.add_argument('--verbose', action='store_true',
            help="maximum verbosity")

    args = parser.parse_args()

    return args

def main(args=sys.argv[1:]):

    args = parse_args(args)

    kwargs = {'strand':args.strand,
              'end':args.end,
              'verbose':args.verbose}

    return end_coverage(args.bamfile, **kwargs)

if  __name__ == '__main__':
    sys.exit(main())

