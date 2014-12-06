#!/usr/bin/env python2

'''
Join forward and reverse reads end-to-end
after reverse-complementing the reverse read.

If a final length is specified, a gap of the
appropriate size will be inserted in the middle.

Note that this script does not look for overlaps
and will fail if the desired final length is
less than the sum of the input lengths.

'''

from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
import itertools
import sys


def join_seqs(s1, s2, length=None):
    if length:
        pad_length = length - len(s1) - len(s2)
        try:
            pad = SeqRecord(
                UnknownSeq(pad_length, character='-'),
                letter_annotations = {'phred_quality': [0] * pad_length},
                )
        except ValueError:
            sys.exit('Total length of the two reads exceeds given length (%s)' % (length))
        else:
            s_joined = s1 + pad + s2.reverse_complement()
    else:
        s_joined = s1 + s2.reverse_complement()

    ## assumes the read ID ends in a 2-char suffix for direction (e.g. _1)
    s_joined.id = s1.id[:-2]
    s_joined.description = '' ## not required for fastq
    return s_joined


def write_joined(ffile, rfile, joined_file, length=None):
    freads = SeqIO.parse(ffile, 'fastq')
    rreads = SeqIO.parse(rfile, 'fastq')

    with open(joined_file, 'w') as outfile:
        for fread, rread in itertools.izip(freads, rreads):
            outfile.write(join_seqs(fread, rread, length=length).format('fastq'))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
            'ffile',
            help='Forward reads in FastQ format'
            )
    parser.add_argument(
            'rfile',
            help='Reverse reads in FastQ format'
            )
    parser.add_argument(
            'joined_file',
            help='Name of output file for joined reads'
            )
    parser.add_argument(
            '-l', '--expected-length',
            type=int,
            help='Final length of joined reads if gaps should be added'
            )
    args = parser.parse_args()

    write_joined(
            args.ffile,
            args.rfile,
            args.joined_file,
            args.expected_length
            )
