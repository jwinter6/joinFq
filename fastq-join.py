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
                UnknownSeq(pad_length, character='N'),
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


ffile = 'reads1.fq'
rfile = 'reads2.fq'
joined_file = 'reads_joined.fq'
expected_length = 20


write_joined(ffile, rfile, joined_file, 50)
