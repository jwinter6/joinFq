#!/usr/bin/env python2

import joinFq
from Bio import SeqIO

ffile = 'reads1.fq'
rfile = 'reads2.fq'
joined_file = 'reads_joined.fq'
final_length = 50

# write out
joinFq.write_joined(
        ffile,
        rfile,
        joined_file,
        final_length
        )

# read in
joined = SeqIO.parse(joined_file, 'fastq')

# check
for read in joined:
    assert len(read) == final_length
