joinFq
======

Join paired FastQ reads without merging
---------------------------------------

Join forward and reverse reads end-to-end after reverse-complementing
the reverse read.

If a final length is specified, a gap of the appropriate size will be
inserted in the middle.

Run `python joinFq.py -h` for the usage message.

Note that this script does not look for overlaps and will fail if the
desired final length is less than the sum of the input lengths.

Installation
------------

First install python modules via PIP

NumPy
`pip install NumPy`

BioPython
`pip install Biopython`
