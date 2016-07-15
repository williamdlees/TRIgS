# Copyright (c) 2015 William Lees

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# List FASTA records that match the specified pattern

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import os
import sys
import argparse
from Bio import SeqIO
import re
import random


def main(argv):
    parser = argparse.ArgumentParser(
        description='Provide a (randomized) percentage sample of the FASTA file.')
    parser.add_argument('percent', help='percentage required (0-100)')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    perc = float(args.percent)
    if perc < 0 or perc > 100:
        print 'Percentage must be between 0 and 100.'
        quit()
    perc = perc/100

    recs = []
    for rec in SeqIO.parse(args.infile, "fasta"):
        recs.append(rec)

    recs = random.sample(recs, int(perc*len(recs)))
    SeqIO.write(recs, args.outfile, "fasta")


if __name__ == "__main__":
    main(sys.argv)

