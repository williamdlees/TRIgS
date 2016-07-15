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


# Remove any FASTA records with duplicated headers

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import os
import sys
import argparse
from Bio import SeqIO

def main(argv):
    parser = argparse.ArgumentParser(description='If records are found with identical IDs, remove all but one on the assumption that they are duplicates.')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default = sys.stdout)
    args = parser.parse_args()

    recs = {}
    for rec in SeqIO.parse(args.infile, "fasta"):
        recs[rec.id] = rec

    SeqIO.write(recs.values(), args.outfile, "fasta")
        
if __name__=="__main__":
    main(sys.argv)

