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
    parser = argparse.ArgumentParser(description='If identical records are found, remove all but the first.')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default = sys.stdout)
    parser.add_argument('--seq', help='remove records with an identical sequemce (default is to remove records with identical IDs)', action="store_true")

    args = parser.parse_args()

    recs = {}
    i = 0
    for rec in SeqIO.parse(args.infile, "fasta"):
        if args.seq:
            if str(rec.seq) not in recs:
                recs[str(rec.seq)] = rec 
        else:
            if rec.id not in recs:
                recs[rec.id] = rec
        i += 1

    SeqIO.write(recs.values(), args.outfile, "fasta")
    print('%d duplicates removed, leaving %d records.' % (i-len(recs), len(recs)))
        
if __name__=="__main__":
    main(sys.argv)

