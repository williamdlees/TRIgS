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

def main(argv):
    parser = argparse.ArgumentParser(description='Copy FASTA records that match the specified pattern to the output file.')
    parser.add_argument('pattern', help='pattern to match (regex')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default = sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default = sys.stdout)
    parser.add_argument('-s', '--seq', help='copy records whose sequences match (default is to match the header)', action='store_true')
    parser.add_argument('-e', '--exc', help='copy records that do not match instead of those that do', action='store_true')
    args = parser.parse_args()
    
    pattern = re.compile(args.pattern.upper())

    recs = []
    for rec in SeqIO.parse(args.infile, "fasta"):
        if args.seq:
            match = re.match(pattern, str(rec.seq).upper())
        else:
            match = re.match(pattern, rec.id)
        if (match and not args.exc) or (args.exc and not match):
            recs.append(rec)
            
    SeqIO.write(recs, args.outfile, "fasta")
            
if __name__=="__main__":
    main(sys.argv)

