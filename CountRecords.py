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


# Count FASTA records accounting for duplicate counts in the header

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import os
import sys
import argparse
import csv
import re
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    parser = argparse.ArgumentParser(description='Count FASTA records, accounting for duplicate counts in the header')
    parser.add_argument('dupheader', help='Prefix for duplicate count, eg "DUPCOUNT=" for Presto, "size=" for usearch')
    parser.add_argument('infile', help='input file (FASTA)')
    args = parser.parse_args()

    record_count = 0
    count_with_dupes = 0
    dupheader_not_seen = False
    
    with open(args.infile, 'r') as fi:
        for line in fi:
            if line[0] == '>':
                record_count += 1
                if args.dupheader in line:
                    spl = line.split(args.dupheader)
                    count = None
                    for i in range(1, len(spl[1])):
                        if spl[1][0:i].isdigit():
                            count = int(spl[1][0:i])
                        else:
                            break
                    if count:
                        count_with_dupes += count
                    else:
                        print 'Warning: record count not identified in %s' % line
                else:
                    count_with_dupes +=1
                    dupheader_not_seen = True

    print 'Total FASTA records: %d' % record_count
    print 'Total reads (including duplicates): %d' % count_with_dupes
    if dupheader_not_seen:
        print 'Warning: one or more records did not include a duplicate count (these were assumed to be singletons).'

if __name__=="__main__":
    main(sys.argv)

