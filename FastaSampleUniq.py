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


# Find the distribution of the number of unique sequences in a sample taken from the input file, accounting for duplicate counts in the header 

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import os
import sys
import argparse
import itertools
import random
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


def main(argv):
    parser = argparse.ArgumentParser(
        description='Find the average number of unique sequences in a set of random samples taken from the input file, accounting for duplicate counts in the header')
    parser.add_argument('sample_size', help='Number of records to sample')
    parser.add_argument('iterations', help='Number of times to sample')
    parser.add_argument('dupheader', help='Prefix for duplicate count, eg "DUPCOUNT=" for Presto')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input file (FASTA)')
    parser.add_argument('--plot', help='plot the distribution of the number in each sample', action="store_true")
    args = parser.parse_args()

    recs = []
    seqs = {}
    fasta_recs = 0

    for rec in SeqIO.parse(args.infile, "fasta"):
        fasta_recs += 1
        count = None
        if args.dupheader in rec.id:
            spl = rec.id.split(args.dupheader)
            for i in range(1, len(spl[1])):
                if spl[1][0:i].isdigit():
                    count = int(spl[1][0:i])
                else:
                    break
        if count is None:
            count = 1
        for _ in range(count):
            recs.append(str(rec.seq))
        seqs[str(rec.seq)] = 1
            
    print 'Input file contains %d FASTA records representing %d reads and %d unique sequences.' % (fasta_recs, len(recs), len(seqs))
    
    sample_size = int(args.sample_size)
    if sample_size > len(recs):
        print 'Error: requested sample size is greater than the total number of reads. '
        quit()
    
    results = []    
    for _ in range(int(args.iterations)):
        #print '.'
        uniques = {}
        sample = random.sample(recs, sample_size)
        for rec in sample:
            uniques[rec] = 1
        results.append(len(uniques))

    res = np.array(results)
    print 'mean number of unique values in sample: %0.2f' % np.mean(res)
    print 'standard deviation: %0.2f' % np.std(res)
    print '95 precent confidence limits: +/- %0.2f' % (1.96*np.std(res))

    if args.plot:
        plt.hist(results)
        plt.show()

if __name__ == "__main__":
    main(sys.argv)

