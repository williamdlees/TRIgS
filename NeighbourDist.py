# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# Plot the minimum distance between sequences in the file, uusing Hamming distance between sequences of the same length. 

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import os.path
import sys
import argparse
import csv
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio import SeqIO
from Bio import Phylo
import Levenshtein
import random

def main(argv):
    parser = argparse.ArgumentParser(description='Plot the minimum distance between sequences in the file, uusing Hamming distance between sequences of the same length.')
    parser.add_argument('infile', help='input file (FASTA)')
    parser.add_argument('outprefix', help='prefix for output files')
    parser.add_argument('-l', '--limit', help='limit to at most this many sequences, drawn at random without replacement')    
    parser.add_argument('-v', '--verbose', help='display progress', action='store_true')    
    parser.add_argument('-i', '--interactive', help='display charts interactively', action='store_true')    
    parser.add_argument('-g', '--length_lims', help='axis limits for CDR length plots <xmin,xmax,ymin,ymax>')
    parser.add_argument('-d', '--dist_lims', help='axis limits for CDR length plots <xmin,xmax,ymin,ymax>')
    parser.add_argument('-c', '--csv', help='produce comma separated variable output instead of plots', action='store_true')
    args = parser.parse_args()

    length_lims = None
    dist_lims = None
    
    def split_lims(lims):
        ll = lims.split(',')
        if len(ll) != 4:
            print('Error: limit specifier "%s" should contain four numbers separated by commas.' % lims)
            exit()
        else:
            for i in range(len(ll)):
                ll[i] = float(ll[i])
            return ll
    
    if args.length_lims:
        length_lims = split_lims(args.length_lims)        

    if args.dist_lims:
        dist_lims = split_lims(args.dist_lims)        

    # move seqs into dict, indexed by length. Eliminate duplicates.
    seq_list = []
    seen_seqs = {}
    for seq_record in SeqIO.parse(args.infile, "fasta"):
        seq = str(seq_record.seq)
        if seq not in seen_seqs:
            seq_list.append(seq)
            seen_seqs[seq] = 1
        
    print '%d unique sequences.' % len(seen_seqs)

    if args.limit and len(seq_list) > int(args.limit):
        if args.verbose:
            print 'Limiting to a sample of %d' % int(args.limit)
        seq_list = random.sample(seq_list, int(args.limit))

    seqs = {}
    for seq in seq_list:
        length = len(seq)
        if length in seqs:
            seqs[length].append(seq)
        else:
            seqs[length] = [seq]
            
    maxval = max(seqs.keys()) + 1
    minval = min(seqs.keys())
    sizes = np.zeros(maxval)
    for k, v in seqs.items():
        sizes[k] = len(v)

    sizes = sizes / sizes.sum()

    # plot length distribution

    if not args.csv:                    
        plt.figure()
        plt.bar(range(maxval), sizes, 1/1.5, color='blue')
        plt.xlabel('CDR3 length')
        plt.ylabel('Frequency')
        if length_lims:
            plt.axis(length_lims)
        plt.savefig(args.outprefix + '_length_distribution.pdf')        
        if args.interactive:
            plt.show()
    else:
        with open(args.outprefix + '_length_distribution.csv', 'wb') as fo:
            writer = csv.writer(fo)
            writer.writerow(['Length'] + range(minval, maxval))
            writer.writerow(['Frequency'] + list(sizes[minval:maxval]))
            
    # Calculate min distances
    
    mindists = []
    for seq_length, seq_list in seqs.items():
        if args.verbose:
            print 'Processing sequences of length %d' % seq_length
        for i in range(len(seq_list)-1):
            mindist = 9999
            for j in range(i+1, len(seq_list)):
                h = Levenshtein.hamming(seq_list[i], seq_list[j])
                mindist = min(mindist, h)
            mindists.append(float(mindist)/seq_length)
            
    if not args.csv:                    
        plt.figure()
        plt.hist(mindists, bins=50)
        plt.xlabel('Minimum distance')
        plt.ylabel('Occurrences')
        if dist_lims:
            plt.axis(dist_lims)
        plt.savefig(args.outprefix + '_min_dist.pdf')        
        if args.interactive:
            plt.show()
    else:
        freq, bins = np.histogram(mindists, 50)
        with open(args.outprefix + '_min_dist.csv', 'wb') as fo:
            writer = csv.writer(fo)
            writer.writerow(['Distance'] + list(bins))
            writer.writerow(['Occurrences'] + list(freq))

if __name__=="__main__":
    main(sys.argv)

