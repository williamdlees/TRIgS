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


# Using BLAST, create a CSV file that lists the % identity of the specified sequence to all sequences from the specified germline 

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
from Bio import pairwise2
from Bio.Alphabet import generic_nucleotide
from Bio import SeqIO
from Bio import Phylo

def main(argv):
    parser = argparse.ArgumentParser(description='Create a CSV file that lists the % identity of the specified sequence to all sequences from the specified germline in an IgBLASTPlus format file, and also to their germline.')
    parser.add_argument('queryfile', help='file containing query sequence (IgBLASTPlus format, nt)')
    parser.add_argument('queryseq', help='sequence ID of query sequence')
    parser.add_argument('datafile', help='file containing sequences to match against (IgBLASTPlus format, nt)')
    parser.add_argument('germline', help='sequences from datafile will be used provided they are productive and match this V-germline')
    parser.add_argument('germlinefile', help='File containing the germline sequence (IgBLASTPlus format, nt) (the sequence ID must be identical to the germline')
    parser.add_argument('outfile', help='output file name (CSV)')
    parser.add_argument('-u', '--unproductive', help='include unproductive sequences', action='store_true')
    parser.add_argument('-v', '--verbose', help='verbose output on a restricted set of sequences', action='store_true')
    args = parser.parse_args()

    outfile = args.outfile
    germseq = extract_seq_record(args.germlinefile, args.germline)
    queryseq = extract_seq_record(args.queryfile, args.queryseq)
    
    if '|' in args.germline:
        args.germline = args.germline.split('|')[1] #assume IMGT naming
 
    if not germseq:
        print "Sequence for germline %s not found in germline file %s." % (args.germline, args.germlinefile)
        quit()

    if not queryseq:
        print "Sequence for query %s not found in query file %s." % (args.queryseq, args.queryfile)
        quit()

    count = 0
    with open(args.datafile, 'r') as fi, open(outfile, 'w') as fo:
        writer = csv.writer(fo)
        writer.writerow(['SequenceId', 'TargetDist', 'GermlineDist'])
        reader = csv.DictReader(fi, delimiter='\t')
        for row in reader:
            if match_allele(args.germline, row['V-GENE and allele']) and (args.unproductive or ('unproductive' not in row['Functionality'] and len(row['FR1-IMGT']) > 0 and len(row['FR4-IMGT']))) > 0:
                if args.verbose:
                    print 'Sequence: %s' % row['Sequence ID']
                q_d = distance(row, queryseq, args.verbose, False)
                g_d = distance(row, germseq, args.verbose, True)
                writer.writerow([row['Sequence ID'], round(q_d, 2), round(100.0 - g_d, 2)])
                
            count +=1
            if args.verbose and count > 10:
                break

def distance(query, template, verbose, vregion):
    (total, identical) = (0, 0)
    for field in ('FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT', 'FR4-IMGT'):
        (q, t) = (query[field], template[field])
        if field == 'FR1-IMGT':
            if len(q) < len(t):
                t = t[0-len(q):]
            elif len(q) > len(t):
                q = q[0-len(t):]
            assert(len(q) == len(t))
        elif field == 'FR3-IMGT':
            if vregion:
                if len(q) < len(t):
                    t = t[:len(q)]
                elif len(q) > len(t):
                    q = q[:len(t)]
                assert (len(q) == len(t))
        elif field == 'FR4-IMGT':
            if vregion:
                q = ''
                t = ''
            else:
                if len(q) < len(t):
                    t = t[:len(q)]
                elif len(q) > len(t):
                    q = q[:len(t)]
                assert(len(q) == len(t))

        (tot, i, alignment) = field_distance(q, t)
        total += tot
        identical += i
        
        if verbose:
            print '%s: total: %d identical: %d' % (field, tot, i)
            if alignment:
                print pairwise2.format_alignment(*alignment)
            if field == 'FR1-IMGT' or field == 'FR4-IMGT':
                print q
                print t

    return (100.0 * identical) / total

def field_distance(f1, f2):
    if len(f1) == 0 or len(f2) == 0:
        return (0, 0, None)
    
    alignment = None
    if len(f1) != len(f2):
        for a in pairwise2.align.globalxs(f1, f2, -5, -2):
            (f1, f2, alignment) = (a[0], a[1], a)
            
    (total, identical) = (0, 0)
    for (c1, c2) in zip(f1, f2):
        total += 1
        if c1 == c2:
            identical += 1 

    return (total, identical, alignment)

def extract_seq_record(filename, id):
    with open(filename, 'r') as fi:
        reader = csv.DictReader(fi, delimiter='\t')
        recs = []
        for row in reader:
            if row['Sequence ID'] == id:
                return row

    return None

def strip_subtype(s):
    s = s.split('*')
    return s[0]

def match_allele(a1, a2):
    a1.replace(' or', '')
    a2.replace(' or', '')
    a1 = a1.split(',')
    a2 = a2.split(',')
    for x in a1:
        for y in a2:
            if strip_subtype(x) == strip_subtype(y):
                return True
            
    return False

            

if __name__=="__main__":
    main(sys.argv)

