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


# Given the id of a cluster member, extract records for all members of that cluster from an IMGT-style file

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import math
import sys
import argparse
from Bio import SeqIO
import LevenshteinMaxDist as ld
import re
import igraph
import subprocess
import os
import time
import csv


def main(argv):
    parser = argparse.ArgumentParser(description='Given the ID of a cluster member, extract records for all members of that cluster from an IMGT-style file.')
    parser.add_argument('id', help='ID of the cluster member')
    parser.add_argument('clstfile', help='cluster file (CD-HIT format)')
    parser.add_argument('imgtfile', help='file from which to extract records (IMGT, IgBLASTPlus format)')
    parser.add_argument('outfile', help='output file (IMGT, IgBLASTPlus format)')
    parser.add_argument('-i', '--ignore_size', help='Ignore size designations in IDs (if there is a semicolon in the ID, it and any following text will be ignored)', action="store_true")
    args = parser.parse_args()

    # Find the wanted cluster

    wanted_id = args.id.split(';')[0] if args.ignore_size else args.id

    cluster_ids = {}
    found = False
    with open(args.clstfile, 'r') as fi:
        for line in fi:
            if line[0] == '>':
                if found:
                    break
                cluster_ids = {}
            else:
                m = re.search('>(.+?)\.\.\.', line)
                if m:
                    id = m.group(1)
                    if args.ignore_size:
                        id = id.split(';')[0]
                    cluster_ids[id] = 1
                    if id == wanted_id:
                        found = True
                        
    if not found:
        print 'ID %s was not found in the cluster file.' % args.id
        exit()
        
    print '%d records found in cluster.' % len(cluster_ids)

    fieldnames = None
    wanted_rows = []
    with open(args.imgtfile, 'r') as fi:
        reader = csv.DictReader(fi, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            id = row['Sequence ID']
            if args.ignore_size:
                id = id.split(';')[0]
            if id in cluster_ids:
                wanted_rows.append(row)
                
    with open(args.outfile, 'wb') as fo:
        writer = csv.DictWriter(fo, fieldnames=fieldnames)
        writer.writeheader()
        for row in wanted_rows:
            writer.writerow(row)
            id = row['Sequence ID']
            if args.ignore_size:
                id = id.split(';')[0]
            del(cluster_ids[id])
            
    if len(cluster_ids) > 0:
        print 'Warning: records corresponding to the following IDs were not found:'
        for id in cluster_ids.keys():
            print id


if __name__ == "__main__":
    main(sys.argv)

