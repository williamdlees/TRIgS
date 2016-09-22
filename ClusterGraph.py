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


# Render a network from the provided sequences, in which records 1LD away from each other are joined. 

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import math
import sys
import argparse
from Bio import SeqIO
import Levenshtein as ld
import re
import csv


def main(argv):
    parser = argparse.ArgumentParser(description='Create an edge list from a FASTA sequence file.')
    parser.add_argument('infile', help='sequence file (FASTA)')
    parser.add_argument('clstfile', help='cluster file (CD-HIT format)')
    parser.add_argument('colourfile', help='color scheme for samples')
    parser.add_argument('outfile_prefix', help='prefix for output edge and vertex CSV files')
    parser.add_argument('idprefix', help='prefix of generic sequence ids (others will be coloured as highlightcolour)')
    parser.add_argument('highlightcolour', help='highlight colour to use for interesting sequences)')
    parser.add_argument('cutoff', help='distance cutoff for adjacency (expressed as a fraction between 0 and 1)')
    parser.add_argument('-l', '--limit', help='limit on number of nodes to add (for testing)')
    args = parser.parse_args()
    limit = int(args.limit) if args.limit else None
    cutoff = float(args.cutoff)
    
    colours = read_colours(args.colourfile)
    seqs = {}
    seq_counts = {}
    first_sample_seen = {}
    read_seqs(seqs, seq_counts, first_sample_seen, colours, args.infile)
    clusters = read_clusters(args.clstfile, limit, colours, args.idprefix)
    
    edges = []
    verts = []
    cluster_count = 0
    for cluster in clusters:
        cluster_count += add_items(edges, verts, colours, args.idprefix, seqs, seq_counts, first_sample_seen, cutoff, args.highlightcolour, cluster)
    print '%d productive clusters written to file.' % cluster_count
        
    write_results(edges, verts, args.outfile_prefix)
        

# Provide a dict of read counts, indexed by sequence, a dict of seqs, indexed by id, and a dict of the first sample (colour) in which the seq was seen, indexed by sequence
def read_seqs(seqs, seq_counts, first_sample_seen, colours, infile):
    rc = re.compile('.*size\=(\d+)')

    for rec in SeqIO.parse(infile, "fasta"):
        m = rc.match(rec.id)
        if m:
            c = int(m.group(1))
        else:
            c = 1

        seqs[rec.id] = str(rec.seq)        
        seq_counts[str(rec.seq)] = seq_counts.get(str(rec.seq), 0) + c
        
        sample = None
        for i in range(len(colours)):
            if colours[i][0] in rec.id:
                sample = i
                break
        
        if sample is not None:        
            first_sample_seen[str(rec.seq)] = min(first_sample_seen.get(str(rec.seq), 99), sample)
        
    return seq_counts

# Provide a list of cluster members
def read_clusters(clstfile, limit, colours, idprefix):
    clusters = []
    cluster_ids = []
    count = 0
    with open(clstfile, 'r') as fi:
        for line in fi:
            if line[0] == '>':
                if len(cluster_ids) > 1:
                    clusters.append(cluster_ids)
                cluster_ids = []
                count += 1
                if limit and count >= limit:
                    break
            else:
                m = re.search('>(.+?)\.\.\.', line)
                if m:
                    id = m.group(1)
                    # only include 'special' sequences, or sequences which have a defined colour
                    include = False
                    if id[:len(idprefix)] == idprefix:
                        for (sampleid, colour) in colours:
                            if sampleid in id:
                                include = True
                                break
                    else:
                        include = True
                    if include:
                        cluster_ids.append(id)

        if len(cluster_ids) > 1:
            clusters.append(cluster_ids)

    return clusters
    
    
# Read the colour file
def read_colours(colourfile):
    colours = []
    with open(colourfile, 'r') as fi:
        for row in fi:
            cs = row.split(',')
            if len(cs) != 2:
                print 'Error in colour file: %s' % row
                quit()
            colours.append((cs[0], cs[1].rstrip()))

    return colours

    
# Add the edge and vertex lists for a cluster
def add_items(edges, verts, colours, idprefix, seqs, seq_counts, first_sample_seen, cutoff, highlightcolour, cluster):
    # Check the cluster spans at least two samples
    samples_hit = []
    for id in cluster:
        for colour in colours:
            if colour[0] in id:
                if colour[0] not in samples_hit:
                    samples_hit.append(colour[0])
                    if len(samples_hit) > 1:
                        break
                        
    if len(samples_hit) < 2:
        return 0

    newedges = len(edges)
    for i in range(len(cluster)):
        id1 = cluster[i]
        for id2 in cluster[i+1:]:
            if len(seqs[id1]) == len(seqs[id2]):
                limit = int(len(seqs[id1]) * cutoff)
                hd = ld.hamming(seqs[id1], seqs[id2])
                if hd <= limit:
                    edge = {}
                    edge['Source'] = id1
                    edge['Target'] = id2
                    edge['Hamming'] = hd
                    if hd == 1:
                        edge['Color'] = 'black'
                    else:
                        edge['Color'] = 'white'
                    edges.append(edge)
                
    for id in cluster:
        found = False
        for edge in edges[newedges:]:
            if edge['Source'] == id or edge['Target'] == id:
                found = True
                break
        if not found:
            print 'Warning: vertex %s is not connected.' % id
            
        vert = {}
        vert['Id'] = id
        if id[:len(idprefix)] == idprefix:
            vert['color'] = 'black'
            if seqs[id] in first_sample_seen:                
                vert['color'] = colours[first_sample_seen[seqs[id]]][1]
        else:
            vert['color'] = highlightcolour

        vert['size'] = 2 + 2*math.log(seq_counts[seqs[id]])
        verts.append(vert)
    
    return 1
        
# Write result files
def write_results(edges, verts, prefix):
    with open(prefix + '_edges.csv', 'wb') as fo:
        writer = csv.DictWriter(fo, fieldnames=['Source', 'Target', 'Hamming', 'Color'])
        writer.writeheader()
        for edge in edges:
            writer.writerow(edge)

    with open(prefix + '_verts.csv', 'wb') as fo:
        writer = csv.DictWriter(fo, fieldnames=['Id', 'color', 'size'])
        writer.writeheader()
        for vert in verts:
            writer.writerow(vert)
    
if __name__ == "__main__":
    main(sys.argv)

