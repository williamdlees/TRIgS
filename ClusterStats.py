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


__author__ = 'William Lees'
__docformat__ = "restructuredtext en"


import sys
import argparse
import re
import numpy as np
import matplotlib.pyplot as plt

# Report statistics of a cluster file produced by CD-HIT or ClusterSeqs.py

def main(argv):
    parser = argparse.ArgumentParser(description="""Report statistics of a cluster file produced by CD-HIT or ClusterSeqs.py.
    """)
    parser.add_argument('clstfile', help='cluster file (CD-HIT format)')
    parser.add_argument('samples', help='comma-separated list of sample names')
    parser.add_argument('-s', '--save', help='Save graphical output to file (as opposed to interactive display)')

    args = parser.parse_args()
    sample_names = args.samples.split(',')
    outfile = args.save if args.save else None

    cluster_count = 0
    cluster_sizes = []
    cluster_spans = []
    unique_clusters = {}
    shared_clusters = {}
    samples_represented = {}
    cluster_sizes_by_sample = {}
    cluster_size_by_sample = {}
    largest_clusters = {}
    for s in sample_names:
        samples_represented[s] = 0
        unique_clusters[s] = 0
        shared_clusters[s] = 0
        cluster_sizes_by_sample[s] = []
        cluster_size_by_sample[s] = 0
        largest_clusters[s] = 0
    cluster_members = 0
    
    with open(args.clstfile, 'r') as fi:
        for line in fi:
            if line[0] == '>':
                if cluster_members > 1:
                    span = 0
                    for k,v in samples_represented.items():
                        if v > 0:
                            span += 1

                    if span > 0:
                        cluster_spans.append(span)
                        cluster_count += 1
                        cluster_sizes.append(cluster_members)
                        if span == 1:
                            for k, v in samples_represented.items():
                                if v > 0:
                                    if v > 1:
                                        unique_clusters[k] += 1
                                    cluster_sizes_by_sample[k].append(cluster_size_by_sample[k])
                                    largest_clusters[k] = max(cluster_size_by_sample[k], largest_clusters[k])
                                    break
                        else:
                            for k, v in samples_represented.items():
                                if v > 0:
                                    if v > 1:
                                        shared_clusters[k] += 1
                                    cluster_sizes_by_sample[k].append(cluster_size_by_sample[k])
                                    largest_clusters[k] = max(cluster_size_by_sample[k], largest_clusters[k])

                cluster_members = 0
                for k, v in samples_represented.items():
                    samples_represented[k] = 0
                    cluster_size_by_sample[k] = 0

            else:
                m = re.search('>(.+?)\.\.\.', line)
                if m:
                    id = m.group(1)
                    sample = find_sample(id, sample_names)
                    if sample:
                        samples_represented[sample] += 1
                        cluster_size_by_sample[sample] += 1
                    cluster_members += 1

    if cluster_members > 1:
        span = 0
        for k, v in samples_represented.items():
            if v > 0:
                span += 1

        if span > 0:
            cluster_spans.append(span)
            cluster_count += 1
            cluster_sizes.append(cluster_members)
            if span == 1:
                for k, v in samples_represented.items():
                    if v > 0:
                        if v > 1:
                            unique_clusters[k] += 1
                        cluster_sizes_by_sample[k].append(cluster_size_by_sample[k])
                        largest_clusters[k] = max(cluster_size_by_sample[k], largest_clusters[k])
                        break
            else:
                for k, v in samples_represented.items():
                    if v > 0:
                        if v > 1:
                            shared_clusters[k] += 1
                        cluster_sizes_by_sample[k].append(cluster_size_by_sample[k])
                        largest_clusters[k] = max(cluster_size_by_sample[k], largest_clusters[k])



    cluster_spans = np.array(cluster_spans)
    
    print 'Number of clusters with more than one member, and at least one member in one of the nominated samples: %d' % cluster_count
    print 'Largest cluster size in that set: %d' % max(cluster_sizes)
    
    print '\nSample\tUnique\tShared\tTotal\tMax Size\tGini Index'
    for s in sample_names:
        print '%s\t%d\t%d\t%d\t%d\t\t%0.2f' % (s, unique_clusters[s], shared_clusters[s], unique_clusters[s] + shared_clusters[s], largest_clusters[s], gini(cluster_sizes_by_sample[s]))        
    
    counts = np.bincount(cluster_spans)
    
    print '\nSpan\tOccurrences'
    for i in range(1, len(counts)):
        print '%d\t%d' % (i, counts[i])
    
    fig, ax = plt.subplots()
    ax.bar(range(1, len(counts)), counts[1:], width=1, align='center')
    ax.set(xticks=range(1, len(counts)), xlim=[0, len(counts)])
    plt.xlabel('Samples spanned')
    plt.ylabel('Occurrences')
    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()
    

def find_sample(id, samples):
    for sample in samples:
        if sample in id:
            return sample
        
# From https://subversion.american.edu/aisaac/notes/incdist.pdf

def gini(x):
    N = len(x)
    x.sort()
    G = sum(x[i] * (N-i) for i in xrange(N))
    G = 2.0 * G/(N * sum(x))
    return (1 + 1./N) - G

    
if __name__ == "__main__":
    main(sys.argv)

