# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# Use single-linkage clustering to create clusters, using either Leveshtein distance (supporting clusters of variable length sequences)
# or hanning distance (fixed length only)

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import sys
import argparse
from Bio import SeqIO
#import LevenshteinMaxDist as ld
import Levenshtein as ld
import random
import itertools
import time
from multiprocessing import Pool
import scipy.cluster.hierarchy as sh
import cPickle as pickle

LIMIT = 10000
pool = None
verbose = False
check = False
hamming = False

def main(argv):
    parser = argparse.ArgumentParser(description='Use single-linkage clustering to create clusters of sequences.')
    parser.add_argument('infile', help='input file (FASTA)')
    parser.add_argument('outfile', help='output file (CD-HIT clstr format)')
    parser.add_argument('timelinefile', help='timeline file tracking cluster activity')
    parser.add_argument('cutoff', help='cutoff threshold (between 1 and 0)')
    parser.add_argument('labels', help='comma-separated list of timeline labels, which must be uniquely present in sequence ids')
    parser.add_argument('-l', '--limit', help='limit to at most this many sequences, drawn at random without replacement')    
    #parser.add_argument('-m', '--hamming', help='cluster by hamming distance (default is Levenshtein)', action='store_true')    
    parser.add_argument('-t', '--threads', help='number of threads to use (default is number of processors)')    
    parser.add_argument('-c', '--check', help='check independently that resulting clusters are correctly formed', action='store_true')    
    parser.add_argument('-d', '--dump', help='dump cluster structure to clusters.pickle for debug purposes', action='store_true')
    parser.add_argument('-v', '--verbose', help='display progress', action='store_true')    
    args = parser.parse_args()

    fulltime = time.time()

    infile = args.infile
    cutoff = float(args.cutoff)
    limit = int(args.limit) if args.limit else None
    global verbose
    verbose = args.verbose
    global hamming
    hamming = True
    global check
    check = args.check
    dump = args.dump

    global pool
    if args.threads:
        pool = Pool(int(args.threads))
    else:
        pool = Pool()

    seq_list = read_seqs(infile)
    if limit:
        seq_list = sample_seqs(seq_list, limit)

    cluster_list = cluster_same_length(sort_seqs_by_length(seq_list), cutoff)
    print('Finished single-length clusters after %d minutes' % ((time.time() - fulltime) / 60))
    
    if hamming:
        all_clusters = list(itertools.chain.from_iterable(cluster_list))
    else:
        all_clusters = merge_lengths(cluster_list, cutoff)
        
    if dump:
        dump_structure(all_clusters, 'clusters.pickle')
    print('Finished cluster merges after %d minutes' % ((time.time()-fulltime)/60))

    # deep_merge is not needed as the current merge_lengths function performs an entire merge
    #deep_merge(all_clusters, cutoff)
    write_clusters(all_clusters, args.outfile)
    write_timeline(all_clusters, args.timelinefile, args.labels)

    print('Written cluster files after %d minutes' % ((time.time()-fulltime)/60))

    if check:
        check_clusters(seq_list, all_clusters, cutoff)

    print('Finished after %d minutes' % ((time.time()-fulltime)/60))


def deep_merge(all_clusters, cutoff):
    # Final repetative merge until no more can be achieved
    t1 = time.time()
    merge_new_clusters(all_clusters, [], cutoff)

    merged = True
    merges = 0
    while merged:
        merged = merge_within_clusters(all_clusters, cutoff)
        merges += 1

    t2 = time.time()
    if verbose:
        print('time for deep merge: %0.2f' % (t2 - t1))

def dump_structure(structure, dumpfile):
    print('dumping cluster structure to %s' % dumpfile)
    with open(dumpfile, 'w') as pf:
        pickle.dump(structure, pf)


def merge_lengths(cluster_list, cutoff):
    all_clusters = []
    rate = 0
    for clusters in cluster_list:
        seq_length = len(clusters[0][0][0][0])
        
        if verbose:
            print("Merging %d clusters of length %d into a total of %d" % (len(clusters), seq_length, len(all_clusters)))
        
        all_seqs = sum([len(c[0]) for c in all_clusters])
        new_seqs = sum([len(c[0]) for c in clusters])
        comps = all_seqs * new_seqs
        
        if verbose:
            print("%d comparisons (%d x %d)" % (comps, all_seqs, new_seqs))
            if rate != 0:
                print("predicted time: %d seconds" % (int(comps / rate)))
                
        t1 = time.clock()
        merge_new_clusters(all_clusters, clusters, cutoff)
        t2 = time.clock()
        rate = float(comps) / (t2 - t1) if t1 != t2 else 0
        
        if verbose:
            print("actual time %0.2f seconds (%0.2f comparisons per second)" % (t2 - t1, rate))
    return all_clusters


def cluster_same_length(seqs, cutoff):
    cluster_list = []
    for seq_length, seq_list in seqs.items():
        if verbose:
            print('Processing sequences of length %d (%d items)' % (seq_length, len(seq_list)))
        if len(seq_list) > 1:
            t1 = time.clock()
            chunks = [seq_list[x:x + LIMIT] for x in range(0, len(seq_list), LIMIT)]
            total_clusters = []
            for chunk in chunks:
                clusters = get_clusters(chunk, cutoff)
                merge_new_clusters(total_clusters, clusters, cutoff)
            t2 = time.clock()
            if verbose:
                print('%d clusters processed in %0.2f seconds' % (len(total_clusters), t2 - t1))
            cluster_list.append(total_clusters)
        elif len(seq_list) == 1:
            cluster_list.append([([seq_list[0]], seq_length, seq_length)])
            
    if check:
        seq_nos = 0
        for seq_length, seq_list in seqs.items():
            seq_nos += len(seq_list)
        l = 0
        for c in cluster_list:
            for cs in c:
                l += len(cs[0])
        if seq_nos != l:
            print('Error: sequence losses in cluster_same_length: entry %d exit %d' % (seq_nos, l))

    return cluster_list


def sort_seqs_by_length(seq_list):
    seqs = {}
    for seq, id in seq_list:
        length = len(seq)
        if length in seqs:
            seqs[length].append((seq, id))
        else:
            seqs[length] = [(seq, id)]
    return seqs


def sample_seqs(seq_list, limit):
    if verbose:
        print('Limiting to a sample of %d' % limit)
    seq_list = random.sample(seq_list, limit)
    #seq_list = seq_list[:5000]
    return seq_list


def read_seqs(infile):
    seq_list = []
    for seq_record in SeqIO.parse(infile, "fasta"):
        seq = str(seq_record.seq)
        seq_list.append((seq, seq_record.id))
    
    if verbose:
        print('%d sequences.' % len(seen_seqs))
    return seq_list

def seqs_in_cluster_list(list):
    l = 0
    for x in list:
            l += len(x[0])
    return l

def get_dists(x):
    (lowrow, highrow, seq_list, hamming) = x
    dists = []
    for i in range(lowrow, highrow):
        for j in range(i+1, len(seq_list)):
            if hamming:
                dists.append(ld.hamming(seq_list[i][0], seq_list[j][0]))
            else:
                dists.append(ld.distance(seq_list[i][0], seq_list[j][0]))
    return (dists)

def get_clusters(seq_list, cutoff):
    # Calculate distances as a compressed distance matrix
    seq_length = len(seq_list[0][0])
    cut = int(seq_length*cutoff)
    dists = []
    chunk_size = int(max(1, len(seq_list)/40))

    t1 = time.time()
    dists = pool.map(get_dists, zip(range(0, len(seq_list), chunk_size), range(chunk_size, len(seq_list)+chunk_size-1, chunk_size), itertools.repeat(seq_list), itertools.repeat(hamming)))
    dists = list(itertools.chain.from_iterable(dists))
    t2 = time.time()
    pooltime = t2 - t1

    t1 = time.time()
    Z = sh.linkage(dists, 'single')
    cs = sh.fcluster(Z, cut, criterion='distance')
    t2 = time.time()
    clustertime = t2 - t1

    if verbose:
        print('pool time %0.2f cluster time %0.2f' % (pooltime, clustertime))

    clusters = []
    for i in range(len(set(cs))):
        clusters.append(([], seq_length, seq_length))

    for i in range(len(cs)):
        clusters[cs[i]-1][0].append(seq_list[i])

    if check:
        if seqs_in_cluster_list(clusters) != len(seq_list):
            print('Error: sequence losses in get_clusters: entry %d exit %d' % (len(seq_list), l))

    return clusters

# find clusters to which this new cluster should be merged, and return their indeces in a list
def find_merge_point(x):
    cs1, c, cutoff, hamming = x
    c2, c2_max_length, c2_min_length = c
    merges = []
    for s2, i2 in c2:
        s2_len = len(s2)
        for i in range(len(cs1)):
            if i not in merges:
                (c1, c1_max_length, c1_min_length) = cs1[i]
                for s1, i1 in cs1[i][0]:
                    cut = int(cutoff * min(len(s1), s2_len))
                    if hamming:
                        if c1_min_length != c1_max_length or c2_min_length != c2_max_length:
                            print 'Error: variable length clusters found.'
                        if c1_min_length == c2_min_length and ld.hamming(s1, s2) <= cut:
                            merges.append(i)
                            break
                    else:
                        if c1_min_length - c2_max_length <= cut and c2_min_length - c1_max_length <= cut and ld.distance(s1, s2, cut) <= cut:
                            merges.append(i)
                            break

    return merges

# merge new clusters in the list cs2 into the list cs1. Return an updated list cs1 
def merge_new_clusters(cs1, cs2, cutoff):
    t1 = time.time()
    
    seq_length = 0
    if check:
        seq_length = seqs_in_cluster_list(cs1) + seqs_in_cluster_list(cs2)
    
    merge_points = pool.map(find_merge_point, zip(itertools.repeat(cs1), cs2, itertools.repeat(cutoff), itertools.repeat(hamming)))
    t2 = time.time()
    
    if verbose:
        print('new_cluster time: %0.2f' % (t2-t1))

    merged_clusters = {}   # indexed by merged cluster, specifies where it has been merged to
    
    for i in range(len(cs1)):
        merged_clusters[i] = i
    
    for points, c in zip(merge_points, cs2):
        (c2, c2_max_length, c2_min_length) = c
        if len(points) > 0:
            p = merged_clusters[points[0]]
            (c1, c1_max_length, c1_min_length) = cs1[p]
            for point in points[1:]:
                q = merged_clusters[point]
                if p != q: # don't merge a cluster into itself...
                    (cp, cp_max_length, cp_min_length) = cs1[q]
                    c1.extend(cp)
                    if hamming:
                        if c1_max_length != c2_max_length or c1_min_length != c2_min_length:
                            print 'Error: merging clusters of different lengths.'
                    cs1[p] = (cs1[p][0], max(c1_max_length, cp_max_length), min(c1_min_length, cp_min_length))
                    cs1[q] = None
                    merged_clusters[q] = p
                    for k,v in merged_clusters.items():
                        if v == q:
                            merged_clusters[k] = p

            cs1[p][0].extend(c2)
            cs1[p] = (cs1[p][0], max(c1_max_length, c2_max_length), min(c1_min_length, c2_min_length))
            if hamming:
                if c1_max_length != c2_max_length or c1_min_length != c2_min_length:
                    print 'Error: merging clusters of different lengths.'
        else:
            cs1.append((c2, c2_max_length, c2_min_length))

    for i in range(len(cs1)-1, -1, -1):
        if cs1[i] == None:
            del(cs1[i])

    if check:
        if seq_length != seqs_in_cluster_list(cs1):
            print('Error: sequence losses in merge_new_clusters: entry %d exit %d' % (seq_length, seqs_in_cluster_list(cs1)))

    return cs1

def merge_within_clusters(cs1, cutoff):
    merged = False
    merge_points = []

    seq_length = 0
    if check:
        seq_length = seqs_in_cluster_list(cs1)

    def upto(data, start=0):
        n = start
        while True:
            yield data[:n]
            n += 1
    
    if verbose:
        print('clusters to merge: %d' % (len(cs1)))
    
    t1 = time.time()
    foo = zip(upto(cs1, 1), cs1[1:], itertools.repeat(cutoff))
    merge_points = pool.map(find_merge_point, zip(upto(cs1, 1), cs1[1:], itertools.repeat(cutoff), itertools.repeat(hamming)))
    # need to adjust merge_points indexing because we started looking from cs1[1]...
    merge_points = [None] + merge_points
    t2 = time.time()
    
    if verbose:
        print('within_cluster time: %0.2f' % (t2-t1))

    remove_items = []    
    for i in range(len(merge_points)-1, -1, -1):
        if merge_points[i] is not None:
            (c1, c1_max_length, c1_min_length) = cs1[merge_points[i]]
            (c2, c2_max_length, c2_min_length) = cs1[i]
            cs1[merge_points[i]][0].extend(c2)
            cs1[merge_points[i]] = (cs1[merge_points[i]][0], max(c1_max_length, c2_max_length), min(c1_min_length, c2_min_length))
            remove_items.append(i)
            merged = True
            
    if len(remove_items) > 0:            
#        remove_items.reverse()
        for i in remove_items:
            del cs1[i]

    if check:
        if seq_length != seqs_in_cluster_list(cs1):
            print('Error: sequence losses in merge_within_clusters: entry %d exit %d' % (seq_length, seqs_in_cluster_list(cs1)))

    return (merged)

def write_clusters(cluster_list, outfile):
    with open(outfile, 'w') as fo:
        clusternum = 0
        for cluster, max_len, min_len in cluster_list:
            fo.write('>Cluster %d\n' % clusternum)
            clusternum += 1
            seqnum = 0
            for seq, id in cluster:
                fo.write('%d %daa, >%s... \n' % (seqnum, len(seq), id))
                seqnum += 1

def write_timeline(cluster_list, timelinefile, labels):
    with(open(timelinefile, 'w')) as fo:
        fo.write('#samples %s\n' % labels)
        labels = labels.split(',')
        times = [str(x) for x in range(len(labels))]        
        fo.write('#times %s\n' % ' '.join(times))
        sums = {}
        index = 1
        for cluster, max_len, min_len in cluster_list:
            row = {}
            for seq,id in cluster:
                for label in labels:
                    if label in id:
                        row[label] = 1 + row.get(label, 0)
                        sums[label] = 1 + sums.get(label, 0)
            if len(row) > 1:
                counts = []
                for label in labels:
                    if label in row:
                        counts.append(str(row[label]))
                    else:
                        counts.append('0')
                fo.write('%s_%d %s\n' % ('C', index, ' '.join(counts)))
                index += 1
        counts = []
        for label in labels:
            counts.append(str(sums[label]) if label in sums else '0')
        fo.write('#sums %s\n' % ' '.join(counts))
        
def check_clusters(seq_list, all_clusters, cutoff):
    # Sequences in seq_list are unique
    print('Checking for unique input sequences')
    seqs = {}
    for seq in seq_list:
        seqs[seq[0]] = seqs.get(seq, 0) + 1
        
    for k, v in seqs.items():
        if v != 1:
            print('Error: sequence %s appears in seq_list more than once.' % k)
            
    # There is a one-to-one correspondence between sequences in seq_list and sequences in all_clusters
    print('Checking for one-to-one correspondence between input and output sequences.')
    for cluster, max_len, min_len in all_clusters:
        for seq, id in cluster:
            if seq not in seqs:
                print('Error: sequence %s is in all_clusters but not in seq_list.' % seq)
            else:
                seqs[seq] += 1

    for k, v in seqs.items():
        if v < 2:
            print('Error: sequence %s appears in seq_list but not in all_clusters.' % k)
        elif v > 2:
            print('Error: sequence %s appears in all_clusters more than once.' % k)
            
    # Each sequence in a cluster has a nearest neighbour within the cutoff distance
    print('Checking cluster membership.')
    for cluster, max_len, min_len in all_clusters:
        if len(cluster) > 1:
            for i in range(len(cluster)):
                seq1, id1 = cluster[i]
                belongs = False
                for j in range(len(cluster)):
                    seq2, id2 = cluster[j]
                    cut = int(cutoff * min(len(seq1), len(seq2)))
                    if hamming:
                        if i != j and ld.hamming(seq1, seq2) <= cut:
                            belongs = True
                            break
                    else:
                        if i != j and ld.distance(seq1, seq2, cut) <= cut:
                            belongs = True
                            break
                if not belongs:
                    print('Error: sequence %s (id %s) has no sufficiently near neighbours in its cluster.' % (cluster[i]))
                
    # No clusters are mergeable
    print('Checking that clusters are distinct.')
    for i in range(len(all_clusters)):
        c1 = all_clusters[i][0]
        for c2, max_len, min_len in all_clusters[i+1:]:
            for s1, i1 in c1:
                for s2, i2 in c2:
                    cut = int(cutoff * min(len(s1), len(s2)))
                    if hamming:
                        if len(s1) == len(s2) and ld.hamming(s1, s2) <= cut:
                            print('Error: sequences %s (id %s) and %s (id %s) are in different clusters but are within the cutoff distance.' % (s1, i1, s2, i2))
                    else:
                        if ld.distance(s1, s2, cut) <= cut:
                            print('Error: sequences %s (id %s) and %s (id %s) are in different clusters but are within the cutoff distance.' % (s1, i1, s2, i2))

if __name__=="__main__":
    main(sys.argv)

