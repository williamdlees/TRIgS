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


# Read an IgBlastPlus file and plot germline usage
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
    parser = argparse.ArgumentParser(description='Read an IgBlastPlus file and plot germline usage.')
    parser.add_argument('infiles', help='input files, names separated by commas. (IgBLASTPlus format, nt or aa)')
    parser.add_argument('field', help='input field (e.g. "V-GENE and allele"')
    parser.add_argument('detail', help='F (family), G (germline) or A (allele)')
    parser.add_argument('-t', '--titles', help='titles for each plot, separated by commas')
    parser.add_argument('-l', '--limit', help='limit to at most this many most frequent categories')
    parser.add_argument('-s', '--save', help='Save output to file (as opposed to interactive display)')
    parser.add_argument('-d', '--dupheader', help='Prefix for duplicate count, eg "DUPCOUNT=" for Presto')
    parser.add_argument('-c', '--cols', help='Number of columns for plot')
    parser.add_argument('-a', '--alpha_sort', help='sort columns alphabetically (default is by decreasing size)', action='store_true')    
    parser.add_argument('-f', '--frequency', help='Express chart in terms of frequency rather than number of reads', action='store_true')    
    parser.add_argument('-y', '--ymax', help='Max y-value to use on all charts')

    args = parser.parse_args()
    infiles = args.infiles.split(',')
    if args.detail not in ['F', 'G', 'A']:
        print 'Error: detail must be one of F, G, A'
        exit()
    detail = args.detail
    field = args.field
    limit = int(args.limit) if args.limit else None
    ncols = int(args.cols) if args.cols else 1
    frequency = args.frequency
    ymax = float(args.ymax) if args.ymax else None
    outfile = args.save if args.save else None
    titles = args.titles.split(',') if args.titles else infiles

    nrows = len(infiles) / ncols
    if len(infiles) % ncols != 0:
        nrows += 1
        
    alpha_sort = args.alpha_sort
    dupheader = args.dupheader

    stats = []
    for infile in infiles:
        stats.append(determine_stats(alpha_sort, detail, dupheader, field, frequency, infile, limit))

    if not outfile or len(outfile) < 5 or outfile[-4:] != '.csv':
        plt.figure(figsize=(16,4*nrows))
        plot_number = 1
        for (stat, title) in zip(stats, titles):
            (heights, legends) = stat
            plot_file(heights, legends, frequency, ymax, nrows, ncols, plot_number, title)
            plot_number += 1
        plt.tight_layout()
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
    else:
        with open(outfile, 'wb') as fo:
            writer = csv.writer(fo)
            for (stat, title) in zip(stats, titles):
                (heights, legends) = stat
                writer.writerow([''])
                writer.writerow([title])
                writer.writerow(['Germline'] + legends)
                writer.writerow(['Occurrences'] + heights)


def plot_file(heights, legends, frequency, ymax, nrows, ncols, plot_number, title):
    y_pos = np.arange(len(heights))
    plt.subplot(nrows, ncols, plot_number)
    plt.bar(y_pos, heights,  alpha=0.5)
    plt.xticks(y_pos + 0.5, legends, rotation=-70)
    plt.xlabel(title)
    if ymax:
        plt.ylim(0, ymax)
    if frequency:
        plt.ylabel('Frequency')
    else:
        plt.ylabel('Reads')
    plt.tight_layout()


def determine_stats(alpha_sort, detail, dupheader, field, frequency, infile, limit):
    germ_usage = {}
    with open(infile, 'r') as fi:
        ln = fi.readline()
        sep = ("\t" if "\t" in ln else ",")
        fi.seek(0)
        reader = csv.DictReader(fi, delimiter=sep)
        for row in reader:
            if 'unproductive' not in row['Functionality'] and field in row and row[field] != '':
                germs = to_germ(row[field], detail)
                for germ in germs:
                    germ_usage[germ] = germ_usage.get(germ, 0) + (
                    get_size(row['Sequence ID'], dupheader) if dupheader else 1)
    vals = {}
    total_reads = 0
    for k, v in germ_usage.items():
        total_reads += v
        if v in vals:
            vals[v].append(k)
        else:
            vals[v] = [k]
    heights = []
    legends = []
    indeces = sorted(vals.keys(), reverse=True)
    if limit:
        indeces = indeces[:limit]
    for i in indeces:
        for val in vals[i]:
            heights.append(i)
            legends.append(val)
    if alpha_sort:
        # Juggle the order, now that we know which values we'll be plotting
        germ_usage = {}
        for z in zip(legends, heights):
            germ_usage[z[0]] = z[1]

        heights = []
        legends = []

        for k in sorted(germ_usage.keys()):
            heights.append(germ_usage[k])
            legends.append(k)
    if frequency:
        for i in range(len(heights)):
            heights[i] = float(heights[i]) / total_reads
    return heights, legends


# Convert germline field to a list with the requested amount of detail

def to_germ(germlines, detail):
    result = []
    for germline in germlines.split(','):
        if detail == 'A':
            g = germline
        elif detail == 'G':
            g = germline.split('*')[0]
        else:
            g = germline.split('-')[0]
            
        if g not in result:
            result.append(g)
    
    return result

# Find duplicate size, or return 1

def get_size(s, dupheader):
    count = None
    if dupheader in s:
        spl = s.split(dupheader)
        for i in range(1, len(spl[1])):
            if spl[1][0:i].isdigit():
                count = int(spl[1][0:i])
            else:
                break

    return count if count else 1    

if __name__=="__main__":
    main(sys.argv)

