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


# Create just enough of a database file to work with Change-o's DefineClones

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
    parser = argparse.ArgumentParser(description='Read an IgBlastPlus file and plot the contained CDR3 lengths.')
    parser.add_argument('-u', '--unique', help='only count unique sequences', action='store_true')
    parser.add_argument('infiles', help='input files, separated by commas (IgBLASTPlus format, nt or aa)')
    parser.add_argument('-t', '--titles', help='titles for each plot, separated by commas')
    parser.add_argument('-c', '--cols', help='Number of columns for plot')
    parser.add_argument('-d', '--dupheader', help='Prefix for duplicate count, eg "DUPCOUNT=" for Presto')
    parser.add_argument('-s', '--save', help='Save output to file (as opposed to interactive display)')
    parser.add_argument('-y', '--ymax', help='Max y-value to use on all charts')
    parser.add_argument('-x', '--xmax', help='Max x-value to use on all charts')
    args = parser.parse_args()
    infiles = args.infiles.split(',')
    titles = args.titles.split(',') if args.titles else None
    ncols = int(args.cols) if args.cols else 1
    xmax = float(args.xmax) if args.xmax else None
    ymax = float(args.ymax) if args.ymax else None
    outfile = args.save if args.save else None
    dupheader = args.dupheader

    if titles and len(titles) != len(infiles):
        print 'Number of titles must match number of files!'
        exit()
    
    nrows = len(infiles) / ncols
    if len(infiles) % ncols != 0:
        nrows += 1

    lengths = []
    for infile in infiles:
        lengths.append(determine_stats(dupheader, infile, args.unique))
        
    if not outfile or len(outfile) < 5 or outfile[-4:] != '.csv':
        plt.figure(figsize=(16,4*nrows))
        plot_number = 1
        for length in lengths:
            title = titles[plot_number-1] if titles else None
            plot_file(length, xmax, ymax, title, nrows, ncols, plot_number)
            plot_number += 1
        plt.tight_layout()
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
    else:
        minlength = 9999
        maxlength = 0
        for length in lengths:
            maxlength = max(maxlength, max([i for i in range(len(length)) if length[i] > 0]))
            minlength = min(minlength, min([i for i in range(len(length)) if length[i] > 0]))

        with open(outfile, 'wb') as fo:
            writer = csv.writer(fo)
            writer.writerow(['Length'] + range(minlength, maxlength))
            ititle = iter(titles) if titles else iter(infiles)
            for length in lengths:
                writer.writerow([next(ititle)] + list(length[minlength:maxlength]))


def plot_file(lengths, xmax, ymax, title, nrows, ncols, plot_number):
    for max in range(599, 0, -1):
        if lengths[max] > 0:
            break
    
    plt.subplot(nrows, ncols, plot_number)
    if ymax:
        plt.ylim(0, ymax)
    if xmax:
        plt.xlim(0, xmax)
    plt.bar(np.arange(max+1), lengths[:max+1])
    if title:
        plt.xlabel(title)
    plt.tight_layout()


def determine_stats(dupheader, infile, unique):
    seen = {}
    with open(infile, 'r') as fi:
        ln = fi.readline()
        sep = ("\t" if "\t" in ln else ",")
        fi.seek(0)
        reader = csv.DictReader(fi, delimiter=sep)
        lengths = np.zeros(600)
        for row in reader:
            if 'unproductive' not in row['Functionality'] and row['CDR3-IMGT'] and (
                not unique or row['CDR3-IMGT'] not in seen):
                lengths[len(row['CDR3-IMGT'])] += (get_size(row['Sequence ID'], dupheader) if dupheader else 1)
                seen[row['CDR3-IMGT']] = 1
    return lengths


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

