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
import matplotlib.colors as mcolors
import matplotlib.mlab as mlab
import itertools

def main(argv):
    parser = argparse.ArgumentParser(description='Read an IgBlastPlus file and plot the contained CDR3 lengths.')
    parser.add_argument('infiles', help='input files, separated by commas (IgBLASTPlus format, nt or aa)')
    parser.add_argument('-b', '--barcolour', help='colour or list of colours for bars')
    parser.add_argument('-c', '--cols', help='Number of columns for plot')
    parser.add_argument('-d', '--dupheader', help='Prefix for duplicate count, eg "DUPCOUNT=" for Presto')
    parser.add_argument('-g', '--gradientfill', help='fill bars with a gradiented colour', action='store_true')    
    parser.add_argument('-gh', '--grid_horizontal', help='horizontal grid lines', action='store_true')
    parser.add_argument('-gv', '--grid_vertical', help='vertical grid lines every n bars')
    parser.add_argument('-ga', '--gauss', help='plot best-fit Gaussian distribution', action='store_true')    
    parser.add_argument('-s', '--save', help='Save output to file (as opposed to interactive display)')
    parser.add_argument('-sz', '--size', help='Figure size (x,y)')
    parser.add_argument('-t', '--titles', help='titles for each plot, separated by commas')
    parser.add_argument('-u', '--unique', help='only count unique sequences', action='store_true')
    parser.add_argument('-w', '--width', help='relative bar width (number between 0 and 1)')
    parser.add_argument('-xmax', '--xmax', help='Max x-value to use on all charts')
    parser.add_argument('-xmin', '--xmin', help='Min x-value to use on all charts')
    parser.add_argument('-y', '--ymax', help='Max y-value to use on all charts')
    args = parser.parse_args()
    infiles = args.infiles.split(',')
    titles = args.titles.split(',') if args.titles else [""]
    ncols = int(args.cols) if args.cols else 1
    xmin = int(args.xmin) if args.xmin else 1
    xmax = int(args.xmax) if args.xmax else None
    ymax = float(args.ymax) if args.ymax else None
    outfile = args.save if args.save else None
    dupheader = args.dupheader
    bar_width = float(args.width) if args.width else 1.0
    mapcolour = args.barcolour if args.barcolour else 'blue'
    mapcolour = mapcolour.split(',')
    grid_vertical = int(args.grid_vertical) if args.grid_vertical else False
    gauss = args.gauss
    (sizex, sizey) = args.size.split(',') if args.size else (8,4)

    nrows = len(infiles) / ncols
    if len(infiles) % ncols != 0:
        nrows += 1

    lengths = []
    for infile in infiles:
        lengths.append(determine_stats(dupheader, infile, args.unique))
        
    if not outfile or len(outfile) < 5 or outfile[-4:] != '.csv':
        plt.figure(figsize=(float(sizex),float(sizey)))
        plot_number = 1
        for (length, title, colour) in zip(lengths, itertools.cycle(titles), itertools.cycle(mapcolour)):
            plot_file(length, xmin, xmax, ymax, title, nrows, ncols, plot_number, colour, bar_width, args.gradientfill, args.grid_horizontal, grid_vertical, gauss)
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


def plot_file(lengths, xmin, xmax, ymax, title, nrows, ncols, plot_number, mapcolour, bar_width, gradientfill, grid_horizontal, grid_vertical, gauss):
    if not xmax:
        for xmax in range(len(lengths)-1, 0, -1):
            if lengths[xmax] > 0:
                break
    
    ax = plt.subplot(nrows, ncols, plot_number)
    ax.tick_params(direction='out', top=False, right=False)

    if title != "":
        plt.xlabel(title)

    bar_pos = np.arange(xmax + 1)
    
    # if xmax is set, all the x axes will be the same, so only put labels on the bottom row
    # if we have a vertical grid, label just one point in each column
    
    if xmax is None or plot_number > (nrows-1)*ncols:
        if grid_vertical:
            labels = []
            if grid_vertical % 2 == 0:
                midpoint = grid_vertical/2
            else:
                midpoint = int(grid_vertical/2)
            for tick in bar_pos:
                if tick % grid_vertical == midpoint:
                    labels.append(str(tick))
                else:
                    labels.append(' ')
        else:
            labels = ['%d' % x for x in bar_pos]
    else:
        labels = []
            
    plt.xticks(bar_pos+0.5, labels)

    if ymax:
        plt.ylim(0, ymax)
    if xmax or xmin != 0:
        plt.xlim(xmin, xmax)

    if bar_width < 1.:
        bar_pos = bar_pos + (1 - bar_width) / 2.

    if grid_horizontal:
        plt.grid(which='major', axis='y', c='black', linestyle='-', alpha=0.6, zorder=1)

    if gradientfill:
        gbar(bar_pos, lengths, mapcolour, width=bar_width)
    else:
        plt.bar(bar_pos, lengths, width=bar_width, color=mapcolour, zorder=10)

    if grid_vertical:
        pos = 0
        while pos < xmax:
            ymin, ymax = plt.ylim()
            plt.plot([bar_pos[pos] - (1 - bar_width)/2, bar_pos[pos] - (1 - bar_width)/2], [ymin, ymax], c='black', linestyle='-', alpha=0.6, zorder=1)
            pos += grid_vertical
            
    # Remove every other y label because we get far too many by default
    
    locs, labels = plt.yticks()
    newlocs = []
    newlabels = []
    
    for i in range(0, len(labels)):
        if i % 2 != 0:
            newlocs.append(locs[i])
            newlabels.append(str(int(locs[i])))
            
    plt.yticks(newlocs, newlabels)
            
    if gauss:
        values = []
        for i in range(len(lengths)-1):
            if lengths[i] > 0:
                foo = [i]*lengths[i]
                values += ([i] * lengths[i])
        values = np.array(values)
        mean = np.mean(values)
        variance = np.var(values)
        sigma = np.sqrt(variance)
        x = np.linspace(min(values), max(values), 100)
        plt.plot(x, mlab.normpdf(x, mean, sigma)*len(values), zorder=20)

    ax.set_aspect('auto')
    plt.tight_layout()


def gbar(x, y, mapcolour, width=1, bottom=0):
    X = [[.6, .6], [.7, .7]]
    c = mcolors.ColorConverter().to_rgb
    cm = make_colormap([c('white'), c(mapcolour)])
    for left, top in zip(x, y):
        if top != bottom:
            right = left + width
            plt.imshow(X, interpolation='bicubic', cmap=cm, extent=(left, right, bottom, top), alpha=1, zorder=10)
            plt.plot([left, left], [bottom, top], color='black', linestyle='-', zorder=20)
            plt.plot([right, right], [bottom, top], color='black', linestyle='-', zorder=20)
            plt.plot([right, left], [top, top], color='black', linestyle='-', zorder=20)


# From http://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


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

