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
import matplotlib.colors as mcolors
import itertools

def main(argv):
    parser = argparse.ArgumentParser(description='Read an IgBlastPlus file and plot germline usage.')
    parser.add_argument('infiles', help='input files, names separated by commas. (IgBLASTPlus format, nt or aa)')
    parser.add_argument('field', help='input field (e.g. "V-GENE and allele"')
    parser.add_argument('detail', help='F (family), G (germline) or A (allele)')
    parser.add_argument('-a', '--alpha_sort', help='sort columns alphabetically (default is by decreasing size)', action='store_true')    
    parser.add_argument('-b', '--barcolour', help='colour or list of colours for bars')
    parser.add_argument('-c', '--cols', help='Number of columns for plot')
    parser.add_argument('-co', '--cons', help='Consolidate stats into a single table (csv output only)', action='store_true')
    parser.add_argument('-d', '--dupheader', help='Prefix for duplicate count, eg "DUPCOUNT=" for Presto')
    parser.add_argument('-f', '--frequency', help='Express chart in terms of frequency rather than number of reads', action='store_true')    
    parser.add_argument('-g', '--gradientfill', help='fill bars with a gradiented colour', action='store_true')    
    parser.add_argument('-gh', '--grid_horizontal', help='horizontal grid lines', action='store_true')
    parser.add_argument('-gv', '--grid_vertical', help='vertical grid lines every n bars')
    parser.add_argument('-l', '--limit', help='limit to at most this many most frequent categories')
    parser.add_argument('-s', '--save', help='Save output to file (as opposed to interactive display)')
    parser.add_argument('-sz', '--size', help='Figure size (x,y)')
    parser.add_argument('-t', '--titles', help='titles for each plot, separated by commas')
    parser.add_argument('-w', '--width', help='relative bar width (number between 0 and 1)')
    parser.add_argument('-y', '--ymax', help='Max y-value to use on all charts')

    args = parser.parse_args()
    infiles = args.infiles.split(',')
    if args.detail not in ['F', 'G', 'A']:
        print 'Error: detail must be one of F, G, A'
        exit()
    detail = args.detail
    field = args.field
    consolidate = args.cons
    limit = int(args.limit) if args.limit else None
    ncols = int(args.cols) if args.cols else 1
    frequency = args.frequency
    ymax = float(args.ymax) if args.ymax else None
    outfile = args.save if args.save else None
    titles = args.titles.split(',') if args.titles else infiles
    bar_width = float(args.width) if args.width else 1.0
    mapcolour = args.barcolour if args.barcolour else 'blue'
    mapcolour = mapcolour.split(',')
    grid_vertical = int(args.grid_vertical) if args.grid_vertical else False

    nrows = len(infiles) / ncols
    if len(infiles) % ncols != 0:
        nrows += 1
    (sizex, sizey) = args.size.split(',') if args.size else (8*ncols,4*nrows)

    alpha_sort = args.alpha_sort
    dupheader = args.dupheader

    stats = []
    for infile in infiles:
        stats.append(determine_stats(alpha_sort, detail, dupheader, field, frequency, infile, limit))

    if consolidate:
        fullstats = []
        for infile in infiles:
            fullstats.append(determine_stats(alpha_sort, detail, dupheader, field, frequency, infile, None))
        all_germlines_required = []
        for stat in stats:
            (_, legends) = stat
            for legend in legends:
                if legend not in all_germlines_required:
                    all_germlines_required.append(legend)
        all_germlines_required.sort()
        fullheightlist = []
        for (fullstat, title) in zip(fullstats, itertools.cycle(titles)):
            stat_lookup = {}
            (heights, legends) = fullstat
            for (height, legend) in zip(heights, legends):
                stat_lookup[legend] = height
            all_heights = []
            for germline in all_germlines_required:
                all_heights.append(stat_lookup[germline] if germline in stat_lookup else 0)
            fullheightlist.append(all_heights)

    if not outfile or len(outfile) < 5 or outfile[-4:] != '.csv':
        plt.figure(figsize=(float(sizex),float(sizey)))
        if not consolidate:
            plot_number = 1
            for (stat, title, colour) in zip(stats, itertools.cycle(titles), itertools.cycle(mapcolour)):
                (heights, legends) = stat
                plot_file(heights, legends, frequency, ymax, nrows, ncols, plot_number, title, colour, bar_width, args.gradientfill, args.grid_horizontal, grid_vertical)
                plot_number += 1
        else:
            plot_multi(fullheightlist, all_germlines_required, frequency, ymax, titles, mapcolour, bar_width, args.gradientfill, args.grid_horizontal, grid_vertical)
        plt.tight_layout()
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
            
    else:
        with open(outfile, 'wb') as fo:
            writer = csv.writer(fo)
            if not consolidate:
                for (stat, title) in zip(stats, itertools.cycle(titles)):
                    (heights, legends) = stat
                    writer.writerow([''])
                    writer.writerow([title])
                    writer.writerow(['Germline'] + legends)
                    writer.writerow(['Occurrences'] + heights)
            else:
                writer.writerow(['Germline'] + all_germlines_required)                    
                for (heights, title) in zip(fullheightlist, itertools.cycle(titles)):
                    writer.writerow([title] + all_heights)


def plot_file(heights, legends, frequency, ymax, nrows, ncols, plot_number, title, mapcolour, bar_width, gradientfill, grid_horizontal, grid_vertical):
    x_pos = np.arange(len(heights))
    ax = plt.subplot(nrows, ncols, plot_number)
    plt.xticks(x_pos+0.5, legends, rotation=-70, ha='center')
    ax.tick_params(direction='out', top=False, right=False)
    plt.xlabel(title)
    if ymax:
        plt.ylim(0, ymax)
    if frequency:
        plt.ylabel('Frequency')
    else:
        plt.ylabel('Reads')

    plt.xlim(0, len(heights))
    
    bar_pos = x_pos
    
    if bar_width < 1.:
        bar_pos = bar_pos + (1-bar_width)/2.
    
    if grid_horizontal:
        plt.grid(which='major', axis='y', c='black', linestyle='-', alpha=0.6, zorder=1)
        
    if grid_vertical:
        pos = grid_vertical
        while pos < len(x_pos):
            plt.plot([x_pos[pos], x_pos[pos]], [0, ymax], c='black', linestyle='-', alpha=0.6, zorder=1)
            pos += grid_vertical

    if gradientfill:
        gbar(bar_pos, heights, mapcolour, width=bar_width)
    else:
        plt.bar(bar_pos, heights, width=bar_width, color=mapcolour, zorder=10)
        
    # Remove every other y label because we get far too many by default
    
    locs, labels = plt.yticks()
    newlocs = []
    newlabels = []
    
    for i in range(0, len(labels)):
        if i % 2 != 0:
            newlocs.append(locs[i])
            if frequency:
                newlabels.append(str(float(locs[i])))
            else:
                newlabels.append(str(int(locs[i])))
            
    plt.yticks(newlocs, newlabels)

    ax.set_aspect('auto')
    plt.tight_layout()


def plot_multi(heightlist, legends, frequency, ymax, titles, mapcolour, bar_width, gradientfill, grid_horizontal, grid_vertical):
    x_pos = np.arange(len(legends))
    ax = plt.subplot(1, 1, 1)
    plt.xticks(x_pos+0.5, legends, rotation=-70, ha='center')
    ax.tick_params(direction='out', top=False, right=False)
    if ymax:
        plt.ylim(0, ymax)
    if frequency:
        plt.ylabel('Frequency')
    else:
        plt.ylabel('Reads')

    plt.xlim(0, len(legends))
    
    bar_pos = x_pos
    
    if bar_width < 1.:
        bar_pos = bar_pos + (1-bar_width)/2.
    
    if grid_horizontal:
        plt.grid(which='major', axis='y', c='black', linestyle='-', alpha=0.6, zorder=1)
        
    if grid_vertical:
        pos = grid_vertical
        while pos < len(x_pos):
            plt.plot([x_pos[pos], x_pos[pos]], [0, ymax], c='black', linestyle='-', alpha=0.6, zorder=1)
            pos += grid_vertical

    bar_width = bar_width/len(heightlist)
    i = 0
    for heights,colour in zip(heightlist, itertools.cycle(mapcolour)):
        if gradientfill:
            gbar(bar_pos + i*bar_width, heightlist[i], colour, width=bar_width)
        else:
            plt.bar(bar_pos + i*bar_width, heightlist[i], width=bar_width, color=colour, zorder=10)
        i += 1
        
    # Remove every other y label because we get far too many by default
    
    locs, labels = plt.yticks()
    newlocs = []
    newlabels = []
    
    for i in range(0, len(labels)):
        if i % 2 != 0:
            newlocs.append(locs[i])
            if frequency:
                newlabels.append(str(float(locs[i])))
            else:
                newlabels.append(str(int(locs[i])))
            
    plt.yticks(newlocs, newlabels)

    ax.set_aspect('auto')
    plt.tight_layout()


def gbar(x, y, mapcolour, width=1, bottom=0):
    X = [[.6, .6], [.7, .7]]
    c = mcolors.ColorConverter().to_rgb
    cm = make_colormap([c('white'), c(mapcolour)])
    for left, top in zip(x, y):
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

