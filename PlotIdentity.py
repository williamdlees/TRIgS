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
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Alphabet import generic_nucleotide
from Bio import SeqIO
from Bio import Phylo
from itertools import izip

def main(argv):
    parser = argparse.ArgumentParser(description='Create an Identity/Divergence plot.')
    parser.add_argument('repertoire', help='file containing repertoire sequence identities (CSV)')
    parser.add_argument('-a', '--adjust', help='Adjust labels to prevent overlap (requires package adjustText)', action='store_true')
    parser.add_argument('-b', '--bar', help='Plot a colour bar', action='store_true')
    parser.add_argument('-c', '--colourmap', help='colourmap')
    parser.add_argument('-g', '--background', help='Set the contour colourwhere the density is zero')
    parser.add_argument('-mx', '--maxx', help='max divergence value to show')
    parser.add_argument('-my', '--miny', help='min identity value to show')
    parser.add_argument('-p', '--points', help='comma-seperated list of identity files and formats')
    parser.add_argument('-s', '--save', help='Save output to file (as opposed to interactive display)')
    args = parser.parse_args()

    if args.adjust:
        from adjustText import adjust_text
        
    colourmap = args.colourmap if args.colourmap else 'hot_r'

    plist = args.points.split(',')
    points = []
    
    repertoire = read_file(args.repertoire)
    
    def pairwise(iterable):
        a = iter(iterable)
        return izip(a, a)
    
    if len(plist) > 0:
        try:
            for file, format in pairwise(plist):
                points.append((read_file(file), format))
        except IOError:
            print 'file "%s" cannot be opened.' % file  
        except:
            print '"points" must consist of pairs of files and formats.'
            quit()

    max_divergence = int(args.maxx) if args.maxx else None
    min_identity = int(args.miny) if args.miny else None
    savefile = args.save if args.save else None

    if not max_divergence:
        max_divergence = max(repertoire['GermlineDist'])
        for point in points:
            max_divergence = max(max_divergence, max(point[0]['GermlineDist']))
        max_divergence = int(max_divergence) + 1
        
    if not min_identity:
        min_identity = min(repertoire['TargetDist'])
        for point in points:
            min_identity = min(min_identity, min(point[0]['TargetDist']))
        min_identity = int(min_identity)

    H, yedges, xedges = np.histogram2d(repertoire['TargetDist'], repertoire['GermlineDist'], bins=[101-min_identity, max_divergence+1], range=[[min_identity, 101], [-1, max_divergence]], normed=False)

    # For alternative interpolations and plots, see http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram2d.html
    # For colour maps, see http://matplotlib.org/examples/color/colormaps_reference.html

    fig = plt.figure()
    
    cm = plt.cm.get_cmap(colourmap)
    
    if args.background:
        cm.set_under(color=args.background)

    ax = fig.add_subplot(1,1,1)
    im = plt.imshow(H, interpolation='bilinear', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], vmin=0.0000001, cmap=cm)
    ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])
    
    if args.bar:
        cb = plt.colorbar(im, shrink=0.8, extend='neither')
        cb.ax.set_ylabel('sequences', rotation=90)

    texts = []
    
    for point in points:
        df, format = point
        markersize = 5
        label = False
        labelcolour = 'black'
        fmt = format.split('/')
        format = fmt[0]
        for f in fmt[1:]:
            if f[0] == 'm':
                markersize = int(f[1:])
            elif f[0] == 'l':
                label = True
                if len(f) > 1:
                    labelcolour = f[1:]
            else:
                print 'Unrecognised format string: %s' % format
        for index, row in df.iterrows():
            if label:
                if args.adjust:
                    texts.append(plt.text(row['GermlineDist'], row['TargetDist'], row['SequenceId'], bbox={'pad':0, 'alpha':0}, fontdict={ 'color': labelcolour}))
                else:
                    texts.append(plt.text(row['GermlineDist'] + 0.2, row['TargetDist'] - 0.2, row['SequenceId'], bbox={'pad':0, 'alpha':0}, fontdict={ 'color': labelcolour}))
            ax.plot(row['GermlineDist'], row['TargetDist'], format, markersize=markersize)

    if args.adjust:
        adjust_text(texts)

    plt.xlabel('Germline Divergence (%)')
    plt.ylabel('Target Ab Identity (%)')
    
    if savefile:
        plt.savefig(savefile)
    else:
        plt.show()


def read_file(file):
    df = pd.read_csv(file, converters={'SequenceId': lambda x: x})
    for key in ("SequenceId", "TargetDist", "GermlineDist"):
        if key not in df.keys():
            print 'File %s does not contain a column "%s"' % (file, key)
            quit()
            
    for index, row in df.iterrows():
        try:
            x = row[1] * row[2]     # check they behave like numbers
        except:
            print 'Error in file %s: malformed row at %s.' % (file, row[0])

    if len(df) < 1:
        print '%s: empty file.' % file
        quit()

    return df


if __name__=="__main__":
    main(sys.argv)

