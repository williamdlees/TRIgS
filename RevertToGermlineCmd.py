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

# Command-line script to infer the germline sequences of a set of chains, given their IMGT nucleotide analysis.
# Each gene is reverted to its germline and aligned with the chain. Other elements (P, N) are assumed not to have
# mutated from the germline. This hypothesis should be checked by examining other members of the clonal
# family.
#
# Usage: python RevertToGermlineCmd.py imgt_nt germline_file species_name outfile option [-mdd]
# where:
#
# imgt_nt is the name of an IMGT nucleotide analysis csv or tab separated file containing analyses of the sequences
# germline_file is the name of the germline library to use (single file in IMGT format). The IMGT germline library
# for all species is downloadable at http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP
# species_name is the species name as used in the germline file, e.g. "Homo sapiens"
# outfile is the name of the output file.
# option is any combination of:
# f - list 'full' germlines where genes are reverted and ns preserved
# v - list germlines where just the v-gene is reverted and other regions gapped
# j - list germlines where all genes are reverted and n-regions gapped
# i - list input sequences
# o - list derived germline for each input sequence
# c - list consensus for derived germlines
# x - write verbose analysis through the report functio 
# -m followed immediately by a number will list V-gene mutations from germline that are seen in all
# sequences of a particular germline, for all germlines where there are at least that number of sequences.
# in all cases, output sequences are aligned on full codon boundaries and incomplete codons are gapped out

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import sys
from GermlineFromIMGT import germline_from_imgt

def main(argv):
    print "RevertToGermlineCmd v1.1"
    if len(argv) not in (6,7):
        print 'usage python RevertToGermlineCmd.py imgt_nt germline_file species_name outfile {ciofvj} [-m<number>].'
        sys.exit(0)
    
    option = argv[5]
    for char in option:
        if char not in 'ciofvjx':
            print 'usage python RevertToGermlineCmd.py imgt_nt germline_file species_name outfile {ciofvj} [-m<number>].'
            sys.exit(0)
            
    fixed_mut = 0
    if len(argv) == 7:
        try:
            fixed_mut = (int)(argv[6].replace("-m", ""))
        except:        
            print 'usage python RevertToGermlineCmd.py imgt_nt germline_file species_name outfile {ciofvj} [-m<number>].'
            
    (sequence_file, germline_lib, species_name, output_file, option) = (argv[1], argv[2], argv[3], argv[4], argv[5])
    germline_from_imgt(sequence_file, germline_lib, species_name, output_file, option, report, error, fixed_mut)

def report(str):
    print str

def error(str):
    print str

if __name__ == "__main__":
    main(sys.argv)

