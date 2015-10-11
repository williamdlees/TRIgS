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

# This is the command-line driver for the annotation pipeline.
# The pipeline is run in the current directory. Usually it will be best to analyse each clonal family in a separate directory.
# The pipeline takes, as principal inputs, a set of aligned nucleotide sequences and an inferred phylogenetic tree. It
# produces as its main outputs:
#   - The inferred ancestral sequence of each internal node in the tree
#   - An annotated tree, in which internal and external node names are supplemented with the amino acid transitions
#     as deduced from the inferred ancestral states
#   - Optionally, in the case that the sequences represent antibody heavy or light chain sequences, an analysis of the
#     variation at each location in the CDRs, showing those that are common to germline throughout the tree, those that
#     share a common residue across all sequences except the germline, and those that exhibit variation across leaf nodes.
#
# Key dependencies (see readme.txt for installation notes):
#
# PHYLIP dnaml (http://evolution.genetics.washington.edu/phylip.html). Dnaml is used for ancestral reconstruction.
# 
# BioPython (http://biopython.org). BioPython is used for sequence manipulation and i/o.
# 
# ETE Toolkit (http://etetoolkit.org). The ETE Toolkit is used to render phylogentic trees. 
#
# Weblogo (http://http://weblogo.berkeley.edu/logo.cgi) and its ghostscrpt dependency. Weblogo must be installed such
# that the command-line seqlogo program is on the path and capable of rendering .png files via ghostscript. The web
# interface is not required. If seqlogo is not available, the logo file will not be produced.
#
# Usage: python AnnotateTreeCmd.py seqnumfile seqfile treefile cdrfile tag wd
#
# Command-line Arguments:
#
# seqnumfile: Defines the amino-acid numbering scheme for the sequences. This can contain insertions (e.g. 99A, 99B) as well
#             as deletions, following the numbering used in PDB files or IMGT files. The numbering scheme is used in the AA 
#             alignment, in annotations to the tree, and in identifying CDRs for analysis if required.
#
#             The first line of the file contains a single number, representing the number of the first AA in the sequence, 
#             or alternatively a number and letter separated by a comma, e.g. 89,A.
#             Subsequent lines indicate deviations from an incremental numbering scheme taking one of the following forms:
#               <number>, - a number followed by a minus sign indicates that the number should be skipped.
#               <number>, <single letter or multiple numbers> - the numbered residue is followed by an insertion
#                - hence 99,A would indicate that residue 99 is followed by residue 99A
#                - 100,11 would indicate that residue 100 is followed by residue 100.11
#                - a subsequent 99, B would indicate that 99A is followed by 99B
#               <number>; <single letter or multiple numbers> - the numbered residue is PRECEDED by an insertion (as in the IMGT junction)
#
#             Example 1 - a single line in the file, containing the number 56. The first AA in the file will be numbered 56, and
#             subsequent AAs will be numbered 57, 58, 59, etc.
#
#             Example 2 - a four line file:
#             1,A
#             3,-
#             6,A
#             6,B
#             The sequence will be numbered 1A,2,4,5,6,6A,6B,7,8,9,....
#
#             Example 3 - showing the use of insertions preceding the ordinal:
#             1
#             2,1
#             2,2
#             3;2
#             3;1
#             3
#             The sequence will be numbered 1,2.1,2.2,3.2,3.1,3,4... 
#
# seqfile;    A FASTA file containing the nucleotide sequences to be analysed. Sequences must all have the same length and must 
#             represent a whole number of valid codons. Whole-codon gaps (represented by ---) are allowed. The number of sequences
#             must match the number of labelled nodes in the treefile (see below) and the FASTA labels must identically match
#             the treefile's node labels. The first sequence in the file must represent the root. In the case of antibody sequences, this
#             would normally be the germline, which may if necessary be inferred using the companion script GermlineFromIMGT.py.
#
# treefile:   The input tree, in Newick format. The pipeline will root the tree on the first node listed in seqfile: the tree as supplied
#             therefore does not need to be rooted.
#
# cdrfile:    A file that optionally determines the starting and ending position of each CDR, using the numbering scheme above. This
#             is represented as six positions on a single line, representing the starting and ending position of CDRs 1,2 and 3 respectively.
#             The sixth number (representing the upper bound of CDR3) is allowed to be higher than the highest position in the sequences, 
#             in which case CDR3 will be taken to run from its lower position (as indicated in the file) through to the last position in
#             the sequence.  The sequences are not required to span the entire range of positions specified in the cdrfile: for example
#             the lowest position in the sequences could be 56, even though CDR1 was specified as spanning positions 27-38. 
#             
#             The file may be empty, in which case CDR analysis is not performed and a warning is raised: other pipeline outputs
#             will however be produced. Use an empty file for annotating non-Ab sequences.
#
#             Example 1: 27,38,56,65,105,140   - this example follows the IMGT numbering scheme
#
#             Example 2: 1,6A,30B,35,95C,120   - this example demonstrates the use of 'insertion' positions
#
# tag:        This is a string, which will be prepended to the filename of each output file produced by the script. The output files 
#             otherwise have fixed file names, listed below.
#
# wd:         The relative or absolute pathaname of the working directory that the pipeline should use. This directory must exist at the time
#             the script is called. Output files and working files will be created in this directory. You are recommended to use a separate 
#             working directory for each analysis, so that working files for each analysis can be examined. Please note that, on the command line, 
#             pathnames to other files (seqfile, treefile etc.) must be specified relative to the directory in which the script is run, not 
#             relative to the working directory.  
#
# Output Files:
#
# tag_aa_alignment.fa:            AA translation of input sequences and inferred ancestral intermediates, in FASTA format.
#
# tag_aa_alignment.txt:           AA translation of input sequences and inferred ancestral intermediates, in pretty print format.
#
# tag_nt_alignment.fa:            The input nt sequences and inferred ancestral intermediates, in FASTA format.
#
# tag_annotated_treefile.new:     The input tree, in Newick format, annotated with AA transitions. The transitions are provided as node labels.
#
# tag_annotated_treefile.png:     The above tree, rendered in .png format.
#
# tag_annotated_treefile.svg:     The above tree, rendered in .svg format.
#
# tag_annotated_treefile_tot.new: The input tree, in Newick format, annotated with the number of AA transitions in each branch.
#
# tag_annotated_treefile_tot.png: The above tree, rendered in .png format.
#
# tag_annotated_treefile_tot.svg:  The above tree, rendered in .svg format.
#
# tag_annotated_treefile_sum.new: The input tree, in Newick format, annotated with a summary of transitions in each FR and CDR*.
#
# tag_annotated_treefile_sum.png: The above tree, rendered in .png format*.
#
# tag_annotated_treefile_sum.svg:  The above tree, rendered in .svg format*.
#
# tag_cdr_analysis.html:          The CDR analysis, as an HTML table*.
#
# tag_intermediates_treefile.new: The input tree, in Newick format, with intermediate nodes labelled. The main purpose of this tree is to indicate
#                                 the relative position of each inferred ancestral intermediate in the tree.
#
# tag_intermediates_treefile.png: The intermediate tree, rendered in .png format.
#
# tag_intermediates_treefile.svg: The intermediate tree, rendered in .svg format.
#
# tag_aa_logo.png:                A logo frequency plot of the AA translation of all input sequences EXCEPT the first (root) sequence (the root sequence
#                                 is assumed to have been added post facto and not part of the sequenced sample).
#
# * - starred files are only produced if the CDR positions are defined in cdrfile.
#
#  Working Files:
#
# dnaml input, output and intermediate files are left in the directory and may be referred to if necessary - for example if a dnaml error is reported.
# In particular, rst.txt contains details of the ancestral reconstruction, and dnaml.txt is the logfile for the dnaml run. Please refer to the dnaml
# documentation for further details.
#
# tag_alignment_for_logo.fa contains the AA sequences used for the logo. 
#
# Configuration Files:
#
# dnaml.ctl (in the same directory as AnnotateTree.py) specifies the input selections for dnaml and may be changed if desired.
  

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import os
import sys
from Dnaml import Dnaml
from Bio import Phylo
from Bio import AlignIO
from Bio import SeqIO
from Alignment import Alignment
from RenderTree import RenderTree
from AnalyseCDR import AnalyseCDR
import gc
import subprocess




def main(argv):
    print "AnnotateTreeCmd v1.0"
    if len(argv) == 2 and argv[1] == '-t':
        conduct_tests()
        exit(0)
    elif len(argv) != 7:
            print 'usage python AnnotateTreeCmd.py seqnumfile seqfile treefile cdrfile tag wd.'
            sys.exit(0)

    for file in argv[1:4]:
        check_file(file)
        
    (seqnumfile, seqfile, treefile, cdrfile, tag, wdir) = argv[1:7]

    if len(cdrfile) > 0:
        check_file(cdrfile)
    else:
        cdrfile = None

    try:
        if not os.path.exists(wdir):
            os.makedirs(wdir)
    except:
        print "Error creating directory %s." % wdir
        sys.exit(0)
        
    try:
        msa = Alignment()
        msa.read_nt(seqfile)    # Check that the sequence comprises a valid set of codons
        for seq in msa:
            if '*' in seq:
                print "Stop codon found in sequence %s." % seq.id
                sys.exit(0)
    except:
        print "Error parsing %s: %s." % (seqfile, sys.exc_info()[1])
        sys.exit(0)
        
    try:
        seq_pos = msa.read_position_numbers(seqnumfile)
    except:
            print "Error parsing %s: %s." % (seqnumfile, sys.exc_info()[1])
            sys.exit(0)        

    if cdrfile is not None:
        try:
            acdr = AnalyseCDR(msa, file_name=cdrfile)
        except:
            print "Error parsing %s: %s." % (cdrfile, sys.exc_info()[1])
            sys.exit(0)

    try:
        seq_align = AlignIO.read(seqfile, "fasta")
    except:
        try:
            seq_align = AlignIO.read(seqfile, "phylip")
        except:
            print "Error parsing %s: %s." % (seqfile, sys.exc_info()[1])
            sys.exit(0)

    try:
        tree = Phylo.read(treefile, "newick")
    except:
        print "Error parsing %s: %s." % (treefile, sys.exc_info()[1])
        sys.exit(0)

    dnaml = Dnaml()
    
    int_aas = dnaml.run_dnaml(seq_align, tree, seq_pos, cdrfile, wdir, report, tag)

    if int_aas is not None:
        try:
            if cdrfile is not None:
                acdr = AnalyseCDR(int_aas, file_name=cdrfile)
                cdr_output = acdr.analyse()
                fo = open(wdir + "/" + tag + "cdr_analysis.html", "w")
                fo.write(cdr_output)
                fo.close()
        except:
            print "Warning: CDRs were not analysed: " + str(sys.exc_info()[1])

        try:
            gc.collect()
            RenderTree.render_annotate(wdir + "/" + tag + "annotated_treefile.new", wdir + "/" + tag + "annotated_treefile.png")
            gc.collect()
            RenderTree.render_annotate(wdir + "/" + tag + "annotated_treefile.new", wdir + "/" + tag + "annotated_treefile.svg")
            gc.collect()
            if cdrfile is not None:
                RenderTree.render_annotate(wdir + "/" + tag + "annotated_treefile_sum.new", wdir + "/" + tag + "annotated_treefile_sum.png")
                gc.collect()
                RenderTree.render_annotate(wdir + "/" + tag + "annotated_treefile_sum.new", wdir + "/" + tag + "annotated_treefile_sum.svg")
                gc.collect()
            RenderTree.render_annotate(wdir + "/" + tag + "annotated_treefile_tot.new", wdir + "/" + tag + "annotated_treefile_tot.png")
            gc.collect()
            RenderTree.render_annotate(wdir + "/" + tag + "annotated_treefile_tot.new", wdir + "/" + tag + "annotated_treefile_tot.svg")
            gc.collect()
            RenderTree.render_annotate(wdir + "/" + tag + "intermediates_treefile.new", wdir + "/" + tag + "intermediates_treefile.png")
            gc.collect()
            RenderTree.render_annotate(wdir + "/" + tag + "intermediates_treefile.new", wdir + "/" + tag + "intermediates_treefile.svg")
            gc.collect()
        except:
            print "Error rendering trees: " + str(sys.exc_info()[1])

        first = True
        orig_recs = []
        for rec in SeqIO.parse(wdir + "/" + tag + "aa_alignment.fa", "fasta"):
            if not first and "node_" not in rec.id:
                orig_recs.append(rec)
            first = False
        
        logo_alignment_file = wdir + "/" + tag + "alignment_for_logo.fa"
        SeqIO.write(orig_recs, wdir + "/" + tag + "alignment_for_logo.fa", "fasta")
        
        with open(wdir + "/" + tag + "weblogo_status.txt", "w") as fo:
            retcode = subprocess.call("seqlogo -f %salignment_for_logo.fa -F PNG -o aa_logo -h 2 -w 20 -acS" % tag, cwd=wdir, shell=True, stdout=fo, stderr=subprocess.STDOUT)
            if retcode == 1:
                fo.write("Trying seqlogo.pl instead.\n")
                retcode = subprocess.call("seqlogo.pl -f %salignment_for_logo.fa -F PNG -o aa_logo -h 2 -w 20 -acS" % tag, cwd=wdir, shell=True, stdout=fo, stderr=subprocess.STDOUT)
            if retcode == 1:
                print "Weblogo not installed: logo plot will not be generated."


def conduct_tests():
    pass

def check_file(filename):
    if not os.path.isfile(filename):
        print "File '%s' does not exist." % filename
        sys.exit(0)

def report(str):
    print str

if __name__=="__main__":
    main(sys.argv)

