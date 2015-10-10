# AnnotateTree #

Given a set of nucleotide sequences and an inferred phylogenetic tree, create a set of reports:

- Phylogenetic trees showing inferred ancestors and amino acid substitutions
- Nucleotide and amino acid alignments of input sequences and inferred ancestors
- Frequency-based AA logo plot, based on input sequences

Sequence numbering can be controlled by the user, allowing both deletions and insertions, using either PDB- or IMGT-style insertion labels. Optionally, if the sequences are Ig variable region sequences, the user may specify the CDR locations, in which case additional reports are created showing the AA variability within each CDR. 

## Pre-Requisites ##

<a href="https://www.python.org/downloads">Python 2.7.7, or a later version of Python 2.</a>
          
<a href="http://biopython.org">BioPython.</a>
          
<a href="http://evolution.genetics.washington.edu/phylip.html">PHYLIP</a>. 
          
<a href="http://etetoolkit.org/">The ETE Toolkit</a> and its dependency <a href="http://www.riverbankcomputing.co.uk/software/pyqt/download">PyQt4.</a>
          
Optionally, <a href="http://weblogo.berkeley.edu">Berkeley Weblogo</a> and <a href="http://www.ghostscript.com">Ghostscript</a> may be installed in order to generate a frequency plot of input sequences.

## Installation Notes ##

The dnaml executable from PHYLIP must be copied to the directory containing the analysis toolset, or its location included in the windows or linux PATH environment variable.

On Linux (but not on Windows), an X server is required to render the phylogenetic trees. Please see the notes on installing an X server <a href="http://pythonhosted.org/ete2/tutorial/tutorial_webplugin.html">here</a>, if ETE reports that it cannot connect to the X server. If necessary, you can provide a standalone X server such as xdm, as suggested in the ETE documentation <a href="http://pythonhosted.org/ete2/tutorial/tutorial_webplugin.html">here</a>. In this case, you will also need to configure the DISPLAY environment variable to point to the X server: set DISPLAY=localhost:0.0

## Usage ##

    python AnnotateTreeCmd.py seqnumfile seqfile treefile cdrfile tag wd

Argument|Meaning
---------|-------
`seqnumfile`|A file defining the sequence numbering (see 'File Formats' below).
`seqfile`|A FASTA file containing the nucleotide sequences to be analysed. Sequences must all have the same length and must represent a whole number of valid codons. Whole-codon gaps (represented by ---) are allowed. The number of sequences must match the number of labelled nodes in the treefile (see below) and the FASTA labels must identically match the treefile's node labels. The first sequence in the file must represent the root. In the case of antibody sequences, this would normally be the germline, which may if necessary be inferred using the companion script [RevertToGermline](RevertToGermline.md).
|`treefile`|The input tree, in Newick format. The pipeline will root the tree on the first node listed in `seqfile`: the tree as supplied therefore does not need to be rooted.
|`cdrfile`|A file listing the locations of the CDRs (see 'File Formats' below)
|`tag`|A string which will be prepended to the filename of each output file produced by the script. The output files otherwise have fixed file names, listed below.
|`wd`|The relative or absolute pathaname of the working directory that the pipeline should use. This directory must exist at the time the script is called. Output files and working files will be created in this directory. You are recommended to use a separate working directory for each analysis, so that working files for each analysis can be examined. Please note that, on the command line, pathnames to other files (seqfile, treefile etc.) must be specified relative to the directory in which the script is run, not relative to the working directory.  

## Output Files ##

Filename|Contents
--------|--------
tag_aa_alignment.fa:|AA translation of input sequences and inferred ancestral intermediates, in FASTA format.|
tag_aa_alignment.txt:|AA translation of input sequences and inferred ancestral intermediates, in pretty print format.
tag_nt_alignment.fa:|The input nt sequences and inferred ancestral intermediates, in FASTA format.
tag_annotated_treefile.new:|The input tree, in Newick format, annotated with AA transitions. The transitions are provided as node labels.
tag_annotated_treefile.png:|The above tree, rendered in .png format.
tag_annotated_treefile.svg:|The above tree, rendered in .svg format.
tag_annotated_treefile_tot.new:|The input tree, in Newick format, annotated with the number of AA transitions in each branch.
tag_annotated_treefile_tot.png:|The above tree, rendered in .png format.
tag_annotated_treefile_tot.svg:|The above tree, rendered in .svg format.
tag_annotated_treefile_sum.new:|The input tree, in Newick format, annotated with a summary of transitions in each FR and CDR*.
tag_annotated_treefile_sum.png:|The above tree, rendered in .png format*.
tag_annotated_treefile_sum.svg:|The above tree, rendered in .svg format*.
tag_cdr_analysis.html:|The CDR analysis, as an HTML table*.
tag_intermediates_treefile.new:|The input tree, in Newick format, with intermediate nodes labelled. The main purpose of this tree is to indicate the relative position of each inferred ancestral intermediate in the tree.
tag_intermediates_treefile.png:|The intermediate tree, rendered in .png format.
tag_intermediates_treefile.svg:|The intermediate tree, rendered in .svg format.
tag_aa_logo.png:|A logo frequency plot of the AA translation of all input sequences EXCEPT the first (root) sequence (the root sequence is assumed to have been added post facto and not part of the sequenced sample).

 *starred files are only produced if the CDR positions are defined in cdrfile.

## File Formats ##

### `seqnumfile` ###

Defines the amino-acid numbering scheme for the sequences. This can contain insertions (e.g. 99A, 99B) as well as deletions, following the numbering used in PDB files or IMGT files. The numbering scheme is used in the AA alignment, in annotations to the tree, and in identifying CDRs for analysis if required.

The first line of the file contains a single number, representing the number of the first AA in the sequence, or alternatively a number and letter separated by a comma, e.g. 89,A. Subsequent lines indicate deviations from an incremental numbering scheme taking one of the following forms:

Line|Meaning
----|-------
`<number>, -`|A number followed by a minus sign indicates that the number should be skipped.
`<number>, <single letter or multiple numbers>`|The numbered residue is followed by an insertion. Hence 99,A would indicate that residue 99 is followed by residue 99A. A subsequent 99, B would indicate that 99A is followed by 99B.
`<number>; <single letter or multiple numbers>`|The numbered residue is PRECEDED by an insertion (as in the IMGT junction)

#### *Example 1*  

    56

The first AA in the file will be numbered 56, and subsequent AAs will be numbered 57, 58, 59, etc.

#### *Example 2*

             1,A
             3,-
             6,A
             6,B
             
The sequence will be numbered 1A,2,4,5,6,6A,6B,7,8,9,....


#### *Example 3 - showing the use of insertions preceding the ordinal, as used by IMGT*
             1
             2,1
             2,2
             3;2
             3;1
             3

The sequence will be numbered 1,2.1,2.2,3.2,3.1,3,4... 

### `cdrfile` ###

A file that optionally determines the starting and ending position of each CDR, using the numbering scheme above. This is represented as six positions on a single line, representing the starting and ending position of CDRs 1,2 and 3 respectively. The sixth number (representing the upper bound of CDR3) is allowed to be higher than the highest position in the sequences, in which case CDR3 will be taken to run from its lower position (as indicated in the file) through to the last position in the sequence. The sequences are not required to span the entire range of positions specified in the cdrfile: for example the lowest position in the sequences could be 56, even though CDR1 was specified as spanning positions 27-38. The file may be empty, in which case CDR analysis is not performed and a warning is raised: other pipeline outputs will however be produced. Use an empty file for annotating non-Ab sequences.

#### *Example 1*

    27,38,56,65,105,140

This example follows the standard IMGT numbering scheme

#### `Example 2`

    1,6A,30B,35,95C,120

This example demonstrates the use of 'insertion' positions

## Configuration Files ##

dnaml.ctl (in the same directory as AnnotateTree.py) specifies the input selections for dnaml and may be changed if desired.

## Working Files ##

dnaml input, output and intermediate files are left in the directory and may be referred to if necessary - for example if a dnaml error is reported. In particular, rst.txt contains details of the ancestral reconstruction, and dnaml.txt is the logfile for the dnaml run. Please refer to the dnaml documentation for further details.

tag_alignment_for_logo.fa contains the AA sequences used for the logo. 

## Testing ##

Before use, verify that the pre-requisites listed are installed. Verfy that dnaml is available by opening a command prompt, changing to the directory containing the utilities, and typing dnaml or ./dnaml. If you wish to create logo plots, verify that Berkely Weblogo is installed by typing seqlogo or ./seqlogo.

Executing the commands in test/test.sh will exercise AnnotateTree and RevertToGermline with sample data. Output AnnotateTree will be created in a subdirectory.

