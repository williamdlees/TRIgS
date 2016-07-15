# ClusterGraph #

ClusterGraph takes a cluster file produced by ClusterSeqs or by [CD-HIT](http://weizhongli-lab.org/cd-hit/ "CD-HIT") and produces edge and vertex csv files suitable for import by [Gephi](https://gephi.org/ "Gephi"). 

## Usage ##

     ClusterGraph.py [-h] [-l LIMIT]
                     infile clstfile outfile_prefix colourfile idprefix
                     highlightcolour cutoff

Argument|Meaning
---------|-------
`infile`|Sequence file name (FASTA format).
`clstfile`|Cluster file name (CD-HIT format).
`outfile_prefix`|Prefix for output files
`colour file`|Colour scheme
`idprefix`|Prefix of 'generic' sequence IDs
`highlightcolour`|Colour to use for non-generic sequences
`cutoff`|Distance cutoff for adjacency (expressed as a fraction between 0 and 1)
`[-l LIMIT]`|Read at most LIMIT sequences from the file (intended for testing purposes)
`[-h]`|Provide a help message and exit

## Description ##

`infile` and `clstfile` must contain an identical set of sequences. The majority of sequences are expected to conform to a naming format such that they all share a common prefix `idprefix` of at least one letter: eg S1_123, S2_4356 etc. These sequences are referred to as 'generic', and will be assigned a colour according to a naming scheme. The naming scheme allows a lower level of granularity: in this example, S1 and S2 could originate from different samples, and be assigned different colours. Sequences whose names do not have the common prefix are considered to be special, and are assigned the special highlight colour.

An edge will be created between two sequences provided that their identities (measured by hamming distance) differ by less than the cutoff (if the sequences are different lengths, an edge will not be created). This allows connected graphs to be drawn that are more discriminating than the clusters defined in the input file.

## Input Files ##

`colourfile` defines the batches and the colours. Example:

    A1,red
    A2,yellow
    A3,green

The first item on each line identifies the class, and the second item the colour of that class. The program searches each sequence ID for each class in turn. The sequence is assigned the colour corresponding to the first matching class. If no classes match, the sequence is assigned the colour black.  

Colours in `colourfile` and on the command line may be specified in any format acceptable to Gephi. 

## Output Files ##

(prefix)_edges.csv, (prefix)_verts.csv - edge and vertex files in Gephi csv format.

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data.