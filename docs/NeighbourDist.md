# NeighbourDist #

NeighbourDist plots the distribution of nearest-neighbour distances across a set of CDR3 sequences. Where these are taken from a clonally expanded sequence set, the distribution will typically show two peaks: the minimum between these peaks indicates a suitable identity threshold for clustering. The field IMGT-CDR3 is used, and is required to be present in the file. It may contain either amino acid or nucleotide sequences. The command also plots the distribution of CDR3 lengths. For large data sets, the analysis can be restricted to a random sample.

## Usage ##

     NeighbourDist.py [-h] [-l LIMIT] [-v] [-i] [-g LENGTH_LIMS]
                      [-d DIST_LIMS]
                      infile outprefix

Argument|Meaning
---------|-------
`infile`|Input file names (IMGT/IgBLASTPlus format).
`outprefix`|Prefix to use on output file names
`[-l LIMIT]`|Randomly select at most LIMIT sequences from the file
`[-v]`|Display verbose status output
`[-i]`|Display the plots interactively as well as writing them to files
`[-g LENGTH_LIMS]`|4 comma-separated numbers specifying ranges for the length distribution plot (xmin,xmax,ymin,ymax)
`[-d DIST_LIMS]`|4 comma-separated numbers specifying ranges for the nearest neighbour plot (xmin,xmax,ymin,ymax)
`[-h]`|Provide a help message and exit

## Output Files ##

(prefix)_CDR3_length.pdf - length distribution
(prefix)_CDR3_min_dist.pdf - nearest neighbour distribution

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data.