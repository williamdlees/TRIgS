# ClusterSeqs #

ClusterSeqs partitions a set of sequences into clusters, such that all sequences in each cluster are the same length, and form a connected network in which each sequence has neighbours that are within the identity threshold. The identity is calculated as the Hamming distance divided by length. 

Typically, clustering methods first create a distance matrix, exhaustively enumerating the pairwise distances between items to be clustered. Once this is in place, there are efficient techniques for clustering, however the distance matrix can take a long time to compute. ClusterSeqs clusters chunks of sequences, and then merges the chunks. Because in NGS analyses the number of clusters is generally much smaller than the number of sequences, this can greatly reduce the number of distances that need to be computed (although the performance of the method will depend on the nature of the clusters).

ClusterSeqs requires [SciPy](http://scipy.org) and [python-Levenshtein](https://pypi.python.org/pypi/python-Levenshtein).

## Usage ##

     ClusterSeqs.py [-h] [-l LIMIT] [-m] [-t THREADS] [-c] [-d] [-u] [-v]
                    infile outfile timelinefile cutoff labels

Argument|Meaning
---------|-------
`infile`|Input file name (IMGT/IgBLASTPlus format).
`timelinefile`|A count, for each cluster, of the number of sequences in each sample.
`cutoff`|A number between 0 and 1: sequences will be clustered if they differ in identity by less than this number.
`labels`|A comma-separated list of strings identifying samples. Each Sequence ID is assumed to match just one label, identifying the sample from which it originated. Labels are used to created the timeline file: this field may be set to a random value such as XXX, in which case the analysis will complete correctly except that the timeline file will not be populated. Alternatively, it may be set to a value present in every Sequence ID, in which case the timline file will list the number of sequences in each cluster. 
`[-l LIMIT]`|Randomly select at most LIMIT sequences from the file
`[-c]`|Perform an independent check after processing, to confirm that the clusters have been correctly formed.
`[-d]`|Dump the cluster structure to clusters.pickle for debug purposes
`[-u]`| Remove duplicate sequences: only the first is retained
`[-v]`|Display verbose status output
`[-h]`|Provide a help message and exit

## Output Files ##

outfile - list of clusters, in [CD-HIT](http://weizhongli-lab.org/cd-hit/ "CD-HIT") format
timelinefile - list of clusters, showing the number of members at each sample timepoint, using sample labels defined by `labels`

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data. The -c option provides independent confirmation that the partitioning of any data set has been successful.