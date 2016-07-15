# ClusterSeqs #

ClusterSeqs partitions a set of sequences into clusters, such that all sequences in each cluster are the same length, and, for each sequence in a cluster, there is at least one other sequence in the cluster within the threshold identity. The identity is calculated as the Hamming distance divided by length. 

## Usage ##

     ClusterSeqs.py [-h] [-l LIMIT] [-m] [-t THREADS] [-c] [-d] [-v]
                    infile outfile timelinefile cutoff labels

Argument|Meaning
---------|-------
`infile`|Input file name (IMGT/IgBLASTPlus format).
`infile`|Input file names (IMGT/IgBLASTPlus format).
`timelinefile`|A count, for each cluster, of the number of sequences in each sample.
`cutoff`|A number between 0 and 1: sequences will be clustered if they differ in identity by less than this number.
`labels`|A comma-separated list of strings identifying samples. Each Sequence ID is assumed to match just one label, identifying the sample from which it originated. Labels are used to created the timeline file: this field may be set to a random value such as XXX, in which case the analysis will complete correctly except that the timeline file will not be populated. Alternatively, it may be set to a value present in every Sequence ID, in which case the timline file will list the number of sequences in each cluster. 
`[-l LIMIT]`|Randomly select at most LIMIT sequences from the file
`[-c]`|Perform an independent check after processing, to confirm that the clusters have been correctly formed.
`[-d]`|Dump the cluster structure to clusters.pickle for debug purposes
`[-v]`|Display verbose status output
`[-h]`|Provide a help message and exit

## Output Files ##

outfile - list of clusters, in [CD-HIT](http://weizhongli-lab.org/cd-hit/ "CD-HIT") format
timelinefile - list of clusters, showing the number of members at each sample timepoint, using sample labels defined by `labels`

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data. The -c option provides independent confirmation that the partitioning has been successful.