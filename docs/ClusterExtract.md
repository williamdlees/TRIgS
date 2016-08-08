# ClusterExtract #

Extract IMGT-style analysis records for all members of a cluster from a cluster file produced by CD-HIT or ClusterSeqs.py.

## Usage ##

    usage: ClusterExtract.py [-h] [-i] id clstfile imgtfile outfile

    Given the ID of a cluster member, extract records for all members of that
    cluster from an IMGT-style file.

    positional arguments:
      id                 ID of the cluster member
      clstfile           cluster file (CD-HIT format)
      imgtfile           file from which to extract records (IMGT, IgBLASTPlus
                         format)
      outfile            output file (IMGT, IgBLASTPlus format)

    optional arguments:
      -h, --help         show this help message and exit
      -i, --ignore_size  Ignore size designations in IDs (if there is a semicolon
                     in the ID, it and any following text will be ignored)

# ClusterStats #

Report statistics of a cluster file produced by CD-HIT or ClusterSeqs.py, broken down by sample.

	positional arguments:
	  clstfile              cluster file (CD-HIT format)
	  samples               comma-separated list of sample names
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -s SAVE, --save SAVE  Save graphical output to file (as opposed to interactive
	                        display)

`samples` is a comma-separated list of sample identifiers: these are strings that match against the FASTA 
ids uniquely, and are used to identify sequences from multiple samples concatenated into the same FASTA file.
Typically they will be suffixes or prefixes.

In the output, `Span` denotes the number of samples that have representatives within a given cluster.

Example output:

	>ClusterStats.py heavy_n_junction.clstr A1,A2,A3,TP,TS,TLN
	
	Number of clusters with more than one member, and at least one member in one of the nominated samples: 28202
	Largest cluster size in that set: 10436
	
	Sample  Unique  Shared  Total   Max Size        Gini Index
	A1      1067    1209    2276    1426            0.58
	A2      613     846     1459    2769            0.64
	A3      1381    2037    3418    2822            0.54
	TP      67      503     570     1127            0.59
	TS      4862    2971    7833    5168            0.66
	TLN     9069    199     9268    1031            0.36
	
	Span    Occurrences
	1       17059
	2       8550
	3       1822
	4       523
	5       196
	6       52
	

![Image](https://rawgit.com/williamdlees/BioTools/master/docs/span.png)