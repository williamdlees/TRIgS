# ClusterExtract #

Extract IMGT-style analysis records for all members of a cluster.

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


