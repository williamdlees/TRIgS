# Tools for FASTA file Manipulation

This is a collection of small utilities that have proved useful in RepSeq analytic pipelines.

## CountRecords


    usage: CountRecords.py [-h] dupheader infile

    Count FASTA records, accounting for duplicate counts in the header

    positional arguments:
      dupheader   Prefix for duplicate count, eg "DUPCOUNT=" for Presto, "size="
                  for usearch
      infile      input file (FASTA)

    optional arguments:
      -h, --help  show this help message and exit

## FastaMatch


    usage: FastaMatch.py [-h] [-s] [-e] pattern [infile] [outfile]

    Copy FASTA records that match the specified pattern to the output file.

    positional arguments:
      pattern     pattern to match (regex
      infile
      outfile

    optional arguments:
      -h, --help  show this help message and exit
      -s, --seq   copy records whose sequences match (default is to match the
                  header)
      -e, --exc   copy records that do not match instead of those that do

## FastaSample

    usage: FastaSample.py [-h] percent [infile] [outfile]

    Provide a (randomized) percentage sample of the FASTA file.

    positional arguments:
      percent     percentage required (0-100)
      infile
      outfile

    optional arguments:
      -h, --help  show this help message and exit

## FastaUniq

    usage: FastaUniq.py [-h] [infile] [outfile]

    If records are found with identical IDs, remove all but one on the assumption
    that they are duplicates.

    positional arguments:
      infile
      outfile

    optional arguments:
      -h, --help  show this help message and exit

## FastaSampleUniq

    usage: FastaSampleUniq.py [-h] [--plot]
                              sample_size iterations dupheader [infile]

    Find the average number of unique sequences in a set of random samples taken from the
    input file, accounting for duplicate counts in the header

    positional arguments:
      sample_size  Number of records to sample
      iterations   Number of times to sample
      dupheader    Prefix for duplicate count, eg "DUPCOUNT=" for Presto
      infile       input file (FASTA)

    optional arguments:
      -h, --help   show this help message and exit
      --plot       plot the distribution of the number in each sample



