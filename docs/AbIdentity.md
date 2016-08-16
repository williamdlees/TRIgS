# AbIdentity #

Identity/Divergence plots show, for a set of sampled sequences, their divergence from the inferred germline and identity to a target sequence of interest. AbIdentity calculates the necessary identities. All sequences are presented in IgBlastPlus format. Only productive sequences that match the specified germline and are full-length (i.e., have at least some sequence in both FR1 and FR4) are processed. For each sequence, the separate regions (FR1, CDR1, etc) are aligned with the target sequence, following which the number of nucleotides in agreement is counted. For FR1 and FR4, sequences are trimmed to the shortest length present in both the target sequence and the sequence under evaluation. The overall identity is expressed as a percentage of nucleotides in agreement. A similar process is followed to determine divergence from the germline, except here the identity is only determined over regions FR1-FR3.


## Usage ##

	AbIdentity.py [-h] [-u] [-v]
                  queryfile queryseq datafile germline germlinefile tag

Argument|Meaning
---------|-------
`queryfile`|File containing the target sequence (IGBLASTPlus format)
`queryseq`|Sequence ID of the target sequence.
`datafile`|File containing sequences to compare against the target (IGBLASTPlus format)
`germline`|Germline of target (only sequences of this germline or its alleles will be compared)'
`germlinefile`|File containing the germline sequence (IgBLASTPlus format)
`outfile`|Output file name (CSV format)
`[-u]`|Include unproductive sequences (the match against germline may be useful, but the match against target will not be)
`[-v]`|Provide field-by-field output to stdout for the first few sequences
`[-h]`|Provide a help message and exit

## Output File ##

The output file is in CSV format with headers SequenceId, TargetDist and GermlineDist. TargetDist gives the percentage identity of the sequence with the target sequence. GermlineDist gives the percentage divergence (100-identity) of the sequence's V-region compared to the germline.

## Notes ##

The germline sequence ID may be specified either in the format used in IMGT germline files, (e.g., M93173|IGHV1S40&#42;01|Oryctolagus), or as a simple germline name (IGHV1S40&#42;01). The sequence ID of the germline in the germline file must match the specification. The germline name in the 'V-GENE and allele' column should always take the simple form. This slightly contorted approach enables IMGT germline files to be used without the need to rewrite their sequence IDs.

When determining which sequences to compare against target and germline, the trailing **nn*s are removed before comparing 'V-GENE and Allele' against the germline - hence all alleles will be matched.

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data.