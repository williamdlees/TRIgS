# IgBLASTPlus #

NCBI's <a href="http://www.ncbi.nlm.nih.gov/igblast/">IgBLAST</a> is a valuable tool for high-speed Ig sequence analysis which has the advantage that it can be run locally, but it provides verbose output that is not amenable to computational processing. Also, although it matches the junction against D and J segments providing best match analysis, it does not identify the 3' end of the junction and hence does not list the CDR3 or the junction sequence. IgBLASTPlus summarises IgBLAST output in an IMGT-style tab-separated format, identifies the conserved W or F that terminates the junction by comparison with the J germline, and hence is able to list the CDR3 and junction sequences.

## Usage ##

IgBLAST must be run with the option '-outfmt 3' to create the analysis file required by IgBLASTPlus.

	python IgBLASTPlus.py [-h] [-v] seqfile jgermfile igblastfile tag

Argument|Meaning
---------|-------
`seqfile`|The sequences analysed by IgBLAST (FASTA format).
`jgermfile`|The J-germline file provided to IgBLAST (FASTA format).
`igblastfile`|The output file provided by IgBLAST (with option '-outfmt 3')
`tag`|A string which will be prepended to the filename of each output file produced by the script. The output files have fixed file names, listed below.
`[-v]`|Provide diagnostic output when parsing IgBLAST alignments.
`[-h]`|Provide a help message and exit

IgBLASTPlus provides two output files whose names are prepended with the `tag` specified on the command line. *tag*_n.txt contains a nucleotide-level analysis, while *tag*_aa.txt provides the equivalent analysis with nucleotide sequences translated to amino acids. Both files are tab-separated. Column names match IMGT names, for example FR1-IMGT, CDR1-IMGT. The output includes two columns that are not provided by IMGT analysis: 'Chain Type', which includes the chain type determined by IgBLAST, and FR4-IMGT, which is deduced by IgBLASTPlus junction analysis, and contains the sequence running from the conserved W or F to the 3' end of the aligned J-gene.

## Testing ##

`test/Test_IgBlastPlus.py` exercises the script with test input and checks the result.