# IgBLASTPlus #

NCBI's <a href="http://www.ncbi.nlm.nih.gov/igblast/">IgBLAST</a> is a valuable tool for high-speed Ig sequence analysis which has the advantage that it can be run locally, but it provides verbose output that is not amenable to computational processing. Also, although it matches the junction against D and J segments providing best match analysis, it does not identify the 3' end of the junction and hence does not list the CDR3 or the junction sequence. IgBLASTPlus summarises IgBLAST output in an IMGT-style tab-separated format, identifies the conserved W or F that terminates the junction by comparison with the J germline, and hence is able to list the CDR3 and junction sequences.

## Usage ##

<b>This version of IgBLASTPlus is compatible with IgBLAST v1.6 onwards. Please check your version of IgBLAST before proceeding and upgrade if necessary.</b>

IgBLAST must be run with the option '-outfmt 3' to create the analysis file required by IgBLASTPlus.

	python IgBLASTPlus.py [-h] [-v] seqfile igblastfile tag

Argument|Meaning
---------|-------
`seqfile`|The sequences analysed by IgBLAST (FASTA format).
`igblastfile`|The output file provided by IgBLAST (with option '-outfmt 3')
`tag`|A string which will be prepended to the filename of each output file produced by the script. The output files have fixed file names, listed below.
`[-v]`|Provide diagnostic output when parsing IgBLAST alignments.
`[-h]`|Provide a help message and exit

## Output Files ##

All files are tab-separated with column headers.

### *tag*_n.txt, *tag*_aa.txt ###

Sequences are provided at nucleotide level in *tag*_n.txt and are translated to protein in *tag*_a.txt.

Fields included:

Field|Content
-----|-------
`Sequence ID`|Query label reported by IgBLAST
`Functionality`|`productive` if a valid, in-frame junction was found that is consistent with IgBLAST's junction analysis and `unproductive` otherwise. This is more rigorous than IgBLAST's 'Productive' summary, in that we require FR3 to have been identified in the alignment, and for the identified FR3 to extend to the first Cysteine of the junction.
`V-GENE and allele`|The top-ranking identified V-gene (multiple values are separated by commas if more than one gene was identified with the same top-ranking score).
`D-GENE and allele`|The top-ranking identified J-gene, as above (blank for light chain sequences).
`J-GENE and allele`|The top-ranking identified J-gene, as above.
`Chain Type`|Chain type as reported by IgBLAST.
`Stop Codon`|Presence in the sequence of a stop codon, as reported by IgBLAST.
`V-J Frame`|As reported by IgBLAST
`Strand`|As reported by IgBLAST
`FR1-IMGT, CDR1-IMGT, FR2-IMGT, CDR2-IMGT, FR3-IMGT`|Sequences as reported in the IgBLAST alignment.
`CDR3-IMGT, JUNCTION-IMGT, FR4-IMGT`|Sequences determined by identifying the conserved terminating W (heavy chain or F (light chain).
`Notes`|Notes and warnings provided by IgBLASTPlus during analysis
`nt sequence`|The full nucleotide sequence comprising the query. This is provided in both the nt and the aa file, making it easy to copy and paste the sequence to a different system for consultation

### *tag*_j.txt ###

Nucleotide analysis of the junction, based on IgBLAST's junction details.

Fields included:

Field|Content
-----|-------
`Sequence ID`|Query label reported by IgBLAST
`Functionality`|'productive' if a valid, in-frame junction was found that is consistent with IgBLAST's junction analysis. This is a little more rigorous than IgBLAST's 'Productive' summary, in that we require FR3 to have been identified in the alignment, and for the identified FR3 to extend to the first Cysteine of the junction.
`V-GENE and allele`|The top-ranking identified V-gene (multiple values are separated by commas if more than one gene was identified with the same top-ranking score).
`D-GENE and allele`|The top-ranking identified J-gene, as above (blank for light chain sequences).
`J-GENE and allele`|The top-ranking identified J-gene, as above.
`JUNCTION`|The junction sequence.
`3'V-REGION, N-REGION, D-REGION, N2-REGION, 5'J-REGION`|Sequences identified from IgBLAST's junction details (`D-REGION` and `N2-REGION` are blank for light chain sequences). `3'V-REGION` and `5'J-REGION` are modified from those reported by IgBLAST so that the entire junction is exactly covered.
`Notes`|Notes as reported in the other files. 


## Testing ##

`test/Test_IgBlastPlus.py` exercises the script with test input and checks the result.