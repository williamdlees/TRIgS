# ExtractFromIMGT #

ExtractFromIMGT extracts specified fields from an IMGT or IgBLASTPlus tab-separated file. Results are provided in FASTA format. Multiple fields can be concatenated, for example CDR1-IMGT+FR2-IMGT+CDR2-IMGT. Options allow the output to be limited to complete fields, and to a specific V-germline or allele.

## Usage ##

	python ExtractFromIMGT.py  [-h] [-g V-GERMLINE] [-d D-GERMLINE] [-j J-GERMLINE] [-p] [-u] [-s] [-t CHAINTYPE] 
	                           [-c] [--fr1_lim FR1_LIM] [--fr4_lim FR4_LIM] [-v] 
                               analysis_file required_fields outfile

Argument|Meaning
---------|-------
`analysis_file`|The tab-separated input file produced by IMGT or IgBLASTPlus. The file may include either nucleotide or AA sequences. In the case of IMGT/HighV-Quest, the file should be the gapped or ungapped nt-sequences or AA-sequences file. In the case of IMGT/V-Quest, the corresponding tabs of the Excel file can be saved as tab-separated variable files.
`required_fields`|Column headers of fields to include in the output, in order, separated by + (e.g. CDR1-IMGT+FR2-IMGT)
`outfile`|Pathname of the output file.
`[-g GERMLINE]`|Only records matched against the specified V-germline will be included. If the analysis file lists several possible V-germlines for a sequence, the record will be included if any of them match. The specified match can be a substring, for example `HV1S45*` will match `IGHV1S45*01, IGHV1S45*02`, etc.
`[-d GERMLINE]`|Only records matched against the specified D-germline will be included. If the analysis file lists several possible V-germlines for a sequence, the record will be included if any of them match. The specified match can be a substring, for example `HV1S45*` will match `IGHV1S45*01, IGHV1S45*02`, etc.
`[-j GERMLINE]`|Only records matched against the specified J-germline will be included. If the analysis file lists several possible V-germlines for a sequence, the record will be included if any of them match. The specified match can be a substring, for example `HV1S45*` will match `IGHV1S45*01, IGHV1S45*02`, etc.
`[-p]`|Restrict the output to sequences that are classed productive in the Functionality column
`[-u]`|Exclude extracts that include an unknown (N) nucleotide (applicable to nucleotide sequence files only)
`[-s]`|Exclude extracts that include stop codons (applicable to amino acid sequence files only)
`[-t CHAINTYPE]`|Only include sequences from the specified Chain Type (this is applicable to IgBLASTPlus files only, as IMGT files are always restricted to a specific Chain Type.
`[-c]`|Only include records with fields that are complete. This is achieved by checking that neighbouring fields are nonblank: for example if CDR2-IMGT is required, records will only be included if FR1-IMGT and FR2-IMGT are nonblank.
`[--fr1_lim FR1_LIM]`|Only include records if they contain at least FR1_LIM nucleotides or residues in FR1-IMGT, and truncate the included FR1-IMGT field to that length if it is specified as a required field. This provides a means to trim sequences at the 5' end while restricting output to sequences that are complete at the 5' end.
`[--fr4_lim FR4_LIM]`|Only include records if they contain at least FR4_LIM nucleotides or residues in FR4-IMGT, and truncate the FR4-IMGT field to that length if it is specified as a required field. This provides a means to trim sequences at the 5' end while restricting output to sequences that are complete at the 5' end. This is only applicable to IgBLASTPlus output, as IMGT analysis does not include FR4-IMGT.
`[-v]`|Provide diagnostic output.
`[-h]`|Provide a help message and exit

`outfile` will be a FASTA file containing the required extracted sequences.

## Testing ##

`test/Test_ExtractFromIMGT.py` exercises the script with test data and checks the result.