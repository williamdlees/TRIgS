# RevertToGermline #

Infer the germline sequences of a set of Ig variable region sequences, given their IMGT nucleotide analysis. Each gene is reverted to its germline and aligned with the chain. Other elements (P, N) are assumed not to have mutated from the germline. This hypothesis should be checked by examining other members of the clonal family.

## Usage ##

    python RevertToGermlineCmd.py imgt_nt germline_file species_name outfile option [-mdd]

Argument|Meaning
--------|-------
imgt_nt|Filename of an IMGT nucleotide analysis csv or tab separated file containing analyses of the sequences
germline_file|Filename of the germline library to use (single file in IMGT format). The IMGT germline library for all species is downloadable [here](http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP)
species_name|The species name as used in the germline file, e.g. "Homo sapiens"
outfile|The name of the output file.
option|Any combination of:
| |f - list 'full' germlines where genes are reverted and ns preserved
| |v - list germlines where just the v-gene is reverted and other regions gapped
| |j - list germlines where all genes are reverted and n-regions gapped
| |i - list input sequences
| |o - list derived germline for each input sequence
| |c - list consensus for derived germlines
| |x - write verbose analysis through the report functio 
|-m*dd*|-m followed immediately by a two digit number will list V-gene mutations from germline that are seen in all sequences of a particular germline, for all germlines where there are at least that number of sequences.

Output sequences are aligned on full codon boundaries and incomplete codons are gapped out

## Testing ##

Executing the commands in test/test.sh will exercise AnnotateTree and RevertToGermline with sample data.