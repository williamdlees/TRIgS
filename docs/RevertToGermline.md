# RevertToGermlineCmd #

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
| |One or more *Output controls*, defining what should be included in the output file:
| |o - include a derived germline corresponding to each input sequence
| |c - include a consensus sequence built from all the derived germlines
| |i - include the input sequences themselves
| |One or more *Content controls*, defining how the germlines should be derived. The resulting sequences include a suffix in their name (shown in brackets below) so that they can be distinguished:
| |f - 'full' germlines where all genes are reverted and n-regions preserved (suffix _germ)
| |v - germlines where just the v-gene is reverted and other regions gapped (suffix _germ_v)
| |j - germlines where all genes are reverted and n-regions gapped (suffix _germ_vdj)
| |*Diagnostic controls:* 
| |x - write verbose analysis through the report function
|-m*dd*|-m followed immediately by a two digit number will list V-gene mutations from germline that are seen in all sequences of a particular germline, for all germlines where there are at least that number of sequences.

Output sequences are aligned on full codon boundaries and incomplete codons are gapped out

## Example Options ##

Option|Resulting Output Sequences
------|--------------------------
| of  |A reverted sequence corresponding to each input sequence. All germlines are reverted, and n-regions preserved.
| ocf |As above, plus a consensus sequence derived from all reverted sequences.
| oicf |As above, and also all input sequences are included.
| cj  |A consensus sequence, derived from reverted sequences where all genes are reverted, and n-regions gapped.



## Testing ##

Executing the commands in test/test.sh will exercise AnnotateTreeCmd and RevertToGermlineCmd with sample data.