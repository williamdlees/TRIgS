# TRIgS - Tools for Rendering Ig Sequences
A collection of bioinformatics tools for sequence analysis, with an emphasis on Next-Generation Sequencing and [Rep-Seq](http://www.ncbi.nlm.nih.gov/pubmed/22043864).

Online versions of some tools are available on our [website](http://cimm.ismb.lon.ac.uk/pat).

[Tools for Clonal Analysis](#tools-for-clonal-analysis)<br>
[Tools for Junction Parsing and Results Manipulation](#tools-for-junction-parsing-and-results-manipulation)<br>
[Tools for FASTA file Manipulation](#tools-for-fasta-file-manipulation)<br>

[This document](docs/Example%20Analysis%20Pipeline.pdf) illustrates the use of Trigs in combination with other tools in a recent analysis.

## Tools for Clonal Analysis

[**ClusterSeqs**](docs/ClusterSeqs.md) partitions sequences into clusters using single-linkage clustering. In a recent analysis, 426,000 junction sequences were clustered in just under 90 minutes. [**NeighbourDist**](docs/NeighbourDist.md) analyses nearest-neighbour distances to guide partitioning, and can down-sample to handle large datasets.

[**ClusterGraph**](docs/ClusterGraph.md) converts output from ClusterSeqs or [CD-HIT](http://weizhongli-lab.org/cd-hit/) into a form that can be imported by [Gephi](https://gephi.org/).

<img src="https://rawgit.com/williamdlees/BioTools/master/docs/clusters.png" width="200">
<br>Example Gephi output, produced from the analysis of an NGS-based heavy-chain repertoire

[**AnnotateTreeCmd**](docs/AnnotateTree.md) creates annotated lineage trees and sequence alignments showing the point at which amino acid substitutions occur. It uses [PHYLIP’s dnaml](http://evolution.genetics.washington.edu/phylip.html) for ancestral reconstruction. Sequence numbering can be defined by the user, for example to match the numbering of a crystal structure, or to match a standard numbering scheme. If the sequences represent a B-cell clonal lineage, additional reports relating to variation in the CDRs can be produced.

![Image](https://rawgit.com/williamdlees/BioTools/master/docs/treediag5.svg)
<br>Part of a tree produced by AnnotateTreeCmd, showing amino acid substitutions

                                                             1         1         1     
                                7        8         9         0         1         2     
                       6784567890124567890123456789012345678901234567890123456789012345
	consensus_germ_vdj YYSDSDKSTAQSVQGRFTASKDSSNLYLHMNQLKTEDSAVYYCA-EW-GAFDYWGKGTMVTVTS
	216                SVG...N..................F...........T......R.RY................
	605                SVG......................F...........T......R.RY..........NGHCHI
	742                PVG......................F...........T......R.ID................
	154                SVG......................F...........T......R.RY................
<br>Part of an alignment produced by AnnotateTreeCmd, showing custom numbering including deletions

[**RevertToGermlineCmd**](docs/RevertToGermline.md) uses a simple approach to infer the germline ancestor of a B-cell variable region sequence, given the IMGT junction analysis. If a clonal lineage is available, the inferred germline can be used to root a phylogenetic tree, from which a more accurate germline can then be inferred using AnnotateTree.

## Tools for Junction Parsing and Results Manipulation

[**IgBLASTPlus**](docs/IgBLASTPlus.md) processes the output from NCBI's [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/), providing a full junction analysis and summarising results in an [IMGT](http://imgt.org)-style tab-separated format. This makes it possible to use an in-house copy of IgBLAST in place of IMGT High-V Quest. IgBLAST analyses of 3-4 million records will typically complete in under an hour. 

The following tools will work on tab- or comma- separated files such as those produced by [IgBLASTPlus](docs/IgBLASTPlus.md), [IMGT](http://imgt.org), or [CHANGE-O](http://clip.med.yale.edu/changeo):

[**ExtractFromIMGT**](docs/ExtractFromIMGT.md) is a flexible tool for extracting sequences in FASTA format. Options allow filtering by germlines, restriction of sequences to specific regions, and so on.

[**PlotGermline**](docs/PlotGermline.md) creates histograms showing germline usage. 

[**Spectratype**](docs/Spectratype.md) creates histograms showing CDR3 length distribution.

[**ClusterExtract**](docs/ClusterExtract.md) uses the sequence IDs in a [**ClusterSeqs**](docs/ClusterSeqs.md) output file to extract all records corresponding to a nominated cluster.

## Tools for FASTA file Manipulation

A collection of small utilities that have proved useful in analysis pipelines. 

[**CountRecords**](docs/FastaTools.md/#countrecords) counts the number of records in a file, accounting for duplicates noted in the header.

[**FastaMatch**](docs/FastaTools.md/#fastamatch) filters records whose sequence or ID match a regular expression.

[**FastaSample**](docs/FastaTools.md/#fastasample) extracts a random sample of records.

[**FastaUniq**](docs/FastaTools.md/#fastauniq) removes duplicate records (identified by ID).

[**FastaSampleUniq**](docs/FastaTools.md/#fastasampleuniq) counts the number of unique sequences in a set of samples, accounting for duplicates noted in the header.

## Installation and Usage

To use the tools, clone this repository or download and unzip the Zip file. All tools require <a href="https://www.python.org/downloads">Python 2 (v 2.7 or later)</a>, and <a href="http://biopython.org">BioPython.</a>. Additional dependencies and usage information are given in the links above for each tool.         

## Further Information

William D. Lees and Adrian J. Shepherd, [&#8220;Utilities for High-Throughput Analysis of B-Cell Clonal Lineages,&#8221;](
http://www.hindawi.com/journals/jir/2015/323506/) <i>Journal of Immunology Research</i>, vol. 2015, Article ID 323506, 9 pages, 2015. doi:10.1155/2015/323506 <li><a href="http://files.hindawi.com/journals/jir/2015/323506.enw">Download citation as EndNote</a>

## Contact

william@lees.org.uk
