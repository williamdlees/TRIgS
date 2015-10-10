# BioTools
A collection of bioinformatics tools for sequence analysis, with an emphasis on Next-Generation Sequencing and [Rep-Seq](http://www.ncbi.nlm.nih.gov/pubmed/22043864 "Rep-Seq").


[**AnnotateTreeCmd**](docs/AnnotateTree.md) creates annotated lineage trees and sequence alignments showing the point at which amino acid substitutions occur. It uses [PHYLIP’s dnaml](http://evolution.genetics.washington.edu/phylip.html) for ancestral reconstruction. Sequence numbering can be defined by the user, for example to match the numbering of a crystal structure, or to match a standard numbering scheme. If the sequences represent a B-cell clonal lineage, additional reports relating to variation in the CDRs can be produced. An [online version](http://cimm.ismb.lon.ac.uk/pat/annotatetree/) is available.

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

[**RevertToGermlineCmd**](docs/RevertToGermline.md) uses a simple approach to infer the germline ancestor of a B-cell variable region sequence, given the IMGT junction analysis. If a clonal lineage is available, the inferred germline can be used to root a phylogenetic tree, from which a more accurate germline can then be inferred using AnnotateTree. An [online version](http://cimm.ismb.lon.ac.uk/pat/germline/) is available.

[**IgBLASTPlus**](docs/IgBLASTPlus.md) processes the output from NCBI's [IgBLAST]("http://www.ncbi.nlm.nih.gov/igblast/"),
providing a full junction analysis and summarising results in an [IMGT]("http://imgt.org")-style tab-separated format. 

[**ExtractFromIMGT**](docs/ExtractFromIMGT.md) is a flexible tool for extracting sequences from IMGT/IgBLASTPlus files.

## Further Information ##

William D. Lees and Adrian J. Shepherd, [&#8220;Utilities for High-Throughput Analysis of B-Cell Clonal Lineages,&#8221;](
http://www.hindawi.com/journals/jir/2015/323506/) <i>Journal of Immunology Research</i>, vol. 2015, Article ID 323506, 9 pages, 2015. doi:10.1155/2015/323506 <li><a href="http://files.hindawi.com/journals/jir/2015/323506.enw">Download citation as EndNote</a>

## Contact ##

william@lees.org.uk
