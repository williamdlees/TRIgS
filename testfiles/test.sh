# Before running this test script, please download the IMGT germline library to testfiles/imgt_germlines.fasta.
# The library can be found at http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP
python ../RevertToGermlineCmd.py Cluster265_IMGT_Nt.csv imgt_germlines.fasta "Danio rerio" test_Cluster265_germlines_germ_of.fasta  of
python ../RevertToGermlineCmd.py Cluster265_IMGT_Nt.csv imgt_germlines.fasta "Danio rerio" test_Cluster265_germlines_germ_ocf.fasta  ocf
python ../RevertToGermlineCmd.py Cluster265_IMGT_Nt.csv imgt_germlines.fasta "Danio rerio" test_Cluster265_germlines_germ_oicf.fasta  oicf
python ../RevertToGermlineCmd.py Cluster265_IMGT_Nt.csv imgt_germlines.fasta "Danio rerio" test_Cluster265_germlines_cj.fasta  cj
python ../AnnotateTreeCmd.py seqnumfile.txt Cluster265_germlines_g_j.fasta Cluster265_germlines_g_j.fasta.treefile cdrfile.txt C265 C265
