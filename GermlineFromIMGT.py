# Copyright (c) 2015 William Lees

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Function to infer the germline sequences of a set of chains, given their IMGT nucleotide analysis.
# Each gene is reverted to its germline and aligned with the chain. Other elements (P, N) are assumed not to have
# mutated from the germline. This hypothesis should be checked by examining other members of the clonal
# family.
#
# sequence_file is the name of an IMGT nucleotide analysis csv or tab separated file containing analyses of the sequences
# germline_lib is the name of the germline library to use (single file in IMGT format). The IMGT germline library
# for all species is downloadable at http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP
# species_name is the species name as used in the germline file, e.g. "Homo sapiens"
# output_file is the name of the output file.
# option is any combination of:
# f - list 'full' germlines where genes are reverted and ns preserved
# v - list germlines where just the v-gene is reverted and other regions gapped
# j - list germlines where all genes are reverted and n-regions gapped
# i - list input sequences
# o - list derived germline for each input sequence
# c - list consensus for derived germlines
# x - write verbose analysis through the report functio 
# fixed_mut - if this parameter is nonzero, the function will list V-gene mutations from germline that are seen in all
# sequences of a particular germline, for all germlines where there are at least fixed_mut sequences.
# in all cases, output sequences are aligned on full codon boundaries and incomplete codons are gapped out
# report is a function called whenever there is status info to report
# error is a function called when there is a (catastrophic) error

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import sys
import csv
import traceback
from Germlib import Germlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment

def germline_from_imgt(sequence_file, germline_lib, species_name, output_file, option, report, error, fixed_mut):
    """
    Test cases mainly for Mutation Analysis:
    Identical sequences with varying lengths at the 5' end
    >>> germline_from_imgt("Testfiles/alleles/varyinglengths.txt", "Testfiles/alleles/imgt_germlines.fasta", "Homo sapiens", "Testfiles/alleles/varyinglengths_out.fasta", "v", doctest_report, doctest_report, 3)
    Processing successfully completed.
    Mutation Analysis:
    IGHV4-34*01 (5 sequences):
     Common mutations: c67g, g88c, t96g, a97c, g103c, t158a, a165g, g166c, c179t, c180g, c181t, g182t, a189g, t207c, c208a, a209c, g210a, a212t, c221t, a229g, g230a, a241g, c248t, t249g, g273a, t274a, g275a, t278c, a280t
    germline:  caggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgctgtctatggtgggtccttcagtggttactactggagctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcatagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcagtagacacgtccaagaaccagttctccctgaagctgagctctgtgaccgccgcggacacggctgtgtattactgtgcgagagg
    consensus: .........ctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgGtgtctatggtgggtccttcaCtggttacGCctggaCctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcaAagtggaGCcaccaactacaaTGTTtccctcGagagtcgagtcaccataCACAtTgacacgtcTaagaaccGAttctccctgaGgctgagTGctgtgaccgccgcggacacggctAAAtaCtTctgtgcgagagg

    Identical sequences with one common deletion and an adjacent deletion that is not present in all sequences
    >>> germline_from_imgt("Testfiles/alleles/somedeletions.txt", "Testfiles/alleles/imgt_germlines.fasta", "Homo sapiens", "Testfiles/alleles/somedeletions_out.fasta", "v", doctest_report, doctest_report, 3)
    Processing successfully completed.
    Mutation Analysis:
    IGHV4-34*01 (5 sequences):
     Common deletions: 156, 157, 158
     Common mutations: c67g, g88c, t96g, a97c, g103c, a165g, g166c, c179t, c180g, c181t, g182t, a189g, t207c, c208a, a209c, g210a, a212t, c221t, a229g, g230a, a241g, c248t, t249g, g273a, t274a, g275a, t278c, a280t
    germline:  caggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgctgtctatggtgggtccttcagtggttactactggagctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcatagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcagtagacacgtccaagaaccagttctccctgaagctgagctctgtgaccgccgcggacacggctgtgtattactgtgcgagagg
    consensus: ...gtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgGtgtctatggtgggtccttcaCtggttacGCctggaCctggatccgccagcccccagggaaggggctggagtggattggggaaatcXXX---agtggaGCcaccaactacaaTGTTtccctcGagagtcgagtcaccataCACAtTgacacgtcTaagaaccGAttctccctgaGgctgagTGctgtgaccgccgcggacacggctAAAtaCtTctgtgcgagagg

    Identical sequences with a common deletion
    >>> germline_from_imgt("Testfiles/alleles/deletions.txt", "Testfiles/alleles/imgt_germlines.fasta", "Homo sapiens", "Testfiles/alleles/deletions_out.fasta", "v", doctest_report, doctest_report, 3)
    Processing successfully completed.
    Mutation Analysis:
    IGHV4-34*01 (5 sequences):
     Common deletions: 153, 154, 155, 156, 157, 158
     Common mutations: c67g, g88c, t96g, a97c, g103c, a165g, g166c, c179t, c180g, c181t, g182t, a189g, t207c, c208a, a209c, g210a, a212t, c221t, a229g, g230a, a241g, c248t, t249g, g273a, t274a, g275a, t278c, a280t
    germline:  caggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgctgtctatggtgggtccttcagtggttactactggagctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcatagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcagtagacacgtccaagaaccagttctccctgaagctgagctctgtgaccgccgcggacacggctgtgtattactgtgcgagagg
    consensus: ...gtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgGtgtctatggtgggtccttcaCtggttacGCctggaCctggatccgccagcccccagggaaggggctggagtggattggggaaatc------agtggaGCcaccaactacaaTGTTtccctcGagagtcgagtcaccataCACAtTgacacgtcTaagaaccGAttctccctgaGgctgagTGctgtgaccgccgcggacacggctAAAtaCtTctgtgcgagagg

    Test with 5 identical sequences
    >>> germline_from_imgt("Testfiles/alleles/identical.txt", "Testfiles/alleles/imgt_germlines.fasta", "Homo sapiens", "Testfiles/alleles/identical_out.fasta", "v", doctest_report, doctest_report, 3)
    Processing successfully completed.
    Mutation Analysis:
    IGHV4-34*01 (5 sequences):
     Common mutations: c67g, g88c, t96g, a97c, g103c, t158a, a165g, g166c, c179t, c180g, c181t, g182t, a189g, t207c, c208a, a209c, g210a, a212t, c221t, a229g, g230a, a241g, c248t, t249g, g273a, t274a, g275a, t278c, a280t
    germline:  caggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgctgtctatggtgggtccttcagtggttactactggagctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcatagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcagtagacacgtccaagaaccagttctccctgaagctgagctctgtgaccgccgcggacacggctgtgtattactgtgcgagagg
    consensus: ...gtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgGtgtctatggtgggtccttcaCtggttacGCctggaCctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcaAagtggaGCcaccaactacaaTGTTtccctcGagagtcgagtcaccataCACAtTgacacgtcTaagaaccGAttctccctgaGgctgagTGctgtgaccgccgcggacacggctAAAtaCtTctgtgcgagagg

    Test with two mutations removed from the first sequence
    >>> germline_from_imgt("Testfiles/alleles/two_dropped.txt", "Testfiles/alleles/imgt_germlines.fasta", "Homo sapiens", "Testfiles/alleles/two_dropped_out.fasta", "v", doctest_report, doctest_report, 3)
    Processing successfully completed.
    Mutation Analysis:
    IGHV4-34*01 (5 sequences):
     Common mutations: t96g, a97c, g103c, t158a, a165g, g166c, c179t, c180g, c181t, g182t, a189g, t207c, c208a, a209c, g210a, a212t, c221t, a229g, g230a, a241g, c248t, t249g, g273a, t274a, g275a, t278c, a280t
    germline:  caggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgctgtctatggtgggtccttcagtggttactactggagctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcatagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcagtagacacgtccaagaaccagttctccctgaagctgagctctgtgaccgccgcggacacggctgtgtattactgtgcgagagg
    consensus: ...gtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgtccctcacctgcgXtgtctatggtgggtccttcaXtggttacGCctggaCctggatccgccagcccccagggaaggggctggagtggattggggaaatcaatcaAagtggaGCcaccaactacaaTGTTtccctcGagagtcgagtcaccataCACAtTgacacgtcTaagaaccGAttctccctgaGgctgagTGctgtgaccgccgcggacacggctAAAtaCtTctgtgcgagagg
    """
    for char in option:
        if char not in 'ciofvjx':
            error('unrecognised option: %s.' % char)
            return
        
    try:
        gl = Germlib(species_name, germline_file=germline_lib)
    except:
        report("Error parsing germline library file: " + str(sys.exc_info()[1]))
        return

    consensus_f = []
    consensus_v = []
    consensus_j = []
    mutated_germs = {}
    
    imgt_nt = {}
    try:
        with open(sequence_file, "r") as sequence_handle:
            ln = sequence_handle.readline()
            sep = ("\t" if "\t" in ln else ",")
            sequence_handle.seek(0)
            reader = csv.DictReader(sequence_handle, delimiter=sep)
            for row in reader:
                imgt_nt[row["Sequence ID"]] = row
    
        outrecs = []
        for id, nt_rec in imgt_nt.iteritems():
            try:
                if "JUNCTION" in nt_rec and nt_rec["JUNCTION"] != None and len(nt_rec["JUNCTION"]) > 0:
                    heavychain = len(nt_rec["V-D-J-REGION"]) > 0
        
                    if heavychain:
                        mAb = (nt_rec["V-REGION"], 
                               nt_rec.get("P3'V", ""), 
                               nt_rec.get("N-REGION", ""), 
                               nt_rec.get("N1-REGION", ""), 
                               nt_rec.get("P5'D", ""),
                               nt_rec.get("D-REGION", ""), 
                               nt_rec.get("P3'D", ""), 
                               nt_rec.get("N2-REGION", ""), 
                               nt_rec.get("P5'J", ""), 
                               nt_rec["J-REGION"])
                    else:
                        mAb = (nt_rec["V-REGION"], 
                               nt_rec.get("P3'V", ""), 
                               nt_rec.get("N-REGION", ""), 
                               nt_rec.get("P5'J", ""), 
                               nt_rec["J-REGION"])
        
                    if 'x' in option:
                        report("%s:" % id)
                        report(" | ".join(mAb))
        
                    # Revert the part of the V-gene that extends to the second Cysteine
                    vregion = nt_rec["V-REGION"]
                    vregion_3prime = nt_rec["3'V-REGION"]
                    vgene_name = Germlib.translate_imgt_name(nt_rec["V-GENE and allele"])
        
                    if len(vregion_3prime) > 0 and vregion[0 - len(vregion_3prime):] != vregion_3prime:
                        report("Error: 3'V-REGION sequence not found at 3' end of V-REGION in sequence %s" % id)
                        continue
        
                    # Remove stray nucleotides from the 5' end of the V-region to give us whole codons (we know the 3' end is aligned)
                    vregion_5prime = vregion[:0 - len(vregion_3prime)] if len(vregion_3prime) > 0 else vregion
                    vregion_5prime = (vregion_5prime if len(vregion_5prime) % 3 == 0 else vregion_5prime[(len(vregion_5prime) % 3):])
        
                    try:
                        vgene_frag1, matchstr_frag1 = gl.match_from_aa(vgene_name, vregion_5prime)
        
                        # For the remaining (3') part, use a global alignment. We use the entire V-region so that the 3prime
                        # region, which might be quite small, aligns against the right part of the sequence
                        vgene_frag2, matchstr_frag2 = gl.match(vgene_name, vregion)[0 - len(vregion_3prime):] if len(vregion_3prime) > 0 else ("", "")
                        
                        if fixed_mut > 0:
                            # Merge the two matchstrings. Starting at the 3' end, we pull matchstring off frag2 until we get beyond vgene_frag2 and are
                            # about to pull the first nt of vgene_frag1. Then we pull the rest off vgene_frag1.
                            
                            mlen = 0
                            matchstr = ""
                            for m in matchstr_frag2[::-1]:
                                if m != 'd':
                                    mlen += 1
                                    if mlen > len(vregion_3prime):
                                        break
                                matchstr += m
                                
                            skip = True
                            for m in matchstr_frag1[::-1]:
                                if skip and m != 'd':
                                    skip = False
                                if not skip:
                                    matchstr += m
                            
                            matchstr = matchstr[::-1]
                            
                            # Sanity check 1 - number of nucleotides in match string should match length of v-region
                            mlen = sum((n != 'd') for n in matchstr)
                            if len(vregion_5prime) + len(vregion_3prime) != mlen:
                                report("Error in match string length for sequence %s" % id)
                                
                            # Sanity check 2 - check matchstring is consistent
                            
                            vgene = str(gl.seq(vgene_name).seq)
                            mismatch = False
                            gt = iter(vgene)
                            vt = iter(vregion)
                            for m in matchstr:
                                if m == 'd':
                                    next(gt)
                                elif m == 'i':
                                    next(vt)
                                elif m == 'm':
                                    if next(gt) != next(vt):
                                        mismatch = True
                                else:
                                    if next(gt) == next(vt):
                                        mismatch = True
                                        
                            if mismatch:
                                report("Error in matchstring for sequence %s:\nvgene: %s\nseq  :  %s\nmatch: %s\n" % (id, vgene, vregion, matchstr))
                            else:
                                en = mutated_germs.get(vgene_name, [])
                                en.append((vregion, matchstr))
                                mutated_germs[vgene_name] = en
                                
                        if nt_rec["J-GENE and allele"] != '':
                            jgene_name = Germlib.translate_imgt_name(nt_rec["J-GENE and allele"])
                            jgene_frag, _ = gl.match(jgene_name, nt_rec["J-REGION"])
                        else:
                            jgene_frag = ''
        
                        if heavychain and nt_rec["D-GENE and allele"] != '':
                            dgene_name = Germlib.translate_imgt_name(nt_rec["D-GENE and allele"])
                            dgene_frag, _ = gl.match(dgene_name, nt_rec["D-REGION"])
                        else:
                            dgene_frag = ''
                    except:
                        report("Error processing sequence " + id + ":")        
                        exc_type, exc_value, exc_traceback = sys.exc_info()
                        report(traceback.format_exception(exc_type, exc_value, exc_traceback, 2))
                        continue
        
                    if heavychain:
                        germline = [
                            vgene_frag1 + vgene_frag2, 
                            nt_rec.get("P3'V", ""), 
                            nt_rec.get("N-REGION", ""), 
                            nt_rec.get("N1-REGION", ""), 
                            nt_rec.get("P5'D", ""),
                            dgene_frag, 
                            nt_rec.get("P3'D", ""), 
                            nt_rec.get("N2-REGION", ""), 
                            nt_rec.get("P5'J", ""), 
                            jgene_frag]
                    else:
                        germline = [
                            vgene_frag1 + vgene_frag2, 
                            nt_rec.get("P3'V", ""), 
                            nt_rec.get("N-REGION", ""), 
                            nt_rec.get("P5'J", ""), 
                            jgene_frag]
        
                    jgene_frag = (jgene_frag if len("".join(germline)) % 3 == 0 else jgene_frag[:0-(len("".join(germline)) % 3)])
                    germline[-1] = jgene_frag
        
                    if 'i' in option:                
                        trunc5 = len(vregion) - len(vregion_5prime + vregion_3prime)
                        if heavychain:
                            trunc3 = (len(nt_rec["V-D-J-REGION"]) - trunc5) % 3
                            if trunc3 != 0:
                                outrecs.append(SeqRecord(Seq(nt_rec["V-D-J-REGION"][trunc5:0-trunc3]), id=id, name=id, description=""))
                            else:
                                outrecs.append(SeqRecord(Seq(nt_rec["V-D-J-REGION"][trunc5:]), id=id, name=id, description=""))
                        else:
                            trunc3 = (len(nt_rec["V-J-REGION"]) - trunc5) % 3
                            if trunc3 != 0:
                                outrecs.append(SeqRecord(Seq(nt_rec["V-J-REGION"][trunc5:0-trunc3]), id=id, name=id, description=""))
                            else:
                                outrecs.append(SeqRecord(Seq(nt_rec["V-J-REGION"][trunc5:]), id=id, name=id, description=""))
        
                    if 'f' in option:
                        if 'x' in option:
                            report("Inferred 'full' germline:")
                            report(" | ".join(germline))
                        sr = SeqRecord(Seq("".join(germline)), id=id + "_germ", name=id + "_germ", description="")
                        consensus_f.append(sr)
                        if 'o' in option:
                            outrecs.append(sr)
        
                    def chunks(l, n):
                        """ Yield successive n-sized chunks from l."""
                        for i in xrange(0, len(l), n):
                            yield l[i:i + n]
        
                    germline = "".join(germline)
                    v_ext = vgene_frag1 + vgene_frag2
        
                    if 'v' in option:
                        g = (v_ext) + '-' * (len(germline) - len(v_ext))
                        germline_v = ""
                        for c in chunks(g, 3):
                            germline_v += c if '-' not in c else '-'*len(c)
            
                        if 'x' in option:
                            report("Inferred germline (v):")
                            report(germline_v)
                        sr = SeqRecord(Seq(germline_v), id=id + "_germ_v", name=id + "_germ_v", description="")
                        consensus_v.append(sr)
                        if 'o' in option:
                            outrecs.append(sr)
        
                    if 'j' in option:
                        if heavychain:
                            g = v_ext + '-' * (
                                len(nt_rec.get("P3'V", "")) + 
                                len(nt_rec.get("N-REGION", "")) + 
                                len(nt_rec.get("N1-REGION", "")) + 
                                len(nt_rec.get("P5'D", ""))) + \
                                dgene_frag + \
                                '-' * (
                                len(nt_rec.get("P3'D", "")) + 
                                len(nt_rec.get("N2-REGION", "")) + 
                                len(nt_rec.get("P5'J", ""))) + \
                                jgene_frag
                        else:
                            g = v_ext + '-' * (len(germline) - len(v_ext) - len(jgene_frag)) + jgene_frag
                            
                        germline_vj = ""
                        for c in chunks(g, 3):
                            germline_vj += c if '-' not in c else '-'*len(c)
            
                        if 'x' in option:
                            report("Inferred germline_vdj:")
                            report(germline_vj)
                        sr = SeqRecord(Seq(germline_vj), id=id + "_germ_vdj", name=id + "_germ_vdj", description="")
                        consensus_j.append(sr)
                        if 'o' in option:
                            outrecs.append(sr)
                else:
                    report("%s: no junction." % id)
            except:
                report("Error processing input record " + id + ":")        
                exc_type, exc_value, exc_traceback = sys.exc_info()
                report(traceback.format_exception(exc_type, exc_value, exc_traceback, 2))
                
        report("Processing successfully completed.")
        
    except:
        report("Error parsing input file: " + str(sys.exc_info()[1]))
        return

                
    if 'c' in option:
        try:
            def checklengths(srs):
                length = -1
                for sr in srs:
                    if length < 0:
                        length = len(sr.seq)
                    elif len(sr.seq) != length:
                        report("Length error in sequence %s" % sr.id)

            if 'f' in option:
                checklengths(consensus_f)
                summary = AlignInfo.SummaryInfo(MultipleSeqAlignment(consensus_f))
                cd = summary.dumb_consensus(ambiguous="-")
                consensus = ""
                for c in chunks(cd, 3):
                    consensus += c if '-' not in c else '-'*len(c)
                report("'Full' germline consensus:")
                report(str(consensus))
                outrecs.insert(0, SeqRecord(consensus, id="consensus_germ_full", name="consensus_germ_full", description=""))
            if 'v' in option:
                checklengths(consensus_v)
                summary = AlignInfo.SummaryInfo(MultipleSeqAlignment(consensus_v))
                cd = summary.dumb_consensus(ambiguous="-")
                consensus = ""
                for c in chunks(cd, 3):
                    consensus += c if '-' not in c else '-'*len(c)
                report("Germline (v) consensus:")
                report(str(consensus))
                outrecs.insert(0, SeqRecord(consensus, id="consensus_germ_v", name="consensus_germ_v", description=""))
            if 'j' in option:
                checklengths(consensus_j)
                summary = AlignInfo.SummaryInfo(MultipleSeqAlignment(consensus_j))
                cd = summary.dumb_consensus(ambiguous="-")
                consensus = ""
                for c in chunks(cd, 3):
                    consensus += c if '-' not in c else '-'*len(c)
                report("Germline vdj consensus:")
                report(str(consensus))
                outrecs.insert(0, SeqRecord(consensus, id="consensus_germ_vdj", name="consensus_germ_vdj", description=""))
        except:
            report("Error generating consensus: %s - %s" % (sys.exc_info()[0], sys.exc_info()[1]))
            
    if fixed_mut > 0:
        try:
            report("Mutation Analysis, showing mutations, insertions and deletions that are common to all sequences from a given germline.")
            report("This will be reported for all germlines for which there are at least %d sequences in the analysis:" % fixed_mut)
            def m_limits(m):
                # Find the upper and lower limits of the matchstr, ignoring leading and trailing deletions
                # limits are expressed as locations relative to the germline (insertions in the matchstr are ignored)
                for i in range(len(m)):
                    if m[i] != 'd':
                        mstart = i
                        break
                for i in range(len(m)-1, -1, -1):
                    if m[i] != 'd' and m[i] != 'i':
                        mend = i
                        break
                        
                loc = 0
                for i in range(len(m)):
                    if i == mstart:
                        start = loc
                    elif i == mend:
                        end = loc
                    if m[i] != 'i':
                        loc += 1
                return (start, end)
            
            for germline, mg in mutated_germs.iteritems():
                if len(mg) >= fixed_mut:
                    # given that the sequences may have different start and end points, compute
                    # the range over which we have coverage from a sufficient number of sequences
                    germseq = gl.seq(germline).seq
                    coverage = [0] * len(germseq)
                    for seq, matchstr in mg:
                        start, end = m_limits(matchstr)
                        for i in range(start, end+1):
                            coverage[i] += 1
                            
                    range_start = 999
                    range_end = -1
                    
                    for i, val in enumerate(coverage):
                        if val >= fixed_mut: 
                            if range_start > i:
                                range_start = i
                            if range_end < i:
                                range_end = i
                        
                    # matches[loc] holds:
                    # 'u' if this location has not as yet been observed in sequences processed
                    # 'm' if it has been observed to match the germline in sequences processed so far
                    # 'c,g,a,t' if it has been observed to be mutated to that value in sequences processed so far
                    # 'm' if it has been observed to be deleted in sequences processed so far
                    # 'x' if if the results at this location are not consistent between sequences
                    
                    matches = ['u'] * len(germseq)
                    insertions = []
                    range_encountered_start = 999
                    range_encountered_end = -1
                    
                    for seq, matchstr in mg:
                        ins = 0
                        loc = 0
                        inserts = []
                        (start, end) = m_limits(matchstr)
                        start = max(start, range_start)
                        end = min(end, range_end)
                        s = iter(seq)
                        for m in matchstr:
                            if m != 'i':
                                ins = 0
                                
                            if m == 'n':
                                sub = next(s)
                                if loc >= start and loc <= end:
                                    if matches[loc] == 'u':
                                        matches[loc] = sub
                                    elif matches[loc] != sub:
                                        matches[loc] = 'x'
                                loc += 1
                            elif m == 'd':
                                if loc >= start and loc <= end:
                                    if matches[loc] == 'u':
                                        matches[loc] = 'd'
                                    elif matches[loc] != 'd':
                                        matches[loc] = 'x'
                                loc += 1
                            elif m == 'i':
                                if loc >= start and loc <= end:
                                    inserts.append((loc, ins))
                                ins += 1
                                next(s)
                            else:
                                if loc >= start and loc <= end:
                                    if matches[loc] == 'u':
                                        matches[loc] = 'm'
                                    elif matches[loc] != 'm':
                                        matches[loc] = 'x'
                                loc += 1
                                next(s)
                                        
                        # Add a new insertion to the consensus list if we see it in this sequence, and it is outside
                        # the range we've encountered so far. 
                        
                        for loc, ins in inserts:
                            if loc < range_encountered_start or loc > range_encountered_end:
                                insertions.append((loc, ins))
                        
                        # Remove insertions from the consensus list if they are in range of this sequence and were not
                        # observed in it
                        
                        for loc, ins in insertions:
                            if loc >= start and loc <= end:
                                if not (loc, ins) in inserts:
                                    insertions.remove((loc, ins))
                        
                        range_encountered_start = min(range_encountered_start, start)        
                        range_encountered_end = max(range_encountered_end, end)        
                
                    report("%s (%d sequences):" % (germline, len(mg)))
                    deletions = []
                    for loc, m in enumerate(matches):
                        if m == 'd':
                            deletions.append(loc)
                    if len(deletions) > 0:
                        report(" Common deletions: %s" % ', '.join([str(n) for n in sorted(deletions)]))
                    
                    if len(insertions) > 0:
                        report(" Common insertions: %s") % ', '.join(["%d.%d" % (loc, ins) for (loc, ins) in sorted(insertions)])
    
                    mutations = []
                    for loc, m in enumerate(matches):
                        if m in ('c', 'a', 'g', 't'):
                            mutations.append("%s%d%s" % (germseq[loc], loc, m))
                    if len(mutations) > 0:
                         report(" Common mutations: %s" % ', '.join([str(n) for n in mutations]))
                        
                    if len(insertions) + len(deletions) + len(mutations) > 0:
                        r_g = ""
                        gi = iter(germseq)
                        for m in matches:
                            r_g += next(gi) if m != 'i' else '-'
                        report( "germline:  %s" % r_g)
                        r_c = ""
                        gi = iter(germseq)
                        for m in matches:
                            if m == 'm':
                                r_c  += next(gi)
                            elif m == 'd':
                                r_c += '-'
                                next(gi)
                            elif m == 'i':
                                r_c += 'i'
                            elif m == 'u':
                                r_c += '.'
                                next(gi)
                            else:
                                r_c += m.upper()
                                next(gi)
                        report( "consensus: %s" % r_c)
                    else:
                        report(" No common insertions, deletions or mutations compared to gertmline")
                else:                            
                    report("%s (%d sequences) - number of sequences is below analysis threshold." % (germline, len(mg)))
        except:
            report("Error creating mutation report:")        
            exc_type, exc_value, exc_traceback = sys.exc_info()
            report(traceback.format_exception(exc_type, exc_value, exc_traceback, 2))
            
    SeqIO.write(outrecs, output_file, "fasta")

def doctest_report(x):
    print(x)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

