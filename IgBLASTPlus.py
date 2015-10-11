# Copyright (c) 2015 William Lees

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the 'Software'), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# A script to condense the IgBlast output format (produced with the option -outfmt 3) into an IMGT-stype tab-separated file.
# Uses the same column headers as IMGT
# Determines the end of the junction by looking for the conserved W (VH) or F (VL) - required to be present at the same
# location in both the query sequence and the J-germline. 

import sys
import csv
import argparse
import re
import traceback
from Bio import SeqIO
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser(description='Summarise each IgBlast analysis on a single tab-separated line.',
                                     epilog='The program creates two output files: <tag>_n.fasta (nucleotide analysis) and <tag>_n.fasta (amino acid analysis)')
    parser.add_argument('seqfile', help='file containing the sequences analysed by IgBlast (FASTA)')
    parser.add_argument('jgermfile', help='the J-germline file used by IgBlast when analysing the sequences')
    parser.add_argument('igblastfile', help='analysis file from IgBlast (produced with option -outfmt 3)')
    parser.add_argument('tag', help='prefix for output files')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()
    
    verbose = args.verbose
    fieldorder = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']

    seqs = {}
    for seq_record in SeqIO.parse(args.seqfile, 'fasta'):
        seqs[seq_record.id] = str(seq_record.seq.upper())

    j_germs = {}
    for seq_record in SeqIO.parse(args.jgermfile, 'fasta'):
        seq_record.seq = seq_record.seq.upper()
        id = seq_record.id.split('|')[1] if len(seq_record.id.split('|')) > 1 else seq_record.id
        j_germs[id] = str(seq_record.seq)

    with open(args.igblastfile, 'r') as fi, open(args.tag + '_n.txt', 'wb') as fo_n, open(args.tag + '_aa.txt', 'wb') as fo_aa:
        fieldnames = ['Sequence ID', 'V-GENE and allele', 'D-GENE and allele', 'J-GENE and allele', 'Chain Type', 'Stop Codon', 'V-J frame', 'Functionality', 
                 'Strand', 'FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT', 'FR4-IMGT', 'JUNCTION-IMGT', 'V 3prime', 'V-D Junction', 'D-region', 'D-J Junction', 'V-J Junction', 'J 5prime', 'Notes', 'nt Sequence']
        writer_n = csv.DictWriter(fo_n, fieldnames, restval='', extrasaction='raise', dialect='excel-tab')
        writer_n.writeheader()
        writer_aa = csv.DictWriter(fo_aa, fieldnames, restval='', extrasaction='raise', dialect='excel-tab')
        writer_aa.writeheader()
        res = {}
        lastline = ''
        alignments = False
        
        for line in fi:
            try:
                rawline = line
                line = line.strip()
                items = line.split('\t')
                
                if 'Query= ' in line:
                    if len(res) > 0:
                        for k,v in res.iteritems():
                            if 'N/A' in v:
                                res[k] = ''         # blank out any N/As from IgBLAST for compatibility with IMGT

                        writer_n.writerow(res) 
                        res = translate_res(res)        
                        writer_aa.writerow(res)
                        
                    res = {}
                    res['Sequence ID'] = line.replace('Query= ', '').strip()
                    res['Notes'] = ''
                    alignments = False
                elif 'V-(D)-J rearrangement summary for query sequence' in lastline:               
                    res['V-GENE and allele'] = items[0]
                    
                    if 'VH' in line:
                        res['D-GENE and allele'] = items[1]
                        bump = 1
                    else:
                        res['D-GENE and allele'] = ''
                        bump = 0
                        
                    res['J-GENE and allele'] = items[bump+1]
                    res['Chain Type'] = items[bump+2]
                    res['Stop Codon'] = items[bump+3]
                    res['V-J frame'] = items[bump+4]                
                    res['Functionality'] = 'productive' if items[bump+5] == 'Yes' else 'unproductive'
                    res['Strand'] = items[bump+6]
    
                    if ',' in res['V-GENE and allele']:
                        selected_vgene = res['V-GENE and allele'].split(',')[0]
                    else:
                        selected_vgene = res['V-GENE and allele']
    
                    if ',' in res['D-GENE and allele']:
                        selected_dgene = res['D-GENE and allele'].split(',')[0]
                    else:
                        selected_dgene = res['D-GENE and allele']
    
                    if ',' in res['J-GENE and allele']:
                        selected_jgene = res['J-GENE and allele'].split(',')[0]
                    else:
                        selected_jgene = res['J-GENE and allele']
    
                elif 'be assigned to either rearranging gene' in lastline:
                    res['V 3prime'] = (items[0].replace('(', '').replace(')', ''))
                    if res['Chain Type'] == 'VH':
                        res['V-D Junction'] = (items[1].replace('(', '').replace(')', ''))
                        res['D-region'] = (items[2].replace('(', '').replace(')', ''))
                        res['D-J Junction'] = (items[3].replace('(', '').replace(')', ''))
                        res['J 5prime'] = (items[4].replace('(', '').replace(')', ''))
                    else:
                        res['V-J Junction'] = (items[1].replace('(', '').replace(')', ''))
                        res['J 5prime'] = (items[2].replace('(', '').replace(')', ''))
                        
                elif 'Alignments' in line:
                    alignments = True
                    al_narrative = ''
                    al_query_seq = ''
                    al_query_start = 0
                    al_query_frag = ''
                    al_query_frag_start = 0
                    al_query_j_start = 0
                    al_query_end = 0
                    al_vmatch_start = 0
                    al_jmatch_start = 0
                    al_jmatch_end = 0
                    raw_fr3 = ''
                    
                    if 'V-GENE and allele' in res:
                        vmatch = res['V-GENE and allele'] if len(res['V-GENE and allele']) == 1 else res['V-GENE and allele'].split(',')[0]
                    else:
                        vmatch = ''
                    if 'J-GENE and allele' in res:
                        jmatch = res['J-GENE and allele'] if len(res['J-GENE and allele']) == 1 else res['J-GENE and allele'].split(',')[0]
                    else:
                        jmatch = ''
                    if 'D-GENE and allele' in res:
                        dmatch = res['V-GENE and allele'] if len(res['V-GENE and allele']) == 1 else res['V-GENE and allele'].split(',')[0]
                    else:
                        dmatch = ''
                        
                elif alignments:
                    items = rawline.split()
                    if len(items) > 1:
                        if 'Query_' in line and len(items) == 4:
                            if len(al_query_seq) == 0:
                                al_query_start = int(items[1])
                            al_query_frag_start = int(items[1])
                            al_query_frag = items[2]
                            al_query_seq += items[2]
                            al_query_end = int(items[3])
                            
                            p = rawline.find(items[2])
                            al_narrative += lastline[p:p + len(items[2])]
                            
                        elif len(items) == 7 and jmatch == items[3]:
                            al_jmatch_end = int(items[6])
                            if al_jmatch_start == 0:
                                al_jmatch_start = int(items[4])
                            if al_query_j_start == 0:
                                query_frag_gaps = 0
                                for i in range(0, len(items[5])):
                                    if al_query_frag[i] == '-':
                                        query_frag_gaps += 1
                                    if items[5][i] != '-':
                                        al_query_j_start = al_query_frag_start + i - query_frag_gaps
                                        break
                            
                        elif 'Lambda' in line:
                            alignments = False
                            if len(al_query_seq) != len(al_narrative) and len(al_narrative) > 0 and len(al_query_seq) > 0:
                                print 'Error: alignment sequence and narrative are out of synch in id %s' % res['Sequence ID']
                            
                            # field labels may be truncated if the fields are short. In the worst case, this can be right down
                            # to a single letter: <C> or <F>. I have not seen any instances of anything shorter than that, but
                            # we should check.
                            
                            def findextents(s_needle, e_needle, haystack):
                                extents = []
                                i = 0
                                haystack_len = len(haystack)
                                while i < haystack_len:
                                    i = haystack.find(s_needle, i)
                                    if i < 0:
                                        return extents
                                    j = haystack.find(e_needle, i + len(s_needle)) + 1
                                    if j < 0:
                                        if verbose:
                                            print "Broken field delimiter in %s: no closing %s." % (haystack, e_needle)
                                        return extents
                                    extents.append((i, j))
                                    i = j
                            
                            extents = findextents('<', '>', al_narrative)
                            
                            if extents is not None:
                                field_values = []
                                field_indeces = []
                                for (start, end) in extents:
                                    fdesc = al_narrative[start:end]
                                    field_values.append(al_query_seq[start:end])
                                    if '<' in fdesc[1:]:
                                        if verbose:
                                            print "Broken field delimiter in %s: unexpected '<" % al_narrative
                                        break
                                    for i, v in enumerate(fieldorder):
                                        fval = -1
                                        if v in fdesc:
                                            fval = i
                                            break
                                    field_indeces.append(fval)
        
                                    
                                # we require at least one unambiguous field
                                start_ind = -1
                                for i, v in enumerate(field_indeces):
                                    if v >= 0:
                                        start_ind = v - i
                                        break
                                        
                                fields_correct = True
                                if start_ind < 0:
                                    if verbose:
                                        print "Can't infer any fields in %s" % al_narrative
                                else:
                                    for i,v in enumerate(field_indeces):
                                        if v < 0:
                                            field_indeces[i] = start_ind + i
                                        elif v != start_ind + i:
                                            if verbose:
                                                print "Fields out of order, or missing field, in %s" % al_narrative
                                            fields_correct = False
                                            break
                                                                                                                
                                if fields_correct:
                                    i = 0
                                    for ind in field_indeces:
                                        res[fieldorder[ind] + '-IMGT'] = field_values[i].replace('-', '')
                                        if fieldorder[ind] == 'FR3':
                                            raw_fr3 = field_values[i]       # store this so that we can locate it in the query, even if gapped
                                        i += 1
                            
                            if res['Sequence ID'] not in seqs:
                                raise ValueError("Sequence %s not found in FASTA file" % res['Sequence ID']) 
                                
                            seq = seqs[res['Sequence ID']]
                            res['nt Sequence'] = seq
                                
                            # Sometimes IgBlast fails to match the entire FR1 region.
                            # Extend FR1 back to cover as much as possible of the v-gene.
                            # Remember indeces reported by IgBlast are 1-based not 0-based.
                            if 'FR1-IMGT' in res and al_query_start > 1 and al_vmatch_start > 1:
                                extent = min(al_query_start, al_vmatch_start)
                                res['FR1-IMGT'] = seq[al_query_start - extent - 1: al_query_start - 1] + res['FR1-IMGT']
                                while (len(res['FR1-IMGT']) % 3) != 0:
                                    res['FR1-IMGT'] = res['FR1-IMGT'][1:]
                                    
                            # Determine CDR3 and Junction
                            # Need to make sure we work with the ungapped query sequence throughout, so that alignment positions
                            # reported by IgBlast are correctly interpreted.
                            if jmatch != '' and 'FR3-IMGT' in res:
                                j_seq = j_germs[jmatch]
                                al_query_seq_ungapped = al_query_seq.replace('-', '')
                                CDR3_start = al_query_seq_ungapped.find(res['FR3-IMGT']) + len(res['FR3-IMGT']) + al_query_start
                                trailer = al_query_seq_ungapped[CDR3_start - al_query_start:]
            
                                # Extend the sequence forward to cover as much as possible of the j-gene.
                                if al_jmatch_start > 0 and al_jmatch_end < len(j_seq) and al_query_end < len(seq):
                                    extent = min(len(j_seq) - al_jmatch_end, len(seq) - al_query_end)
                                    trailer += seq[al_query_end:al_query_end + extent]
                                    
                                # Make sure the sequence doesn't extend beyond the j-gene.
                                if len(trailer[al_query_j_start - CDR3_start:]) > len(j_seq[al_jmatch_start - 1:]):
                                    trailer = trailer[:len(j_seq[al_jmatch_start - 1:]) - len(trailer[al_query_j_start - CDR3_start:])]
                                    
                                # Trim to a whole number of codons
                                while len(trailer) % 3 != 0:
                                    trailer = trailer[:-1]
                                    
                                # Scan the j-region for the conserved W or F, checking that it is present both in the sequence and
                                # in the germline J-gene.
            
                                def chunks(l, n):
                                    '''Yield successive n-sized chunks from l.'''
                                    for i in xrange(0, len(l), n):
                                        yield l[i:i + n]
            
                                if 'VH' in res['Chain Type']:
                                    sig = 'TGG'
                                else:
                                    sig = 'TT[CT]'
            
                                # Start at the first codon boundary that is into the J-gene
                                start_q_pos = al_query_j_start - CDR3_start
                                start_j_pos = al_jmatch_start - 1
                                while start_q_pos % 3 != 0:
                                    start_q_pos += 1
                                    start_j_pos += 1
                                    
                                j_q_trailer = trailer[start_q_pos:]
                                j_s_trailer = j_seq[start_j_pos:]
    
                                found_CDR3_end = False
                                CDR3_end = start_q_pos
                                for q, s in zip(chunks(j_q_trailer, 3), chunks(j_s_trailer, 3)):
                                    if re.match(sig, q) and re.match(sig, s):
                                        found_CDR3_end = True
                                        break
                                    CDR3_end += 3
            
                                if found_CDR3_end:
                                    res['CDR3-IMGT'] = trailer[:CDR3_end]
                                    res['JUNCTION-IMGT'] = res['FR3-IMGT'][-3:] + trailer[:CDR3_end + 3]
                                    res['FR4-IMGT'] = trailer[CDR3_end:]
                                else:
                                    res['Notes'] += 'Closing CDR3 F/W not found.'
                                    res['Functionality'] = 'unproductive'
                        
                lastline = rawline
            
            except:
                id = res['Sequence ID'] if 'Sequence ID' in res else '<unknown>'
                print "Error parsing sequence %s:\n%s" % (id, traceback.format_exc())

        for k,v in res.iteritems():
            if 'N/A' in v:
                res[k] = ''         # blank out any N/As from IgBLAST for compatibility with IMGT
        
        writer_n.writerow(res)
        res = translate_res(res)        
        writer_aa.writerow(res)        

def translate_res(res):
    fieldorder = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT', 'FR4-IMGT']
    for k, v in res.iteritems():
        if 'FR' in k or 'CDR' in k or 'JUNCTION' in k:                        
            if len(v) % 3 != 0:
                # Probably we have the tail-end of a sequence. See whether we have stuff upstream or downstream.
                if k in fieldorder:
                    ind = fieldorder.index(k)
                    if ind == 0 or fieldorder[ind-1] not in res:
                        while len(v) % 3 != 0:
                            v = v[1:]
                    else:
                        while len(v) % 3 != 0:
                            v = v[:-1]
                else:  # junction
                    if 'FR3-IMGT' not in res:
                        while len(v) % 3 != 0:
                            v = v[1:]
                    else:
                        while len(v) % 3 != 0:
                            v = v[:-1]
                res['Notes'] += '%s incomplete or out-of-frame. ' % k   
            
            try:    
                res[k] = str(Seq(v).translate())
            except:
                print "Error translating field %s: (%s)\n%s." % (k, v, traceback.format_exc())
                
    return res


if __name__=='__main__':
    main()
