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

# A script to extract fields or series of fields from IMGT-style tab-separated files, storing the results as FASTA sequences.

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import argparse


def main():
    parser = argparse.ArgumentParser(description='Extract specified sequences from an IMGT-style CSV file.')
    parser.add_argument("analysis_file", help='analysis file from (IMGT, IgBLASTPlus format)')
    parser.add_argument("required_fields", help='the fields to extract, separated by + (e.g. CDR1-IMGT+FR2-IMGT+CDR2-IMGT). Specified fields are concatenated.')
    parser.add_argument("outfile", help='output file (FASTA)')
    parser.add_argument("-g", "--germline", help='Restrict to v-germlines containing the specified string')
    parser.add_argument("-d", "--dgermline", help='Restrict to d-germlines containing the specified string')
    parser.add_argument("-j", "--jgermline", help='Restrict to j-germlines containing the specified string')
    parser.add_argument("-p", "--productive", help="Restrict to productive sequences", action="store_true")
    parser.add_argument("-u", "--unknown_nuc", help="Do not include sequences including N (nt file only)", action="store_true")
    parser.add_argument("-s", "--stopcodon", help="Do not include sequences including stop codons (AA file only)", action="store_true")
    parser.add_argument("-t", "--chaintype", help="Only include sequences from the specified chain type")
    parser.add_argument("-c", "--complete", help="Ensure that only complete fields are included, by checking that neighbouring fields are populated", action="store_true")
    parser.add_argument("--fr1_lim", help="truncate FR1 to fr1_lim nucleotides or amino acids, and do not include sequences with a shorter FR1 length", type=int)
    parser.add_argument("--fr4_lim", help="truncate FR4 to fr4_lim nucleotides or amino acids, and do not include sequences with a shorter FR4 length", type=int)
    parser.add_argument("-v", "--verbose", help="verbose output", action="store_true")
    args = parser.parse_args()
    
    fieldorder = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT', 'FR4-IMGT']
    verbose = args.verbose
    germline = args.germline
    jgermline = args.jgermline
    dgermline = args.dgermline
    chaintype = args.chaintype
    wanted_fields = args.required_fields.split('+')
    
    firstline = True
    with open(args.analysis_file, 'r') as fi, open(args.outfile, 'wb') as fo:
        reader = csv.DictReader(fi, delimiter='\t')
        ln = fi.readline()
        sep = ("\t" if "\t" in ln else ",")
        fi.seek(0)
        reader = csv.DictReader(fi, delimiter=sep)
        for row in reader:
            errors = False
            if firstline:
                # check the specified fields are valid
                for wanted_field in wanted_fields:
                    if wanted_field not in reader.fieldnames:
                        print 'Error: field %s is not included in the analysis file.' % wanted_field
                        errors = True
                if germline and 'V-GENE and allele' not in row:
                    print "Error: specific v-gene requested, but 'V-GENE and allele' field is not included in the analysis file."                    
                if dgermline and 'D-GENE and allele' not in row:
                    print "Error: specific d-gene requested, but 'D-GENE and allele' field is not included in the analysis file."                    
                if jgermline and 'J-GENE and allele' not in row:
                    print "Error: specific j-gene requested, but 'J-GENE and allele' field is not included in the analysis file."                    
                if errors:
                    exit(1)
                    
                firstline = False

            if germline is not None:
                g = row['V-GENE and allele'].split(",")
                match = False
                for gg in g:
                    if germline in gg:
                        match = True
                        break
                if not match:
                    continue

            if dgermline is not None:
                g = row['D-GENE and allele'].split(",")
                match = False
                for gg in g:
                    if dgermline in gg:
                        match = True
                        break
                if not match:
                    continue

            if jgermline is not None:
                g = row['J-GENE and allele'].split(",")
                match = False
                for gg in g:
                    if jgermline in gg:
                        match = True
                        break
                if not match:
                    continue

            if chaintype is not None:
                g = row['Chain Type']
                if chaintype != g:
                    continue
                    
            if args.productive and 'unproductive' in row['Functionality'] or 'productive' not in row['Functionality']:
                continue
                    
            include_row = True
            for wanted_field in wanted_fields:
                if row[wanted_field] is None or row[wanted_field] == '' or row[wanted_field] == 'null':
                    include_row = False

            if args.complete:
                for k in wanted_fields:
                    if k in fieldorder:
                        ind = fieldorder.index(k)
                        if ind > 0 and (row[fieldorder[ind-1]] == '' or row[fieldorder[ind-1]] == 'null'):
                            include_row = False
                        if ind < len(wanted_fields) - 1 and (row[fieldorder[ind+1]] == '' or row[fieldorder[ind+1]] == 'null'):
                            include_row = False

            if args.fr1_lim != None and len(row['FR1-IMGT']) < args.fr1_lim:
                include_row = False

            if args.fr4_lim != None and len(row['FR4-IMGT']) < args.fr4_lim:
                include_row = False

            if include_row:
                seq = ''
                for wanted_field in wanted_fields:
                    if wanted_field == 'FR1-IMGT':
                        if args.fr1_lim != None:
                            seq += row['FR1-IMGT'][0-args.fr1_lim:]
                        else:
                            seq += row['FR1-IMGT']
                    elif wanted_field == 'FR4-IMGT':
                        if args.fr4_lim != None:
                            seq += row['FR4-IMGT'][0:args.fr4_lim]
                        else:
                            seq += row['FR4-IMGT']
                    else:
                        seq += row[wanted_field]
                
                if args.unknown_nuc and ('N' in seq or 'n' in seq):
                    include_row = False
                if args.stopcodon and ('*' in seq or 'X' in seq or 'x' in seq):
                    include_row = False
                
            if include_row:
                SeqIO.write(SeqRecord(Seq(seq), id=row['Sequence ID'], description = ''), fo, 'fasta')

if __name__ == '__main__':
    main()
 