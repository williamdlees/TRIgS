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
#
# Tests for ExtractFromIMGT.py

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import subprocess
import os
import csv

MODEL_OUTPUT_FILE = 'summarise_model_n.txt'
TEST_IGBLAST_OUTPUT = 'summarise_test.out'
TEST_SEQUENCES = 'summarise_test.fasta'
JGERM_FILE = 'rabbit_IGj.fa'
TEST_OUTPUT_FILE = 'summarise_test_n.txt'

try:
    os.remove(TEST_OUTPUT_FILE)
except:
    pass

subprocess.call("python ../IgBLASTPlus.py %s %s %s summarise_test" % (TEST_SEQUENCES, JGERM_FILE, TEST_IGBLAST_OUTPUT), shell=True, cwd=os.getcwd())

with open(MODEL_OUTPUT_FILE, 'r') as model, open(TEST_OUTPUT_FILE, 'r') as test:
    mreader = csv.DictReader(model, delimiter='\t')
    treader = csv.DictReader(test, delimiter='\t')
    for mrow in mreader:
        trow = treader.next()
        
        print "testing record %s --> %s" % (mrow['Sequence ID'], trow['Sequence ID'])
        for mk, mv in mrow.iteritems():
            if trow[mk] != mv:
                print "   Difference in field %s: model %s  test %s" % (mk, mv, trow[mk])
        


