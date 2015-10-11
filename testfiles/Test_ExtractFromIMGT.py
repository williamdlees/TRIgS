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

NN_FILE = 'extract_test_n.txt'
AA_FILE = 'extract_test_aa.txt'
OUTFILE = 'extract_test_out.fasta'
IMGT_AA_FILE = 'extract_imgt_test_aa.txt'

    
    
    
    
def call_extract(infile, options):
    """
    No options
    >>> call_extract(AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT"')
    1: GFTFSSYWMSWVRQAPGKGLQWVANIKQDGRES
    3: GYTFTSYYMHWVRQAPGQGLEWMGIINPSGGST
    4: GFTVRDNHLSWVRQAPGTGPECVSVIYSXGST
    5: SITSGGYSWTWIRQPPGKGLEWIGYIYHSGTI
    6: GYTFSSYGISWVRQAPGQGLEWMGWISAYNGNT
    7: GFSISSHWMSWVRQAPGKGPQWVANIKQDGSEK
    8: GFTVSSNYMTWVRQAPGKGLELVSVIYSGGST
    9: GYTFTSYNMHWVRQAPGQGLEWMGMIYPAGGGI
    10: GGTFDSFAISWLRQAPGQGPEWMGGIAPALDST
    11: GYTFTSFGISWVRQAPGQGLEWMGWSNGYNTDT
    12: GGNFSSYAIGWVRQAPGQGPEWMGGIIPILGIP
    13: GFTVSSNYMSWVRQAPGKGLEWVSVIYSGGST
    
    IGHV1-69*11 only
    >>> call_extract(AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT" -g "IGHV1-69*11"')
    10: GGTFDSFAISWLRQAPGQGPEWMGGIAPALDST
    
    IGHV1-69 only (not present)
    >>> call_extract(AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT" -g "IGHV1-69"')
    10: GGTFDSFAISWLRQAPGQGPEWMGGIAPALDST
    12: GGNFSSYAIGWVRQAPGQGPEWMGGIIPILGIP
    
    Productive only
    >>> call_extract(AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT" -p')
    1: GFTFSSYWMSWVRQAPGKGLQWVANIKQDGRES
    3: GYTFTSYYMHWVRQAPGQGLEWMGIINPSGGST
    5: SITSGGYSWTWIRQPPGKGLEWIGYIYHSGTI
    6: GYTFSSYGISWVRQAPGQGLEWMGWISAYNGNT
    7: GFSISSHWMSWVRQAPGKGPQWVANIKQDGSEK
    8: GFTVSSNYMTWVRQAPGKGLELVSVIYSGGST
    9: GYTFTSYNMHWVRQAPGQGLEWMGMIYPAGGGI
    10: GGTFDSFAISWLRQAPGQGPEWMGGIAPALDST
    11: GYTFTSFGISWVRQAPGQGLEWMGWSNGYNTDT
    12: GGNFSSYAIGWVRQAPGQGPEWMGGIIPILGIP
    13: GFTVSSNYMSWVRQAPGKGLEWVSVIYSGGST

    No stop codons
    Note: this looks for an X or * in the sequences, not the Yes or No in the Stop Codons column
    >>> call_extract(AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT" -s')
    1: GFTFSSYWMSWVRQAPGKGLQWVANIKQDGRES
    3: GYTFTSYYMHWVRQAPGQGLEWMGIINPSGGST
    5: SITSGGYSWTWIRQPPGKGLEWIGYIYHSGTI
    6: GYTFSSYGISWVRQAPGQGLEWMGWISAYNGNT
    7: GFSISSHWMSWVRQAPGKGPQWVANIKQDGSEK
    8: GFTVSSNYMTWVRQAPGKGLELVSVIYSGGST
    9: GYTFTSYNMHWVRQAPGQGLEWMGMIYPAGGGI
    10: GGTFDSFAISWLRQAPGQGPEWMGGIAPALDST
    11: GYTFTSFGISWVRQAPGQGLEWMGWSNGYNTDT
    12: GGNFSSYAIGWVRQAPGQGPEWMGGIIPILGIP
    13: GFTVSSNYMSWVRQAPGKGLEWVSVIYSGGST


    No unknown nucleotides (note: running on AA file will reject any containing Asn)
    >>> call_extract(AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT" -g "IGHV1-69" -u')
    10: GGTFDSFAISWLRQAPGQGPEWMGGIAPALDST

    Only complete sequences
    >>> call_extract(AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT" -c')
    1: GFTFSSYWMSWVRQAPGKGLQWVANIKQDGRES
    3: GYTFTSYYMHWVRQAPGQGLEWMGIINPSGGST
    4: GFTVRDNHLSWVRQAPGTGPECVSVIYSXGST
    6: GYTFSSYGISWVRQAPGQGLEWMGWISAYNGNT
    7: GFSISSHWMSWVRQAPGKGPQWVANIKQDGSEK
    8: GFTVSSNYMTWVRQAPGKGLELVSVIYSGGST
    9: GYTFTSYNMHWVRQAPGQGLEWMGMIYPAGGGI
    10: GGTFDSFAISWLRQAPGQGPEWMGGIAPALDST
    11: GYTFTSFGISWVRQAPGQGLEWMGWSNGYNTDT
    12: GGNFSSYAIGWVRQAPGQGPEWMGGIIPILGIP
    13: GFTVSSNYMSWVRQAPGKGLEWVSVIYSGGST
            
    Truncate FR1
    >>> call_extract(AA_FILE, '"FR1-IMGT" --fr1_lim 5')
    3: SCKAS
    4: SCGAS
    6: SCKAS
    7: SCVAS
    8: SCAVS
    9: SCKAS
    10: SCKTS
    11: SCKAS
    12: SCKAS
    13: SCVAS
                
    Truncate FR4
    >>> call_extract(AA_FILE, '"FR4-IMGT" --fr4_lim 5')
    1: GQGTM
    2: GQGTL
    3: GQGTL
    5: GQGTL
    6: GQGTL
    7: SQGTL
    8: GQGTT
    9: GQGTK
    10: GQGTT
    11: YFDLW
    12: GQGTL
    13: GQGTP
    
    No unknown nucleotides (note: running on AA file will reject any containing Asn)
    >>> call_extract(IMGT_AA_FILE, '"CDR1-IMGT+FR2-IMGT+CDR2-IMGT" -g "IGHV1-69*06"')
    1: GDTFSNYAFSWVRQAPGQGLEWMGIIIPIFGPA
    2: GDTFSNYAFSWVRQVPGQGLEWMGIIIPIFGPA
    4: GGTFSGYGLSWLRQAPGQGPEWMGIIIPIFGPP
    5: GDTFSNYAFSWVRQAPGQGLEWMGIIIPIFGPA
    7: GDTFSNYAFSWVRQAPGQGLEWMGIIIPIFGPA
    """
    args = infile + " " + options + " " + OUTFILE
    
    
    try:
        os.remove(OUTFILE)
    except:
        pass
    
    subprocess.call("python ../ExtractFromIMGT.py " + args, shell=True, cwd=os.getcwd())

    for seq_record in SeqIO.parse(OUTFILE, 'fasta'):
        print "%s: %s" % (seq_record.id, str(seq_record.seq))
        
    return
    


if __name__ == '__main__':
    import doctest
    doctest.testmod()
