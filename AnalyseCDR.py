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

# A class to manage and interpret CDR definitions specified against an Alignment.

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

import os
import sys
from Alignment import Alignment

class AnalyseCDR():
    """
    Create the class against the specified alignment. Cdr positions are specified either as a list of ranges, or 
    in a file. If in a file, the file should contain a single line consisting of six positions, separated by commas.
    >>> msa = Alignment("testfiles/test_cdr.fasta")
    >>> pos = msa.read_position_numbers(file_name="testfiles/test_cdr_seqnum.txt")
    >>> acdr = AnalyseCDR(msa, file_name="testfiles/test_cdr_range.txt")
    >>> print acdr.get_cdrs(msa[0])
    ['QEEEEK', 'QGGGGK', 'QLLLLK']
    """
    def __init__(self, alignment, cdrs=None, file_name=None):
        self.cdrs = cdrs
        if file_name:
            self.__read_cdrs(file_name)    
        self.alignment = alignment
        self.__adjust_cdrs()
        
        
    def __adjust_cdrs(self):
        if self.cdrs and self.alignment:
            # Adjust cdr3 upper bound, to allow for the possibility that in the file it might be set higher than the end of the alignment we have
    
            if self.cdrs[2][1] != "0" and not self.alignment.within_positions(self.cdrs[2][1]):
                self.cdrs[2][1] = self.alignment.max_position()
    
            # Check all lower bounds, and adjust if the lower bound is outside the alignment, but the upper bound is inside
    
            for i in range(3):
                if self.alignment.within_positions(self.cdrs[i][1]) and not self.alignment.within_positions(self.cdrs[i][0]):
                    self.cdrs[i][0] = self.alignment.min_position()
    
            self.range_cdr = []
            for i in range(3):
                if self.alignment.within_positions(self.cdrs[i][0]) and self.alignment.within_positions(self.cdrs[i][1]):
                    self.range_cdr.append(range(self.alignment.index_of(str(self.cdrs[i][0])), self.alignment.index_of(str(self.cdrs[i][1]))+1))
                else:
                    self.range_cdr.append(range(-1, -1))


    def analyse(self):
        """
        Return an html table showing which CDR positions in the alignment set are common to germline (the first record),
        common across all sequences except the germline, and which vary between sequences
        >>> msa = Alignment("testfiles/test_cdr.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/test_cdr_seqnum.txt")
        >>> acdr = AnalyseCDR(msa, file_name="testfiles/test_cdr_range.txt")
        >>> print acdr.analyse()
        <table class='table' border=1>
        <tr><th></th><th>Conserved to Germline</th><th>Common to All</th><th>Variations Across Samples</th></tr>
        <tr><td>CDR1</td><td>20Q,23E,24E,25K</td><td>21H,22H</td><td></td></tr>
        <tr><td>CDR2</td><td>30Q,32G,33G,34G,35K</td><td></td><td>31:AS</td></tr>
        <tr><td>CDR3</td><td>40Q,41L,42L,43L,44L,45K</td><td></td><td></td></tr>
        </table>
        <BLANKLINE>
        """
        common_to_germline = []
        for i in range(0, len(self.alignment[0])):
            value = None
            common = True
            for rec in self.alignment:
                if not "node #" in rec.id:
                    if value is None:
                        value = rec.seq[i]
                    elif rec.seq[i] != value:
                        common = False
                        break
            common_to_germline.append(common)

        common_to_descendents = []
        for i in range(0, len(self.alignment[0])):
            value = None
            common = True
            for (ind, rec) in enumerate(self.alignment):
                if ind != 0 and not "node #" in rec.id:
                    if value is None:
                        value = rec.seq[i]
                    elif rec.seq[i] != value:
                        common = False
                        break
            common_to_descendents.append(common)


        ctg = [[], [], []]
        ctd = [[], [], []]
        leaves = [[], [], []]

        for i in range(3):
            for j in range(0, len(self.alignment[0])):
                if j in self.range_cdr[i]:
                    if common_to_germline[j]:
                        ctg[i].append(self.alignment.position_of(j) + str(self.alignment[0].seq[j]))
                    elif common_to_descendents[j]:
                        ctd[i].append(self.alignment.position_of(j) + str(self.alignment[1].seq[j]))
                    else:
                        variants = []
                        for (ind, rec) in enumerate(self.alignment):
                            if ind != 0 and not "node #" in rec.id:
                                if rec.seq[j] not in variants:
                                    variants.append(rec.seq[j])
                        variants.sort()
                        leaves[i].append(self.alignment.position_of(j) + ":" + "".join(variants))

        output = "<table class='table' border=1>\n<tr><th></th><th>Conserved to Germline</th><th>Common to All</th><th>Variations Across Samples</th></tr>\n"

        for i in range(3):
            output += "<tr><td>CDR" + str(i+1) + "</td><td>" + ",".join(ctg[i]) + "</td><td>" + ",".join(ctd[i]) + "</td><td>" + ",".join(leaves[i]) + "</td></tr>\n"

        output += "</table>\n"
        return output

    def __read_cdrs(self, cdrfile):
        fi = open(cdrfile, "r")
        line = fi.readline().strip().replace(" ", "")
        
        if len(line) == 0:
            return            
        
        self.cdrs = [[], [], []]
        incdrs = line.split(",")
        if len(incdrs) != 6:
            raise AttributeError("CDR file must contain exactly 6 positions on one line, separated by commas")

        for i in range(3):
            a = str(incdrs[i*2])
            b = str(incdrs[i*2 + 1])
            self.cdrs[i].append(a)
            self.cdrs[i].append(b)

        return

    # Return a vector containing the CDRs in the specified record of the alignment, given the cdr position numbers
    def get_cdrs(self, rec):
        # Be flexible about the end point of CDR3
        if self.cdrs[2][1] != "0" and not self.alignment.within_positions(self.cdrs[2][1]):
            self.cdrs[2][1] = self.alignment.max_position()

        ret = []
        for i in range(3):
            if self.alignment.within_positions(self.cdrs[i][0]) and self.alignment.within_positions(self.cdrs[i][1]):
                a = self.alignment.index_of(self.cdrs[i][0])
                b = self.alignment.index_of(self.cdrs[i][1])
                ret.append(str(rec.seq)[a:b+1])
            else:
                ret.append("")

        return ret
    
    
    def categorize_index(self, index):
        """
        Categorise an index into the alignment 
        returns:
        1..3 if the position is in CDR1..3
        11..14 if it is in FR1..4
        >>> msa = Alignment("testfiles/test_cdr.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/test_cdr_seqnum.txt")
        >>> acdr = AnalyseCDR(msa, file_name="testfiles/test_cdr_range.txt")
        >>> print acdr.categorize_index(0)
        11
        >>> print acdr.categorize_index(10)
        1
        >>> print acdr.categorize_index(35)
        3
        >>> print acdr.categorize_index(36)
        14
        """
        for i in range(3):
            if index in self.range_cdr[i]:
                return i+1
            
        for i in range(3):
            if len(self.range_cdr[i]) > 0 and index < self.range_cdr[i][0]:
                return i+11
            
        return 14
    
    def category_diff(self, id1, id0=0):
        """
        Express the difference between sequence id1 and the baseline sequence id0 as a text string showing the
        count of differences in each CDR and FR. If id0 is not specified, the first sequence in the set is used
        >>> msa = Alignment("testfiles/test_cdr.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/test_cdr_seqnum.txt")
        >>> acdr = AnalyseCDR(msa, file_name="testfiles/test_cdr_range.txt")
        >>> print acdr.category_diff("1", "2")
        1(2)0(1)0(0)0
        >>> print acdr.category_diff("1", "1")
        <BLANKLINE>
        """
        totals = {}
        seq1 = self.alignment[0].seq if not id0 else self.alignment.get_record(id0).seq
        seq2 = self.alignment.get_record(id1).seq
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                cat = self.categorize_index(i)
                totals[cat] = totals.get(cat, 0) + 1
                
        if len(totals) > 0:
            return "%d<%d>%d<%d>%d<%d>%d" % (totals.get(11,0), totals.get(1,0), totals.get(12,0),
                                             totals.get(2,0), totals.get(13,0),
                                             totals.get(3,0), totals.get(14,0))
        else:
            return ""


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    


