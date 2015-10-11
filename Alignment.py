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

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO
from Bio.Alphabet import IUPAC

class Alignment(MultipleSeqAlignment):
    """A class to manage a set of aligned sequences. This wraps BioPython's MultipleSeqAlignment class."""
    def __init__(self, file_name=None, data = None, format='fasta'):
        if file_name:
            super(Alignment, self).__init__(AlignIO.read(file_name, format))
        elif data:
            super(Alignment, self).__init__(AlignIO.parse(StringIO(data), format))
        else:
            super(Alignment, self).__init__([])

    def read_nt(self, file_name, format='fasta'):
        """Read a set of nucleotide sequences from a file, translate to amino acid sequences and store in the
        Alignment.

        Existing records in the Alignment, if any, are preserved.

        The nt sequence must represent a whole number of codons and must start on a codon boundary. Whole codon gaps
        (represented by ---) are permitted. Stop codons will be represented by * in the translation.

        >>> msa = Alignment()
        >>> msa.read_nt("testfiles/dna_with_gaps_and_noncoding.fasta")
        >>> print len(msa)
        8
        >>> msa = Alignment()
        >>> msa.read_nt("testfiles/rna_with_gaps_and_noncoding.fasta")
        >>> print len(msa)
        1
        >>> msa = Alignment()
        >>> msa.read_nt("testfiles/dna_with_ambiguous_codons.fasta")
        >>> print len(msa)
        1
        >>> msa = Alignment()
        >>> msa.read_nt("testfiles/rna_with_unaligned_insertion.fasta")  # doctest: +IGNORE_EXCEPTION_DETAIL
        >>> msa = Alignment()
        >>> msa.read_nt("testfiles/rna_with_incomplete_codon.fasta")
        Traceback (most recent call last):
        ...
        AttributeError: Sequence must represent a whole number of codons.
        """
        def chunks(l, n):
            """ Yield successive n-sized chunks from l."""
            for i in xrange(0, len(l), n):
                yield l[i:i+n]

        nsa = AlignIO.read(file_name, format)
        for rec in nsa:
            if len(rec.seq) % 3 != 0:
                raise AttributeError("Sequence must represent a whole number of codons.")

            t = ""
            for codon in chunks(rec.seq, 3):
                if '-' in codon:
                    t += '-'
                else:
                    t += codon.translate()
            self.append(SeqRecord(t, rec.id))

    def read_position_numbers(self, file_name):
        """Read position number definitions from a file.

        Position file format:

        The first line of the file should contain a single number, which is the number of the first residue.
        or <number>,<letter>, signifying that the first residue should be numbered as an insertion, eg 85,A.

        Subsequent lines take one of the following formats:
        <number>, - a number followed by a minus sign indicates that the number should be skipped.
        <number>, <letters or numbers> - the numbered residue is followed by an insertion
          - hence 99,A would indicate that residue 99 is followed by residue 99A
          - 100,11 would indicate that residue 100 is followed by residue 100.11
          - a subsequent 99, B would indicate that 99A is followed by 99B
        <number>; <letters or numbers> - the numbered residue is preceded by an insertion (as in the IMGT junction)

        Numbers in the file must be in ascending order.

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> msa.read_position_numbers(file_name="testfiles/seqnum_with_insertions.txt")
        ['46', '48,A', '49,-']
        """
        first = True
        self.position_numbers = []

        f = open(file_name, "r")
        for line in f:
            line = line.strip()
            if line != '':
                rec = line.split(";") if ";" in line else line.split(",")
                if(first and len(rec) == 1):
                    self.position_numbers.append(str(rec[0]))
                elif len(rec) == 2:
                    self.position_numbers.append(line)
                else:
                    raise AttributeError('Bad file format')
                first = False
        f.close()
        self.__apply_position_numbers()
        return self.position_numbers

    def write_position_numbers(self, file_name):
        """Save a set of position numbers to a file. The format is that used by read_position_numbers.

        >>> msa = Alignment()
        >>> msa.set_position_numbers(position_numbers = ['46', '48,-', '48,A', '49,-'])
        >>> msa.write_position_numbers("testfiles/test_positions.txt")
        >>> msa.read_position_numbers("testfiles/test_positions.txt")
        ['46', '48,-', '48,A', '49,-']
        >>> msa = Alignment()
        >>> msa.set_position_numbers(position_numbers = ['100', '100,1', '100,2', '101;2', '101;1'])
        >>> msa.write_position_numbers("testfiles/test_positions.txt")
        >>> msa.read_position_numbers("testfiles/test_positions.txt")
        ['100', '100,1', '100,2', '101;2', '101;1']
        """
        f = open(file_name, "w")

        for pos in self.position_numbers:
            f.write(pos + "\n")

        f.close()

    def set_position_numbers(self, first=1, position_numbers=None):
        """Define position numbers for the alignment. If a list is specified, position specifications from the list are
        applied (the format is that created by read_position_numbers). Otherwise positions are numbered sequentially
        starting from 'first'. 

        >>> msa = Alignment()
        >>> msa.set_position_numbers(10)
        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> msa.set_position_numbers(10)
        >>> print(msa.seqnum)
        ['10', '11', '12', '13', '14']

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> msa.set_position_numbers(position_numbers = ['46', '48,-', '48,A', '49,-'])
        >>> print(msa.seqnum)
        ['46', '47', '48A', '50', '51']
        """
        if position_numbers:
            self.position_numbers = position_numbers
        else:
            self.position_numbers = [str(first)]
        self.__apply_position_numbers()

    def __apply_position_numbers(self, first=1):
        """Apply the position numbers definition to the current alignment.
       """
        if len(self) < 1:
            return  # no alignment

        length = len(self[0])

        if not hasattr(self, 'position_numbers') or not self.position_numbers:
            self.seqnum = [str(i) for i in range(first, length+first)]
            return

        first = True
        self.seqnum = []
        for line in self.position_numbers:
            rec = line.split(";") if ";" in line else line.split(",")
            if first:
                if len(rec) == 1:
                    self.seqnum.append(str(rec[0]))
                elif len(rec) == 2:
                    self.seqnum.append(str(rec[0]) + str(rec[1]))
                else:
                    raise AttributeError("Bad syntax in sequence number file: " + line)
                lastnum = int(rec[0])
                first = False

            elif len(rec) == 2 and rec[0].isdigit():
                if int(rec[0]) < lastnum:
                    raise AttributeError("Numbers in position file must be in ascending order")
                if rec[1] == "-":
                    while lastnum < int(rec[0]) - 1:
                        lastnum += 1
                        self.seqnum.append(str(lastnum))
                    lastnum = int(rec[0])
                else:
                    if len(rec[1]) != 1 and not ('.' in rec[1]):
                        raise AttributeError("Bad syntax in sequence number file: " + line)                        
                    if "," in line:                    
                        while lastnum < int(rec[0]):
                            lastnum += 1
                            self.seqnum.append(str(lastnum))
                        self.seqnum.append(str(lastnum) + str(rec[1]))
                    elif ";" in line:
                        while lastnum < int(rec[0])-1:
                            lastnum += 1
                            self.seqnum.append(str(lastnum))
                        self.seqnum.append(str(lastnum+1) + str(rec[1]))
                    else:
                        raise AttributeError("Bad syntax in sequence number file: " + line)
            else:
                raise AttributeError("Bad syntax in sequence number file: " + line)

        while len(self.seqnum) < length:
            lastnum += 1
            self.seqnum.append(str(lastnum))

    def position_of(self, ind):
        """Return the position of a residue in the string (as indexed in the alignment), as described by the
        position numbers.

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers("testfiles/seqnum_with_insertions.txt")
        >>> print(msa.position_of(3))
        48A
        """
        return self.seqnum[ind]

    def index_of(self, pos):
        """Return the index of a residue in the string, given its position (as defined in the position file)

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers("testfiles/seqnum_with_insertions.txt")
        >>> print(msa.index_of('48A'))
        3
        """
        return self.seqnum.index(pos)

    def min_position(self):
        """Return the lowest position used in the range

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/seqnum_with_insertions.txt")  # doctest: +NORMALIZE_WHITESPACE
        >>> print msa.min_position()
        46
        """
        return self.position_of(0)

    def max_position(self):
        """Return the highest position used in the range

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/seqnum_with_insertions.txt")  # doctest: +NORMALIZE_WHITESPACE
        >>> print msa.max_position()
        50
        """
        return self.position_of(self.get_alignment_length() -1)

    def within_positions(self, pos):
        """Return True if pos is a valid position, False otherwise

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/seqnum_with_insertions.txt")  # doctest: +NORMALIZE_WHITESPACE
        >>> print msa.within_positions("29A")
        False
        >>> print msa.within_positions("50")
        True
        """
        return pos in self.seqnum

    def _format_rep_header(self, offset=0):
        """Internal function to create numbered header for pretty print output

        Note: this test depends on there being whitespace at the end of lines. Make sure the editor does not remove it.

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/seqnum_with_insertions.txt")  # doctest: +NORMALIZE_WHITESPACE
        >>> print msa._format_rep_header(3)
               5
           67880
              A 
        <BLANKLINE>
        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/seqnum_with_imgt numb.txt")  # doctest: +NORMALIZE_WHITESPACE
        >>> print msa._format_rep_header(3)
           66777
            111 
             2  
        <BLANKLINE>
        """
        # Find the largest number ending in zero, ie the largest that needs a full label
        maxlab = 0    # the highest sequence number that gets a full label
        sublab = 0    # the number of rows to allow for suffixes
        for num in self.seqnum:
            if "." in num:
                subnum = str(num).split(".")
                if subnum[0][-1:] == "0":
                    maxlab = num
                sublab = max(sublab, len(subnum[1]))
            elif not str(num)[-1:].isdigit():
                sublab = max(sublab, 1)
            elif str(num)[-1:] == "0":
                maxlab = num

        labrows = len(str(maxlab))
        ret = ""

        for row in range(labrows-1, 0-sublab-1, -1):
            if offset > 0:
                ret += " ".ljust(offset)
            if row > 0:
                for num in self.seqnum:
                    if str(num)[-1:] == "0" and len(str(num)) > row and not "." in str(num):
                        ret += str(num)[len(str(num))-row-1]
                    else:
                        ret += " "
            elif row == 0:
                for num in self.seqnum:
                    if "." in num:                  #100.1
                        subnum = str(num).split(".")
                        ret += subnum[0][-1:]
                    elif num[-1:].isdigit():        #100
                        ret += str(num)[-1:]
                    else:                            #100A
                        ret += str(num)[-2:-1]
            else:
                for num in self.seqnum:
                    if "." in num:
                        subnum = str(num).split(".")
                        if len(subnum[1]) >= 0-row:
                            ret += subnum[1][len(subnum[1]) + row -1]
                        else:
                            ret += " "
                    elif row == -1 and not (num[-1:]).isdigit():
                        ret += str(num)[-1:]
                    else:
                        ret += " "
                    
            ret += "\n"
        return ret

    def report(self, cols=0):
        """Generate a pretty print report, showing the numbering of residues and identifying substitutions
            cols if specified, denotes the maximum width of the report, which will be wrapped if necessary

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> print(msa.report())
          12345
        1 IATCT
        2 ..AD.
        3 ...D.
        4 ....S
        5 ...D.
        <BLANKLINE>

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> print(msa.report(4))
          12
        1 IA
        2 ..
        3 ..
        4 ..
        5 ..
        <BLANKLINE>
          34
        1 TC
        2 AD
        3 .D
        4 ..
        5 .D
        <BLANKLINE>
          5
        1 T
        2 .
        3 .
        4 S
        5 .
        <BLANKLINE>
        """
        if not hasattr(self, 'seqnum'):
            self.__apply_position_numbers()
        maxidlen = 0
        for rec in self:
            maxidlen = max(len(rec.id), maxidlen)

        ret = self._format_rep_header(maxidlen+1)
        fmt = "%(id)-" + str(maxidlen+1) + "s"

        for recno, rec in enumerate(self):
            ret += fmt % {"id": rec.id}
            for i, c in enumerate(rec.seq):
                if recno != 0 and rec.seq[i] == self[0].seq[i]:
                    ret += "."
                else:
                    ret += str(rec.seq[i]).upper()
            ret += "\n"

        if cols > 0:
            ret = self.__splitlines(ret, cols, maxidlen+1)

        return ret

    @staticmethod
    def __splitlines(report, maxlength, label_cols):
        """
        Split the report (which is assumed to consist of lines of equal length) into a longer report in which each
        line is maxlength or less. name_cols specifies the width of the label field, which is repeated at the start
        of each line.
        """

        def chunks(l, n):
            """ Yield successive n-sized chunks from l."""
            for i in xrange(0, len(l), n):
                yield l[i:i+n]

        inlines = report.split("\n")[:-1]
        labels = [line[:label_cols] for line in inlines]
        data = [line[label_cols:] for line in inlines]
        outlines = []

        for chunk in chunks(zip(*data), maxlength-label_cols):
            a = ["".join(line) for line in zip(*chunk)]
            outlines.extend(["".join(line) for line in zip(labels, a)])
            outlines.extend(" ")

        return "\n".join(["".join(line) for line in outlines])

    def seqdiff(self, id1, id0=None):
        """
        Express the difference between sequence id1 and the baseline sequence id0 as a set of mutations
        If id0 is not specified, the first sequence in the set is used

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/seqnum_with_insertions.txt")  # doctest: +NORMALIZE_WHITESPACE
        >>> msa.seqdiff('2')
        set(['T48A', 'C48AD'])

        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> pos = msa.read_position_numbers(file_name="testfiles/seqnum_with_insertions.txt")  # doctest: +NORMALIZE_WHITESPACE
        >>> msa.seqdiff('2','3')
        set(['T48A'])
        """
        diffs = []

        seq1 = self[0].seq if not id0 else self.get_record(id0).seq
        seq2 = self.get_record(id1).seq
        for i, c in enumerate(seq1):
            if seq1[i] != seq2[i] and seq1[i] != '-' and seq2[i] != '-':
                diffs.append("%s%s%s" % (seq1[i], self.position_of(i), seq2[i]))

        return set(diffs)

    def get_record(self, id):
        """
        Get a sequence record given its id.
        >>> msa = Alignment("testfiles/test_short_aa.fasta")
        >>> print(msa.get_record('2').seq)
        IAADT
        """
        for rec in self:
            if rec.id == id:
                return rec

        return None

if __name__ == "__main__":
    import doctest
    doctest.testmod()

