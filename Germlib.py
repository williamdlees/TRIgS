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

# A class to perform operations on a set of IMGT germline library files.


__author__ = 'William Lees'
__docformat__ = "restructuredtext en"


from StringIO import StringIO
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class Germlib:
    """A class to manage germline libraries downloaded from IMGT.

        The constructor takes as argumenta:
        - the pathname of a file containing the germline definitions, in FASTA
          format, with IMGT format headers. The IMGT germline library is downloadable in this form from
          http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP
        - the species name required (as used in the third field of the FASTA header)
        .
        >>> gl = Germlib("Oryctolagus cuniculus", germline_file="testfiles/imgt_germlines.fasta")
        >>> print len(gl.records)
        6581
    """
    def __init__(self, species_name, germline_file=None, germline_data=None):
        if species_name is None:
            raise AttributeError('No species name specified')

        self.species_name = species_name
        
        if germline_file:
            self.records = list(SeqIO.parse(germline_file, "fasta"))
        elif germline_data:
            self.records = list(SeqIO.parse(StringIO(germline_data), "fasta"))
        else:
            raise AttributeError('No library file or libray data provided')

    def seq(self, name):
        """Find the sequence of an allele given its name

        >>> gl = Germlib("Oryctolagus cuniculus", germline_file="testfiles/imgt_germlines.fasta")
        >>> print gl.seq("IGHD3-2*01").seq
        cctatggggtcctggttcccatggctatggggtt
        """
        for rec in self.records:
            if rec.description.split("|")[1] == name and rec.description.split("|")[2].capitalize() == self.species_name.capitalize():
                return rec

        return None
    
    def enumerate_species(self):
        """Return a list of all species referenced in the germline file.
        
        >>> gl = Germlib("", germline_file="testfiles/imgt_germlines.fasta")
        >>> print gl.enumerate_species()[0]
        Bos taurus
        """
        species = []
        for rec in self.records:
            sn = rec.description.split("|")[2]
            if not sn in species:
                species.append(sn)
                
        return species

    def match(self, name, matchseq):
        """Find the best match of matchseq in the specified allele. Return the corresponding germline fragment.
        The alignment score is calculated by counting +5 for each nucleotide match and -4 for each mismatch, which are
        the parameters used by IMGT. Gaps, extensions score -10 (this is my score since I don't know what IMGT uses).
        
        match returns two strings. The first is the best match to matchseq, as above. The second is a 'match string',
        aligned to the target, containing the characters m (match) n (nonmatch) d (deletion) i (insertion), illustrating
        how the target sequence differs to matchseq at each position.

        >>> gl = Germlib("Oryctolagus cuniculus", germline_file="testfiles/imgt_germlines.fasta")
        >>> print gl.match("IGHV1S45*01", "tgtgcgagag")
        ('tgtgcgagag', 'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddmmmmmmmmmmd')
        >>> print gl.match("IGHJ4*01", "ttgtggggcccgggcaccctggtcaccgtctcctcag")
        ('ttgtggggcccaggcaccctggtcaccgtctcctcag', 'dddddddddddmmmmmmmmmmmnmmmmmmmmmmmmmmmmmmmmmmmmm')
        >>> print gl.match("IGHV1S40*01", "tgtgcga")
        ('tgtgcga', 'dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddmmmmmmmdddd')
        >>> print gl.match("IGHV1S40*01", "tgtgcgag")
        ('tgtgcgag', 'dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddmmmmmmmmddd')
        """
        try:
            targseq = str(self.seq(name).seq)
        except:
            raise ValueError('Germline not found')

        matchseq = str(matchseq)

        # Return an alignment that matches the length of the match sequence, and is best matched to it
        # We work back from the 3' end so that we know we have a good alignment at that end, which is important
        # for the V region.
        alignments = pairwise2.align.localms(matchseq, targseq, 5, -4, -10, -10, one_alignment_only=False)

        best = None
        best_alignment = None
        best_matchstr = ""
        max_match = 0
        for alignment in alignments:
            (align1, align2, score, begin, end) = alignment

            seq = ""
            # Skip over leading gaps in the target - but keep them in case we need some to achieve
            # the correct length match after taking account of deletions.
            skip = True
            skipped = []
            for a,b in reversed(zip(align1, align2)):   # a is in matchseq, b is in targ
                if a != '-':
                    skip = False
                if skip:
                    skipped.append(b)
                if not skip and len(seq) < len(matchseq):
                    if a != '-': 
                        if b != '-':
                            seq = b + seq
                        else:
                            seq = a + seq

            while len(seq) < len(matchseq) and len(skipped) > 0:
                seq += skipped.pop()

            match = self.matches(seq, matchseq)
            if match > max_match and len(seq) == len(matchseq):
                best = seq
                best_alignment = alignment
                max_match = match
                for a,b in zip(align1, align2):
                    if b == '-':
                        best_matchstr += 'i'
                    elif a == '-':
                        best_matchstr += 'd'
                    elif a == b:
                        best_matchstr += 'm'
                    else:
                        best_matchstr += 'n'

        if best is None:
            raise ValueError('Unable to match sequence with germline %s' % name)

        #(align1, align2, score, begin, end) = best_alignment
        #print(format_alignment(align1, align2, score, begin, end))
        #print matchseq
        #print name
        #print matchseq
        #print best
        return best, best_matchstr

    def match_from_aa(self, name, matchseq):
        """
        This function is intended for use in matching V-gene segments that may include indels.

        As for match(), the function find the best match of matchseq to the specified allele and returns the
        corresponding germline fragment. It differs in that it takes advantage of knowledge that the nt sequences
        provided from complete codons. It translates the sequences first, aligns them as protein and then performs
        the reverse translation. If insertions or substitutions are found, the germline sequence is modified for best
        fit by matching the deletion or substitution in the target sequence. This mirrors the assumption that the
        indel occurred during junction formation. If more than one AA alignment has the same maximum score, the alignment
        that has the greatest number of matches at the nucleotide level is used.

        matchseq is expected to start on a codon boundary and to contain a whole number of codons. It is intended to
        represent the V-gene from position 1 to the second Cysteine (position 104 in an IMGT alignment).

        ref: Abascal F, Zardoya R, Telford MJ (2010) TranslatorX: multiple alignment of nucleotide sequences guided by
        amino acid translations
        
        match_from_aa returns two strings. The first is the best match to matchseq, as above. The second is a 'match string',
        aligned to the target, containing the characters m (match) n (nonmatch) d (deletion) i (insertion), illustrating
        how the target sequence differs to matchseq at each position.


        >>> gl = Germlib("Oryctolagus cuniculus", germline_file="testfiles/imgt_germlines.fasta")
        >>> print gl.match_from_aa("IGHV1S45*01", "caggagcaactggaggagtccgggggaggcctggtcaagcctgcgggatccctgacactcacctgcaaagcctctggattcgacctcagtggctactggtacatgtgctgggtccgccaggctccagggaagggcctggagtggatcggatgcattaatactggtgctggcgatcctgactacgcgaattgggcgaaaggccgattcaccgtctccaaagcctcgtcgaccacggtgactctgcaaatgaccagtctgacagccgcggacacggccacctatttctgtgcgaga")
        ('caggagcagctggaggagtccgggggagacctggtcaagcctgagggatccctgacactcacctgcacagcctctggattctccttcagtagctactggtacatatgctgggtccgccaggctccagggaaggggctggagtggatcgcatgcatttatgctggtagtggtagcacttactacgcgagctgggcgaaaggccgattcaccatctccaaaacctcgtcgaccacggtgactctgcaaatgaccagtctgacagccgcggacacggccacctatttctgtgcgaga', 'mmmmmmmmnmmmmmmmmmmmmmmmmmmmnmmmmmmmmmmmmmmnmmmmmmmmmmmmmmmmmmmmmmmnmmmmmmmmmmmmmnnmnmmmmmnmmdddmmmmmmiiimmnmmmmmmmmmmmmmmmmmmmmmmmmmmmmmnmmmmmmmmmmmmmnmmmmmmmnmmnmmmmmnnmdddmmnnnnnmmnmmmmmmmmmnnmmmmmmmmmmmmmmmmmmmmmnmmmmmmmmnmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm')
        """
        def chunks(l, n):
            """ Yield successive n-sized chunks from l."""
            for i in xrange(0, len(l), n):
                yield l[i:i+n]

        def map_nt_to_aa(targ_nt, match_nt, targ_res, match_res):
            """Map codons to an alignment as described in the preamble to match_from_aa()
               Incorporates code from https://www.biostars.org/p/89741/
            """
            matchstr = ""
            res_codons = []
            t_codon = iter([codon for codon in chunks(targ_nt,3)])
            m_codon = iter([codon for codon in chunks(match_nt,3)])
            for t_res, m_res in zip(targ_res, match_res):
                if t_res == '-':
                    res_codons.append(next(m_codon))  # this is an insertion during junction formation
                    matchstr += "iii"
                elif m_res == '-':
                    next(t_codon)           # this is a deletion during junction formation
                    matchstr += "ddd"
                else:
                    tn = next(t_codon)
                    mn = next(m_codon)
                    res_codons.append(tn)
                    matchstr += ''.join('m' if t == m else 'n' for t, m in zip(tn, mn))

            return ''.join(res_codons), matchstr

        if self.seq(name) is None:
            raise ValueError('Germline %s not found' % name)

        targ_nt = str(self.seq(name).seq)
        if targ_nt is None:
            raise ValueError('Germline %s not found' % name)

        # Trim any incomplete codon at the 3' end of the germline sequence. This will always be above the second Cys,
        # provided the sequence is complete.
        targ_nt = (targ_nt if len(targ_nt) % 3 == 0 else targ_nt[:0-(len(targ_nt) % 3)])

        match_nt = str(matchseq)

        if len(match_nt) %3 != 0:
            raise ValueError('Match sequence must be a whole number of codons in length')

        targ_aa = str(Seq(targ_nt, IUPAC.unambiguous_dna).translate())
        match_aa = str(Seq(match_nt, IUPAC.unambiguous_dna).translate())
        
        if "*" in match_aa:
            raise ValueError('Match sequence contains stop codon(s)')
        
        if "*" in targ_aa:
            raise ValueError('Target sequence contains stop codon(s)')
        
        alignments = pairwise2.align.localms(match_aa, targ_aa, 5, -4, -10, -10, one_alignment_only=False)

        best = None
        best_matchstr = ""
        max_match = 0
        for alignment in alignments:
            (match_res, targ_res, score, begin, end) = alignment
            res, matchstr = map_nt_to_aa(targ_nt, match_nt, targ_res, match_res)
            match = self.matches(res, targ_nt)
            if match > max_match:
                best = res
                best_matchstr = matchstr
                max_match = match

        return best, best_matchstr

    @staticmethod
    def translate_imgt_name(name):
        """Parse the allele name found in an IMGT results file (CSV or XLS) and translate it to the allele name
         found in an IMGT germline file. For the time being, as far as I know, the allele name is always the second
         word in the results file and this is assumed here.

        >>> print Germlib.translate_imgt_name("Orycun IGHV1S45*01 F")
        IGHV1S45*01

        If multiple results are shown, the first one is taken:
        >>> print Germlib.translate_imgt_name("Orycun IGHV1S26*01 F, or Orycun IGHV1S52*01 F (see comment)")
        IGHV1S26*01
        """
        return name.split(" ")[1]

    @staticmethod
    def matches(seq1, seq2):
        """Return the number of matches in two sequences of identical length"""
        return sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2))



if __name__ == "__main__":
    import doctest
    doctest.testmod()

