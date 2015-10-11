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


from Bio import AlignIO, SeqIO
from Bio import Phylo
from StringIO import StringIO
from Alignment import Alignment
from AnalyseCDR import AnalyseCDR
import os
import subprocess
import copy
import re
import sys

class Dnaml:
    """A class to manage preparation of files for dnaml, invocation, and parsing of results. At the moment this is
    restricted to the use of dnaml for ancestral reconstruction, based on an existing tree."""

    def __init__(self):
        pass

    def run_dnaml(self, seq_align, ptree, seqpattern, cdrfile, wdir, rep, tag=""):
        """Run dnaml. Arguments are:
           seq_align: the input nt sequences (MultipleSequenceAlignment)
           ptree: phylogenetic tree (Bio.Phylo)
           seqpattern: A list of sequence number directives, in the format accepted by Alignment.set_position_numbers
           wdir: the name of a directory that run_paml should use. This must exist already.
           rep: a function that takes a string as an argument. This will be called should an error or warning be
                necessary (may be called multiple times in one invocation).
           tag: an optional tag to prefix filenames with

           Sequences in seq_align must be the same length, must start on a codon boundary, and be an integral number
           of codons in length. The first sequence must be the ancestral sequence or outgroup. Exactly he same sequence
           names must occur in the alignment and the tree. Sequence name format is pretty flexible (sequences are
           mapped to names acceptable to PAML and remapped after PAML has run).
        """
        root_id = seq_align[0].id

        # Translate clade names to something safe
        namedict = {}
        serial = 1

        for seq in seq_align:
            namedict[seq.id] = "N%09d" % serial
            seq.id = namedict[seq.id]
            serial += 1

        qtree = copy.deepcopy(ptree)

        for clade in qtree.get_terminals():
            if clade.name and clade.name in namedict:
                clade.name = namedict[clade.name]

        # Root the tree on the first record

        first = "N%09d" % 1
        
        try:
            qtree.root_with_outgroup(qtree.find_clades(name=re.escape(first)).next())
        except:
            raise ValueError("Error: root sequence not found in tree.")
                    
        try:
            inv_dict = {v: k for k, v in namedict.items()}
                
            ptree.root_with_outgroup(ptree.find_clades(name=re.escape(inv_dict[first])))
            Phylo.write(ptree, wdir + "/" + "input_treefile.new", "newick", plain=False)
        except:
            raise ValueError("Error rooting trees: check for corrupt tree file or duplicated sequences.")    

        # Write the sequences, in PHYLIP format (real PHYLIP format, as used by PHYLIP!)

        with open(wdir + "/" + "infile", "w") as f:
            f.write("  %d  %d\n" % (len(seq_align), len(seq_align[0])))
            for seq in seq_align:
                f.write("%10s%s\n" % (seq.id, seq.seq.upper()))

        # Write the tree file

        Phylo.write(qtree, wdir + "/" + "intree", "newick")

        if os.path.exists(wdir + "/" + "outfile"):
            os.remove(wdir + "/" + "outfile")
        if os.path.exists(wdir + "/" + "outtree"):
            os.remove(wdir + "/" + "outtree")

        # The path to the module may reference either a .py or a .pyc file...
        
        ctlfile = os.path.abspath(__file__).replace(".pyc", ".ctl") if ".pyc" in os.path.abspath(__file__) \
            else os.path.abspath(__file__).replace(".py", ".ctl")
        
        # Check for dnaml in the current directory
        
        dnamlfile = os.path.abspath(__file__).replace("Dnaml.pyc", "dnaml") if ".pyc" in os.path.abspath(__file__) \
            else os.path.abspath(__file__).replace("Dnaml.py", "dnaml")
        
        if not os.path.exists(dnamlfile):
            dnamlfile = "dnaml" # must be on the path somewhere
        
        with open(wdir + "/" + "dnaml.txt", "w") as o, open(ctlfile, "r") as i:
            subprocess.call(dnamlfile, cwd=wdir, stdin = i, stdout=o)

        if not os.path.isfile(wdir + "/" + "outfile"):
            rep("No output returned by dnaml: please check the logs for the issue.")
            return None

        if os.path.isfile(wdir + "/" + "outfile.txt"):
            os.remove(wdir + "/" + "outfile.txt")
        os.rename(wdir + "/" + "outfile", wdir + "/" + "outfile.txt")

        intseqs = self.__parse_outfile(wdir + "/" + "outfile.txt")

        if not intseqs:
            rep("Unexpected output returned by dnaml: please check the logs for the issue.")
            return None

        # Custom sort function to put the root record first, then others supplied by the user, then intermediate nodes
        def key_ids(rec):
            if rec.id == "N%09d" % 1:
                return 'a__' + rec.id
            elif 'node_' in rec.id:
                return 'z__' + "%04d" % (int)(rec.id.split("_")[1])
            else:
                return 'l__' + rec.id

        labelled_tree = Phylo.read(wdir + "/" + "outtree", "newick")
        intseqs.seek(0)
        int_seqs = Alignment(file_name=intseqs, format="fasta")
        int_seqs.sort(key=key_ids)
        intseqs.seek(0)
        int_aas = Alignment()
        int_aas.read_nt(intseqs, "fasta")
        int_aas.sort(key=key_ids)
        int_aas.set_position_numbers(position_numbers = seqpattern)

        # Put back the original names in all our collections

        for seq in int_seqs:
            if seq.id in inv_dict:
                seq.id = inv_dict[seq.id]
            seq.name = ""
            seq.description = ""

        for seq in int_aas:
            if seq.id in inv_dict:
                seq.id = inv_dict[seq.id]
            seq.name = ""
            seq.description = ""

        nodeid = 1
        for clade in labelled_tree.find_clades(order="preorder"):
            if clade.name is None:
                clade.name = "node_%d" % nodeid            # This relies on our traversal using the same order as dnaml
                nodeid += 1
            else:
                if clade.name in inv_dict:
                    clade.name = inv_dict[clade.name]

        # Now we need to map the labelling of the nodes in the labelled tree to the nodes in the original tree

        self.__map_names(ptree, labelled_tree)
        Phylo.write(ptree, wdir + "/" + tag + "intermediates_treefile.new", "newick", plain=False)

        cladenames = []
        new_int_aas = Alignment()

        for clade in ptree.find_clades():
            if clade.name is not None:
                cladenames.append(clade.name)

        for rec in int_aas:
            if rec.id in cladenames:
                new_int_aas.append(rec)

        int_aas = new_int_aas
        int_aas.set_position_numbers(position_numbers = seqpattern)

        copy_tree = copy.deepcopy(ptree)
        # Calculate AA diffs between each node and its parent, and write to the tree

        labels = {}

        def diffkey(diff):
            return int_aas.index_of(diff[1:-1])

        for clade in ptree.find_clades():
            if clade.name is not None:
                parent = self.__get_parent(ptree, clade)

                if parent is None:
                    path = ptree.get_path(clade)
                    if len(path) == 1 and clade.name != first:
                        fname = inv_dict[first]
                        parent = ptree.find_clades(name = re.escape(fname)).next()

                if parent is not None and parent.name is not None:
                    diffs = list(int_aas.seqdiff(clade.name, parent.name))
                    diffs.sort(key = diffkey)
                    diffs = "+".join(diffs)
                    if "node_" in clade.name:
                        labels[clade.name] = diffs
                    else:
                        labels[clade.name] = str(clade.name) + " " + diffs

        for clade in ptree.find_clades():
            if clade.name is not None and clade.name in labels:
                clade.name = labels[clade.name]

        Phylo.write(ptree, wdir + "/" + tag + "annotated_treefile.new", "newick", plain=False)

        # Now write a tree with summary CDR/FR total changes

        if cdrfile is not None:
            ptree = copy.deepcopy(copy_tree)
            acdr = AnalyseCDR(int_aas, file_name=cdrfile)
            labels = {}
    
            for clade in ptree.find_clades():
                if clade.name is not None:
                    parent = self.__get_parent(ptree, clade)
    
                    if parent is None:
                        path = ptree.get_path(clade)
                        if len(path) == 1 and clade.name != first:
                            fname = inv_dict[first]
                            parent = ptree.find_clades(name = re.escape(fname)).next()
    
                    if parent is not None and parent.name is not None:
                        diffs = acdr.category_diff(clade.name, parent.name)
                        if "node_" in clade.name:
                            labels[clade.name] = diffs
                        else:
                            labels[clade.name] = str(clade.name) + " " + diffs
    
            for clade in ptree.find_clades():
                if clade.name is not None and clade.name in labels:
                    clade.name = labels[clade.name]
    
            Phylo.write(ptree, wdir + "/" + tag + "annotated_treefile_sum.new", "newick", plain=False)

        # And write a tree with counts of total AA changes

        ptree = copy.deepcopy(copy_tree)
        labels = {}

        for clade in ptree.find_clades():
            if clade.name is not None:
                parent = self.__get_parent(ptree, clade)

                if parent is None:
                    path = ptree.get_path(clade)
                    if len(path) == 1 and clade.name != first:
                        fname = inv_dict[first]
                        parent = ptree.find_clades(name = re.escape(fname)).next()

                if parent is not None and parent.name is not None:
                    diffs = list(int_aas.seqdiff(clade.name, parent.name))
                    if "node_" in clade.name:
                        labels[clade.name] = str(len(diffs)) if len(diffs) > 0 else ""
                    else:
                        labels[clade.name] = str(clade.name) + (" " + str(len(diffs)) if len(diffs) > 0 else "")

        for clade in ptree.find_clades():
            if clade.name is not None and clade.name in labels:
                clade.name = labels[clade.name]

        Phylo.write(ptree, wdir + "/" + tag + "annotated_treefile_tot.new", "newick", plain=False)

        f = open(wdir + "/" + tag + "aa_alignment.txt", "w")
        f.write(int_aas.report(100))
        f.close()
        
        f = open(wdir + "/" + tag + "nt_alignment.txt", "w")
        f.write(int_seqs.report(100))
        f.close()

        for rec in int_aas:
            rec.description = ""

        AlignIO.write(int_aas, wdir + "/" + tag + "aa_alignment.fa", "fasta")
        AlignIO.write(int_seqs, wdir + "/" + tag + "nt_alignment.fa", "fasta")
        return int_aas

    def __parse_outfile(self, filename):
        """Internal method to parse the dnaml output file created after ancestral reconstruction."""
        #Fish out the tree with node labels, and the ancestral sequences

        seqs = {}

        with open(filename, "r") as f:
            for line in f:
                if "Reconstructed sequence" in line or not line:
                    break

            if not line:
                return

            for line in f:
                if len(line) > 10:
                    id = line[:11].replace(" ", "")
                    if 'N' not in id:
                        id = "node_" + id
                    seq = line[11:].strip().replace(" ", "")
                    seqs[id] = seqs.get(id, "") + seq

        intseqs = StringIO()
        for id,seq in seqs.iteritems():
            intseqs.write(">%s\n" % id)
            intseqs.write("%s\n" % seq)

        return intseqs


    def __get_parent(self, tree, child_clade):
        """Internal method to find the parent of a clade"""
        node_path = tree.get_path(child_clade)
        if len(node_path) > 1:
            return node_path[-2]
        else:
            return None

    def __map_names(self, ptree, reftree):
        """Map the names of intermediate nodes across from reftree to ptree"""
        for clade in ptree.find_clades(order = 'postorder'):
            if clade.name is None:
                childname = clade.clades[0].name
                if childname is not None:
                    refchild = reftree.find_clades(name=re.escape(childname))
                    refp = self.__get_parent(reftree, refchild.next())
                    if refp is not None:
                        clade.name = refp.name
                    elif clade != ptree.clade:
                        clade.name = reftree.root.name
