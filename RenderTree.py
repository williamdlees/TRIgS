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

# Various methods to render trees for the website

__author__ = 'William Lees'
__docformat__ = "restructuredtext en"

from ete2 import Tree, TreeStyle, TextFace, NodeStyle

class RenderTree():
    def __init__(self):
        pass

    @staticmethod
    def render_annotate(newick_path, output_path):
        """Render the annotated tree, showing internal node names.
           The output_path should end in .PNG, .PDF or .SVG: this will determine the format.

           To aid in testing, if output_path is None, the tree is shown rather than rendered.
        """
        tree = Tree(newick_path, format=1)

        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.branch_vertical_margin = 15

        ns = NodeStyle()
        ns["size"] = 1

        edge = 0
        for node in tree.traverse():
            node.name = node.name.replace("'", "")
            node.name = node.name.replace("+", ",")
            if not node.is_leaf() and node.name != 'NoName':
                f = TextFace(node.name)
                f.margin_top = 5
                f.margin_bottom = 5
                f.margin_right = 10
                f.margin_left = 10
                edge += 1
                node.add_face(f, column=0, position="branch-top")

        if output_path is None:
            tree.show(tree_style=ts)
        else:
            tree.render(output_path)

if __name__ == '__main__':
    RenderTree.render_annotate("C:\\Users\\william\\Dropbox (Personal)\\Research\\TreeAnnotateApp\\temp\\tmpviamgd\\annotated_treefile.new", None)
