# PlotIdentity #

Identity/Divergence plots show, for a set of sampled sequences, their divergence from the inferred germline and identity to a target sequence of interest:

<img src="https://rawgit.com/williamdlees/BioTools/master/docs/identity.png" width="200">
<br>Example Identity/Divergence plot, produced from the analysis of an NGS-based heavy-chain repertoire.


## Usage ##

    PlotIdentity.py [-h] [-a] [-b] [-c COLOURMAP] [-g BACKGROUND]
                    [-mx MAXX] [-my MINY] [-p POINTS] [-s SAVE]
                    repertoire

Argument|Meaning
---------|-------
`repertoire`|File containing the repertoire identities and divergences, in the format produced by AbIdentity
`[-a]`|Attempt to arrange sequence labels so that they do not overlap. This option requires the Python package adjustText to be installed.
`[-b]`|Plot a key (colour bar) showing colours and associated sequence density.
`[-c colourmap`|Use the specified colour map. The list of supported options is [here](http://matplotlib.org/examples/color/colormaps_reference.html "here"). Note that in most cases the map order can be reversed by appending _r to the name.
`[-g background]|Specifies the colour to use where the sequence density is zero (this overrides the colourmap). Available options are [http://matplotlib.org/api/colors_api.html](http://matplotlib.org/api/colors_api.html "here")
`[-mx MAXX]`|The maximum germline divergence to be plotted (sequences with higher divergences will not be shown). If no value is specified, the plot will extend as needed to cover all sequences.
`[-my MINY]`|The minimum Ab identity to be plotted (sequences with lower identity will not be plotted). If no value is specified, the plot will extend as needed to cover all sequences.'
`[-p POINTS]`|Additional files containing individual sequences to be plotted as distinct points (see notes below for the format of this option)
`[-s SAVE]`|Output file (.png, pdf, .jpg). The format will be determined by the extension.
`[-h]`|Provide a help message and exit

## Notes ##

The repertoire dataset will typically include all sequences matching a given germline. Because of its size, it is shown as a contour plot, denoting the density of sequences at each point in the plot. Specific sequences of interest, such as that of an mAb of interest, its inferred germline sequence, and possibly additional sequences such as intermediates, are shown as distinct points. The -p option is used to enable multiple sets od such sequences to be rendered, with control over format.

The format of the option is:

-p "file1,specifier1,file2,specifier2,...filen,specifien"

At least one file and specifier must be provided.

The file(s) each contain one or more sequence IDs, together with their identity and divergence values, in the format produced by AbIdentity.py. The format specifier has three sections, of which the first is mandatory and the others are optional. The three sections are divided by the character /.

Section|Format
-------|------
Display|Defines the colour and shape of the mark to be displayed, using Matplotlib format (see documentation for [ax.plot()](http://matplotlib.org/api/axes_api.html "ax.plot()") for a full description). 
MarkSize|Consists of the letter m followed by an integer. If present, defines the size of the point to be displayed (in pixels).
Label|Consists of the letter l followed by a colour specifier, in Matplotlib format (see documentation for [ax.plot()](http://matplotlib.org/api/axes_api.html "ax.plot()") for a full description). If present, the SequenceID will be displayed beside the point, in the colour specified. 

Example Specifiers

Specifier|Meaning
---------|-------
ro|A red circle of default size
b+/m20|A blue cross of width 20 pixels
yv/lr|A yellow triangle, with the SequenceID displayed as a label, in red
gDv/m20/lr|A large green diamond, with the SequenceID displayed as a label, in red


## Example ##

`testfiles/Test_Commands.bat` contains the command that generates the sample image at the head of this section.