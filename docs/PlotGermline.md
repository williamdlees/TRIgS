# PlotGermline #

PlotGermline will read one or more IMGT/IgBLASTPlus format files, and plot a histogram for each one, showing relative usage of the requested germline gene. If duplicate records have been consolidated, it can take account of a count in the ID field indicating the number of duplicates.

PlotGermline requires [NumPy](http://www.numpy.org) and [matplotlib](http://matplotlib.org).

<img src="https://rawgit.com/williamdlees/BioTools/master/docs/germlines.png" width="800">

## Usage ##

    usage: PlotGermline.py [-h] [-a] [-b BARCOLOUR] [-c COLS] [-co] [-d DUPHEADER] [-f]
                           [-g] [-gh] [-gv GRID_VERTICAL] [-l LIMIT] [-s SAVE]
                           [-sz X,Y] [-t TITLES] [-w WIDTH] [-y YMAX]
                           infiles field detail

Argument|Meaning
---------|-------
`infiles`|Input file names, separated by commas (IMGT/IgBLASTPlus format.
`field`|The germline field to analyse, e.g. 'V-GENE and allele'
`detail`|F (family), G (germline) or A (allele) (this assumes naming similar to human germline genes: otherwise only A is supported)
`[-a]`|Sort genes alphabetically in the plot (default is by decreasing usage)
`[-b BARCOLOUR]`|Colours to use for plots (one or more colours, separated by commas). If there are more plots than colours, the colours will repeat in a cycle.
`[-c COLS]`|Number of plots to place on each row of the output (default 1)
`[-co]`|Consolidate statistics into a single plot or table
`[-d DUPHEADER]`|Prefix of the duplicate count in the Sequence ID, e.g. 'DUPCOUNT='
`[-f]`|Use frequency on the y-axis rather than number of reads
`[-g]`|Fill bars with a gradiented colour (as in the example above)
`[-gh]`|Show horizontal grid lines
`[-gv GRID_VERTICAL]`|Show vertical grid lines after every GRID_VERTICAL number of bars
`[-l LIMIT]`|Limit to the specified number of most frequently used genes
`[-s SAVE]`|Save the plot to the specified file, otherwise display it interactively
`[-sz X,Y]`|Set the total figure size to X by Y inches (note that the smaller the overall figure is, the larger the text is with respect to the plot, which gives some control over font size)
`[-t TITLES]`|Comma-separated list of titles to use for plots (if there are fewer titles than plots, the titles will repeat in cycle)
`[-w WIDTH]`|Column width, between 0 and 1 (0.8 in the above example)
`[-y YMAX]`|Max y-value to use on all plots
`[-h]`|Provide a help message and exit

## Output Files ##

If the -s option is specified, the format will be determined by the extension of the specified output file(.pdf, .png, .jpg). If the extension is .csv, a csv file will be produced instead of a plot.

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data.