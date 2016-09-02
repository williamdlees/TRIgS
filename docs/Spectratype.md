# Spectratype #

Spectratype will read one or more IMGT/IgBLASTPlus format files, and plot a histogram for each one, showing the relative frequency of each CDR3 length. The field IMGT-CDR3 is used, and is required to be present in the file. It may contain either amino acid or nucleotide sequences. If duplicate records have been consolidated, the script can take account of a count in the ID field indicating the number of duplicates.

SpectraType requires [NumPy](http://www.numpy.org) and [matplotlib](http://matplotlib.org).


## Usage ##

    Spectratype.py [-h] [-b BARCOLOUR] [-c COLS] [-d DUPHEADER] [-g] [-gh]
                   [-gv GRID_VERTICAL] [-s SAVE] [-t TITLES] [-u]
                   [-w WIDTH] [-xmax XMAX] [-xmin XMIN] [-y YMAX]
                   infiles

Argument|Meaning
---------|-------
`infiles`|Input file names, separated by commas (IMGT/IgBLASTPlus format.
`[-b BARCOLOUR]`|Colours to use for plots (one or more colours, separated by commas). If there are more plots than colours, the colours will repeat in a cycle
`[-c COLS]`|Number of plots to place on each row of the output (default 1)
`[-d DUPHEADER]`|Prefix of the duplicate count in the Sequence ID, e.g. 'DUPCOUNT='
`[-g]`|Fill bars with a gradiented colour (as in the example above)
`[-gh]`|Show horizontal grid lines
`[-gv GRID_VERTICAL]`|Show vertical grid lines after every GRID_VERTICAL number of bars
`[-s SAVE]`|Save the plot to the specified file, otherwise display it interactively
`[-t TITLES]`|Comma-separated list of titles to use for plots (must have the same number of items as there are files). Default is to use filenames.
`[-u]`|Only count unique sequences
`[-w WIDTH]`|Column width, between 0 and 1 (0.8 in the above example)
`[-y YMAX]`|Max y-value (read count) to use on all plots
`[-xmax XMAX]`|Max x-value (CDR3 length) to use on all plots
`[-xmin XMAX]`|Max x-value (CDR3 length) to use on all plots
`[-h]`|Provide a help message and exit

## Output Files ##

If the -s option is specified, the format will be determined by the extension of the specified output file(.pdf, .png, .jpg). If the extension is .csv, a csv file will be produced instead of a plot.

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data.