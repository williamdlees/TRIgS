# PlotGermline #

PlotGermline will read one or more IMGT/IgBLASTPlus format files, and plot a histogram for each one, showing relative usage of the requested germline gene. If duplicate records have been consolidated, it can take account of a count in the ID field indicating the number of duplicates.

## Usage ##

     PlotGermline.py [-h] [-t TITLES] [-l LIMIT] [-s SAVE] [-d DUPHEADER]
                     [-c COLS] [-a] [-f] [-y YMAX]
                     infiles field detail

Argument|Meaning
---------|-------
`infiles`|Input file names, separated by commas (IMGT/IgBLASTPlus format.
`field`|The germline field to analyse, e.g. 'V-GENE and allele'
`detail`|F (family), G (germline) or A (allele) (this assumes naming similar to human germline genes: otherwise only A is supported)
`[-t TITLES]`|Comma-separated list of titles to use for plots (must have the same number of items as there are files)
`[-l LIMIT]`|Limit to the specified number of most frequently used genes
`[-s SAVE]`|Save the plot to the specified file, otherwise display it interactively
`[-d DUPHEADER]`|Prefix of the duplicate count in the Sequence ID, e.g. 'DUPCOUNT='
`[-c COLS]`|Number of plots to place on each row of the output (default 1)
`[-a]`|Sort genes alphabetically in the plot (default is by decreasing usage)
`[-f]`|Use frequency on the y-axis rather than number of reads
`[-y YMAX]`|Max y-value to use on all plots
`[-h]`|Provide a help message and exit

## Output Files ##

If the -s option is specified, the format will be determined by the extension of the specified output file(.pdf, .png, .jpg).

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data.