# Spectratype #

Spectratype will read one or more IMGT/IgBLASTPlus format files, and plot a histogram for each one, showing the relative frequency of each CDR3 length. The field IMGT-CDR3 is used, and is required to be present in the file. It may contain either amino acid or nucleotide sequences. If duplicate records have been consolidated, the script can take account of a count in the ID field indicating the number of duplicates.

## Usage ##

     Spectratype.py [-h] [-u] [-t TITLES] [-c COLS] [-d DUPHEADER] [-s SAVE]
                    [-y YMAX] [-x XMAX]
                    infiles

Argument|Meaning
---------|-------
`infiles`|Input file names, separated by commas (IMGT/IgBLASTPlus format.
`[-t TITLES]`|Comma-separated list of titles to use for plots (must have the same number of items as there are files). Default is to use filenames.
`[-c COLS]`|Number of plots to place on each row of the output (default 1)
`[-d DUPHEADER]`|Prefix of the duplicate count in the Sequence ID, e.g. 'DUPCOUNT='
`[-s SAVE]`|Save the plot to the specified file, otherwise display it interactively
`[-y YMAX]`|Max y-value to use on all plots
`[-x XMAX]`|Max x-value to use on all plots
`[-u]`|Only count unique sequences
`[-h]`|Provide a help message and exit

## Output Files ##

If the -s option is specified, the format will be determined by the extension of the specified output file(.pdf, .png, .jpg).

## Testing ##

`testfiles/Test_Commands.bat` exercises the script with sample data.