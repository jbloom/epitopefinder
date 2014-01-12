#!python

"""Script for plotting linear density of epitopes as a function of primary sequence.

Uses matplotlib. Designed to analyze output of epitopefinder_getepitopes.py.

Written by Jesse Bloom."""


import os
import sys
import epitopefinder.plot


def main():
    """Main body of script."""
    if not epitopefinder.plot.PylabAvailable():
        raise ImportError("Cannot import matplotlib / pylab, which are required by this script.")
    # output is written to out, currently set to standard out
    out = sys.stdout
    out.write("Beginning execution of epitopefinder_plotlineardensity.py\n")
    # read input file and parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    lines = [line.strip() for line in open(infilename).readlines() if not line.isspace() and line[0] != '#']
    if len(lines) < 4:
        raise IOError("Failed to read enough lines from %s" % infilename)
    entries = lines[0].split(None, 1)
    if not (len(entries) == 2 and entries[0] == 'plotfile'):
        raise ValueError("First line does not validly specify plotfile:\n%s" % lines[0])
    plotfile = entries[1].strip()
    entries = lines[1].split(None, 1)
    if not (len(entries) == 2 and entries[0].strip() == 'title'):
        raise ValueError("Third line does not validly specify title:\n%s" % lines[1])
    title = entries[1].strip()
    if title.upper() == 'FALSE':
        title = False
    entries = lines[2].split()
    if not (len(entries) == 2 and entries[0].strip() == 'fixymax'):
        raise ValueError("Second line does not validly specify fixymax:\n%s" % lines[2])
    fixymax = entries[1].strip()
    if fixymax.upper() == 'FALSE':
        fixymax = False
    else:
        fixymax = float(fixymax)
    datalist = []
    for line in lines[3 : ]:
        entries = line.split(None, 1)
        if len(entries) != 2:
            raise ValueError("Line does not validly specify data file with data label:\n%s" % lines[0])
        (datafile, datalabel) = (entries[0].strip(), entries[1].strip())
        out.write('\nReading data from data file %s' % datafile)
        if not os.path.isfile(datafile):
            raise IOError("Cannot find datafile %s" % datafile)
        if os.path.splitext(datafile)[1].upper() != '.CSV':
            raise ValueError("datafile name %s does not end in .csv" % datafile)
        data = []
        for line in open(datafile).readlines()[1 : ]:
            if line[0] == '#':
                continue
            entries = line.split(',')
            if len(entries) != 2:
                raise ValueError("Line in %s does not contain two entries\n%s" % (datafile, line))
            (x, y) = (float(entries[0]), float(entries[1]))
            data.append((x, y))
        datalist.append((datalabel, data))
    out.write('\n\nNow plotting this data in %s\n' % plotfile)
    epitopefinder.plot.PlotLinearDensity(datalist, plotfile, 'residue number', 'number of epitopes', title=title, fixymax=fixymax)
    out.write("\nScript is complete.\n")


if __name__ == '__main__':
    main() # run the script
