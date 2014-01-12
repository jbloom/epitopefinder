#!python

"""Script for plotting correlation between two properties of protein sites.

Uses matplotlib. Designed to analyze output of epitopefinder_getepitopes.py.

Written by Jesse Bloom."""


import os
import sys
import epitopefinder.io
import epitopefinder.plot


def main():
    """Main body of script."""
    if not epitopefinder.plot.PylabAvailable():
        raise ImportError("Cannot import matplotlib / pylab, which are required by this script.")
    # output is written to out, currently set to standard out
    out = sys.stdout
    out.write("Beginning execution of epitopefinder_plotcorrelation.py\n")
    # read input file and parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    d = epitopefinder.io.ParseInfile(open(infilename))
    out.write("\nRead input arguments from %s\n" % infilename)
    out.write('Read the following key / value pairs:\n')
    for (key, value) in d.iteritems():
        out.write("%s %s\n" % (key, value))
    plotfile = epitopefinder.io.ParseStringValue(d, 'plotfile').strip()
    xdatafile = epitopefinder.io.ParseStringValue(d, 'xdatafile').strip()
    ydatafile = epitopefinder.io.ParseStringValue(d, 'ydatafile').strip()
    for f in [xdatafile, ydatafile]:
        if not os.path.isfile(f):
            raise IOError("Failed to find file %s" % f)
    title = epitopefinder.io.ParseStringValue(d, 'title').strip()
    if title == 'False':
        title = False
    xlabel = epitopefinder.io.ParseStringValue(d, 'xlabel').strip()
    ylabel = epitopefinder.io.ParseStringValue(d, 'ylabel').strip()
    correlation = epitopefinder.io.ParseStringValue(d, 'correlation').strip()
    if correlation == 'None':
        correlation = None
    elif correlation not in ['Pearson', 'Spearman']:
        raise ValueError('correlation must be Pearson or Spearman, instead read %s' % correlation)
    xs = []
    ys = []
    for (data, fname) in [(xs, xdatafile), (ys, ydatafile)]:
        for line in open(fname).readlines()[1 : ]:
            if line.isspace() or line[0] == '#':
                continue
            entries = line.split(',')
            if len(entries) != 2:
                raise ValueError("Line in %s does not contain two entries:\n%s" % fname, line)
            (r, z) = (float(entries[0]), float(entries[1]))
            data.append((r, z))
    xshared = []
    yshared = []
    ydict = dict(ys)
    for (r, x) in xs:
        if r in ydict:
           xshared.append(x)
           yshared.append(ydict[r])
    if len(xshared) < 2:
        raise ValueError("data files %s and %s contain less than %d shared data points (properties for the same site)" % (xdatafile, ydatafile))
    out.write('\nNow plotting the %d data points in the file %s\n' % (len(xshared), plotfile))
    epitopefinder.plot.CorrelationPlot(xshared, yshared, plotfile, xlabel, ylabel, corr=correlation, title=title)
    out.write("\nScript is complete.\n")


if __name__ == '__main__':
    main() # run the script
