#!python

"""Script for plotting distributions of epitopes per site for two sets of sites.

Uses matplotlib. Designed to analyze output of epitopefinder_getepitopes.py.

Written by Jesse Bloom."""


import os
import sys
import random
import epitopefinder.io
import epitopefinder.plot


def main():
    """Main body of script."""
    random.seed(1) # seed random number generator in case P values are being computed
    if not epitopefinder.plot.PylabAvailable():
        raise ImportError("Cannot import matplotlib / pylab, which are required by this script.")
    # output is written to out, currently set to standard out
    out = sys.stdout
    out.write("Beginning execution of epitopefinder_plotdistributioncomparison.py\n")
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
    epitopesbysite1_list = []
    epitopesbysite2_list = []
    for (xlist, xf) in [(epitopesbysite1_list, 'epitopesfile1'), (epitopesbysite2_list, 'epitopesfile2')]:
        epitopesfile = epitopefinder.io.ParseFileList(d, xf)
        if len(epitopesfile) != 1:
            raise ValueError("%s specifies more than one file" % xf)
        epitopesfile = epitopesfile[0]
        for line in open(epitopesfile).readlines()[1 : ]:
            if not (line.isspace() or line[0] == '#'):
                (site, n) = line.split(',')
                (site, n) = (int(site), int(n))
                xlist.append(n)
        if not xlist:
            raise ValueError("%s failed to specify information for any sites" % xf)
    set1name = epitopefinder.io.ParseStringValue(d, 'set1name')
    set2name = epitopefinder.io.ParseStringValue(d, 'set2name')
    title = epitopefinder.io.ParseStringValue(d, 'title').strip()
    if title.upper() in ['NONE', 'FALSE']:
        title = None
    pvalue = epitopefinder.io.ParseStringValue(d, 'pvalue')
    if pvalue.upper() in ['NONE', 'FALSE']:
        pvalue = None
        pvaluewithreplacement = None
    else:
        pvalue = int(pvalue)
        pvaluewithreplacement = epitopefinder.io.ParseBoolValue(d, 'pvaluewithreplacement')
        if pvalue < 1:
            raise ValueError("pvalue must be >= 1")
        if len(epitopesbysite2_list) >= len(epitopesbysite1_list):
            raise ValueError("You cannot use pvalue since epitopesbysite2_list is not a subset of epitopesbysite1_list -- it does not contain fewer sites with specified epitope counts.")
    ymax = None
    if 'ymax' in d:
        ymax = epitopefinder.io.ParseFloatValue(d, 'ymax')
    out.write('\nNow creating the plot file %s\n' % plotfile)
    epitopefinder.plot.PlotDistributionComparison(epitopesbysite1_list, epitopesbysite2_list, set1name, set2name, plotfile, 'number of epitopes', 'fraction of sites', title, pvalue, pvaluewithreplacement, ymax=ymax)
    out.write("\nScript is complete.\n")


if __name__ == '__main__':
    main() # run the script
