#!python

"""Script for analyzing epitope density at selected sites.

Designed to analyze output of epitopefinder_getepitopes.py.

Written by Jesse Bloom."""


import os
import sys
import epitopefinder.io


def main():
    """Main body of script."""
    # output is written to out, currently set to standard out
    out = sys.stdout
    out.write("Beginning execution of epitopefinder_selectsites.py\n")
    # read input file and parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    d = epitopefinder.io.ParseInfile(open(infilename))
    epitopesbysitefile = epitopefinder.io.ParseStringValue(d, 'epitopesbysitefile')
    if not os.path.isfile(epitopesbysitefile):
        raise IOError("Cannot find epitopesbysitefile %s" % epitopesbysitefile)
    selectsitesfile = epitopefinder.io.ParseStringValue(d, 'selectsitesfile')
    sites = epitopefinder.io.ParseStringValue(d, 'sites')
    sites = [int(site) for site in sites.split()]
    retainmultiple = epitopefinder.io.ParseBoolValue(d, 'retainmultiple')
    if not retainmultiple:
        sites = dict([(site, True) for site in sites])
        sites = sites.keys()
    out.write('\nWill look for epitopes at the following sites: %s\n' % ', '.join([str(site) for site in sites]))
    out.write('\nThese epitopes will be written to %s.\n' % selectsitesfile)
    epitopesbysite = dict([(int(line.split(',')[0]), int(line.split(',')[1])) for line in open(epitopesbysitefile).readlines()[1 : ]])
    counts = [(epitopesbysite[site], site) for site in sites]
    counts.sort()
    counts.reverse()
    f = open(selectsitesfile, 'w')
    f.write('Site,NumberUniqueEpitopes\n')
    for (count, site) in counts:
        f.write("%d,%d\n" % (site, count))
    f.close()
    out.write("\nScript is complete.\n")


if __name__ == '__main__':
    main() # run the script
