#!python

"""Script for analyzing epitopes in a protein.

This script is designed for analyzing epitopes in a protein, where the
epitopes are obtained from the Immune Epitope Database.

Written by Jesse Bloom."""


import os
import re
import sys
import epitopefinder.epitopes
import epitopefinder.io
import epitopefinder.sequtils
import epitopefinder.align


def main():
    """Main body of script."""
    # output is written to out, currently set to standard out
    out = sys.stdout
    out.write("Beginning execution of epitopefinder_getepitopes.py\n")
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
    iedbfile = epitopefinder.io.ParseStringValue(d, 'iedbfile').strip()
    if not os.path.isfile(iedbfile):
        raise IOError("Failed to find iedbfile of %s" % iedbfile)
    supertypesfile = epitopefinder.io.ParseStringValue(d, 'supertypesfile').strip()
    if os.path.isfile(supertypesfile):
        supertype_d = epitopefinder.epitopes.ReadSuperTypeClassification(supertypesfile, unclassified_to_none=True)
    else:
        raise IOError("Failed to find supertypesfile of %s" % supertypesfile)
    mhcclass = epitopefinder.io.ParseStringValue(d, 'mhcclass').strip().upper()
    if mhcclass not in ['I', 'II']:
        raise ValueError("Invalid mhcclass, must be I or II but found: %s" % mhcclass)
    epitopelength = epitopefinder.io.ParseStringValue(d, 'epitopelength').split()
    try:
        epitopelength = (int(epitopelength[0]), int(epitopelength[1]))
    except ValueError, TypeError:
        raise ValueError('epitopelength failed to specify two numbers giving min and max length: %s' % epitopelength)
    if epitopelength[0] > epitopelength[1]:
        raise ValueError("minimum epitope length exceeds maximum in %s" % str(epitopelength))
    targetprotsfile = epitopefinder.io.ParseStringValue(d, 'targetprotsfile').strip()
    if not os.path.isfile(targetprotsfile):
        raise IOError("Cannot find targetprotsfile %s" % targetprotsfile)
    targetprots = epitopefinder.sequtils.Read(targetprotsfile)
    seqlength = len(targetprots[0][1])
    for (head, prot) in targetprots:
        if len(prot) != seqlength:
            raise ValueError("The proteins in targetprotsfile do not appear to be aligned, as they do not all have the same length as the first sequence.")
    maxmismatches = epitopefinder.io.ParseIntValue(d, 'maxmismatches')
    musclepath = epitopefinder.io.ParseStringValue(d, 'musclepath')
    purgeredundant = epitopefinder.io.ParseStringValue(d, 'purgeredundant').strip()
    if purgeredundant == 'None':
        purgeredundant = None
    elif purgeredundant in ['MHCgene', 'MHCgroup']:
        pass
    elif purgeredundant == 'MHCsupertype':
        if mhcclass != 'I':
            raise ValueError("purgeredundant can only be set to MHCsupertype if mhcclass is set to I")
    else:
        raise ValueError("Invalid purgeredundant value of %s" % purgeredundant)
    purgeredundantoverlap = epitopefinder.io.ParseIntValue(d, 'purgeredundantoverlap')
    if not os.path.isdir(musclepath):
        raise IOError('musclepath not found: %s' % musclepath)
    epitopeslistfile = epitopefinder.io.ParseStringValue(d, 'epitopeslistfile')
    epitopesbysitefile = epitopefinder.io.ParseStringValue(d, 'epitopesbysitefile')
    # read epitopes
    out.write("\nReading epitopes from iedbfile %s...\n" % (iedbfile))
    out.flush()
    epitopes = epitopefinder.epitopes.ReadCompactIEDB(iedbfile)
    out.write("Read a total of %d epitopes.\n" % len(epitopes))
    out.write("\nNow assigning alleles to class / gene / supertype / group.")
    out.flush()
    unassigned = []
    for ep in epitopes:
        epitopefinder.epitopes.AssignAlleleInfo(ep, supertype_d)
    epitopefinder.epitopes.PrintAlleleSummary(epitopes, out)
    epitopes = [ep for ep in epitopes if ep.mhcclass == mhcclass]
    out.write("\nWe are searcing for MHC class %s epitopes. We have %d such epitopes.\n" % (mhcclass, len(epitopes)))
    epitopes = [ep for ep in epitopes if epitopelength[0] <= len(ep.sequence) <= epitopelength[1]]
    out.write("\nAfter restricting to epitopes with lengths >= %d and <= %d, we have %d epitopes remaining.\n" % (epitopelength[0], epitopelength[1], len(epitopes)))
    # align epitopes to target
    out.write('\nNow searching for matches for each epitope in the target protein sequences. These sequences are:\n')
    for (head, prot) in targetprots:
        out.write('  %s\n' % head.strip())
    matchedepitopes = []
    for ep in epitopes:
        matchlist = epitopefinder.align.GetEpitopeAlignments(ep.sequence, targetprots, maxmismatches, musclepath)
        if matchlist:
            matchedepitopes.append((ep, matchlist))
    epitopes = [tup[0] for tup in matchedepitopes]
    positions = [(tup[1][0][0], tup[1][0][1]) for tup in matchedepitopes]
    out.write("Found %d epitopes that aligned with no gaps and <= %s mismatches.\n" % (len(epitopes), maxmismatches))
    # purge redundant epitopes
    if purgeredundant:
        out.write("\nNow purging as redundant epitopes that have the same %s and >= %d residues of overlap in their alignment to the target protein." % (purgeredundant, purgeredundantoverlap))
        (unique_epitopes, unique_positions, redundant_epitopes, redundant_positions) = epitopefinder.epitopes.PurgeRedundantEpitopes(epitopes, positions, purgeredundantoverlap, purgeredundant)
        (epitopes, positions) = (unique_epitopes, unique_positions)
        out.write("\nFound a total of %d unique epitopes.\n" % len(epitopes))
    else:
        out.write("\nNot performing any purging of redundant epitopes.\n")
        redundant_epitopes = [] * len(epitopes)
        redundant_positions = [] * len(epitopes)
    assert len(epitopes) == len(positions) == len(redundant_epitopes) == len(redundant_positions)
    # write epitopes to epitopeslistfile
    # first sort by the first alignment position
    maxredundant = max([len(xlist) for xlist in redundant_epitopes]) 
    firstpositions = [(positions[i][0], i) for i in range(len(positions))]
    firstpositions.sort()
    out.write('\nWriting list of the unique and redundant epitopes to %s\n' % epitopeslistfile)
    f = open(epitopeslistfile, 'w')
    for title in ['unique_epitope'] + ['redundant-%s_epitope' % (i + 1) for i in range(maxredundant)]:
        f.write('%s_sequence,%s_position,%s_allele,%s_information,%s_reference' % (title, title, title, title, title))
        if (str(maxredundant) not in title) and (maxredundant or 'unique' not in title):
            f.write(',')
    f.write('\n')
    for (x, i) in firstpositions:
        line = []
        line.append(epitopefinder.epitopes.EpitopeSummaryString(epitopes[i], positions[i]))
        nredundant = len(redundant_epitopes[i])
        assert nredundant == len(redundant_positions[i])
        for (jep, jpos) in zip(redundant_epitopes[i], redundant_positions[i]):
            line.append(epitopefinder.epitopes.EpitopeSummaryString(jep, jpos))
        line.append(',' * 4 * (maxredundant - nredundant))
        f.write('%s\n' % (','.join(line)))
    f.close()
    # write epitopes to epitopesbysitefile
    out.write('\nWriting list of epitopes per site to %s\n' % epitopesbysitefile)
    nsites = len(targetprots[0][1])
    sitecounts = dict([(i, 0) for i in range(1, nsites + 1)])
    for (istart, iend) in positions:
        for i in range(istart, iend + 1):
            sitecounts[i] += 1
    f = open(epitopesbysitefile, 'w')
    f.write('Site,NumberUniqueEpitopes\n')
    for i in range(1, nsites + 1):
        f.write('%d,%d\n' % (i, sitecounts[i]))
    f.close()
    out.write("\nScript is complete.\n")



if __name__ == '__main__':
    main() # run the script
