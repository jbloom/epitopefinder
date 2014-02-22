"""Runs the analysis to identify CTL epitopes.

Written by Jesse Bloom."""


import os
import random
import scipy.stats


def ProbDiffer(fullset1, subset1, fullset2, subset2, nrandom=100000):
    """Computes P-value that values in one subset exceed another.

    *fullset1*, *subset1*, *fullset2*, and *subset2* are all lists
    of numbers, with *subset1* being a subset of *fullset1*
    and *subset2* being a subset of *fullset2*. 

    This function draws *nrandom* random paired subsets of equal size to 
    *subset1* from *fullset1* and of equal size to *subset2* from
    *fullset2*, and returns the fraction of such draws in which
    the mean of the random *subset1* exceeds the mean of the
    random *subset2* by at least as much as thea ctual subsets.
    """
    n1 = len(subset1)
    n2 = len(subset2)
    actual_diff = sum(subset1) / float(n1) - sum(subset2) / float(n2)
    nge = 0
    for irandom in range(nrandom):
        isubset1 = random.sample(fullset1, n1)
        isubset2 = random.sample(fullset2, n2)
        i_diff = sum(isubset1) / float(n1) - sum(isubset2) / float(n2)
        if i_diff >= actual_diff:
            nge += 1
    return nge / float(nrandom)


def main():
    """Main body of script."""

    # input files and settings
    random.seed(1)
    epitopefile = 'IEDB_Influenza_Tcell_compact_2013-12-27.csv' # contains epitopes
    supertypefile = 'supertype_classification.txt' # contains supertype classifications
    epitopelength = (8, 12) # restrict epitopes to these lengths
    maxmismatches = 1 # maximum number of allowed mismatches with epitope
    musclepath = '/Users/jbloom/muscle3.8/' # find MUSCLE here
    # sequencesets is a dictionary keyed by the sequence set name.
    # Entries are dictionaries which are initially keyed by
    # "targetprotsfile" and "substitutedsites" (and in some cases "epistaticsites"), 
    # but other keys are added as the script proceeds.
    sequencesets = {
            'human_H3N2':{'targetprotsfile':'targetprots_humanH3N2.fasta', 'substitutedsites':'substitutedsites_humanH3N2.txt', 'epistaticsites':'epistaticsites_humanH3N2.txt'}, 
            'swine':{'targetprotsfile':'targetprots_swine.fasta', 'substitutedsites':'substitutedsites_swine.txt'}, 
            }

    # find the epitopes
    script = 'epitopefinder_getepitopes.py'
    for sequenceset in sequencesets:
        print "\nFinding epitopes for %s..." % sequenceset
        sequencesets[sequenceset]['epitopeslistfile'] = 'epitopeslist_%s.csv' % sequenceset
        sequencesets[sequenceset]['epitopesbysitefile'] = 'epitopesbysite_%s.csv' % sequenceset
        infile = '%s_%s_infile.txt' % (os.path.splitext(script)[0], sequenceset)
        open(infile, 'w').write('\n'.join([\
                "# input file for %s" % script,
                "iedbfile %s" % epitopefile,
                "supertypesfile %s" % supertypefile,
                "mhcclass I",
                "epitopelength %d %d" % epitopelength,
                "targetprotsfile %s" % sequencesets[sequenceset]['targetprotsfile'],
                "maxmismatches %d" % maxmismatches,
                "musclepath %s" % musclepath,
                "purgeredundant MHCgroup",
                "purgeredundantoverlap 8",
                "epitopeslistfile %s" % sequencesets[sequenceset]['epitopeslistfile'],
                "epitopesbysitefile %s" % sequencesets[sequenceset]['epitopesbysitefile'],
            ]))
        os.system("%s %s" % (script, infile))


    # get the epitope counts in the substituted sites and epistatic sites (where applicable)
    script = 'epitopefinder_selectsites.py'
    for sequenceset in sequencesets:
        print "\nGetting epitope counts for selected sites for %s..." % sequenceset
        for subsettype in ['substitutedsites', 'epistaticsites']:
            if subsettype not in sequencesets[sequenceset]:
                continue
            sequencesets[sequenceset]['epitopecounts_%s' % subsettype] = 'epitopecounts_%s_%s.csv' % (subsettype, sequenceset)
            sites = ' '.join([line.strip() for line in open(sequencesets[sequenceset][subsettype]) if not line.isspace()])
            infile = '%s_%s_%s_infile.txt' % (os.path.splitext(script)[0], sequenceset, subsettype)
            open(infile, 'w').write('\n'.join([\
                    '# input file for %s' % script,
                    'epitopesbysitefile %s' % sequencesets[sequenceset]['epitopesbysitefile'],
                    'selectsitesfile %s' % sequencesets[sequenceset]['epitopecounts_%s' % subsettype],
                    'sites %s' % sites,
                    'retainmultiple False',
                    ]))
            os.system('%s %s' % (script, infile))


    # plot the linear densities of epitopes
    script = 'epitopefinder_plotlineardensity.py'
    print "\nPlotting linear epitope density..."
    infile = '%s_infile.txt' % (os.path.splitext(script)[0])
    open(infile, 'w').write('\n'.join([\
                '# input file for %s' % script,
                'plotfile epitopelineardensity.pdf',
                'title False',
                'fixymax False',
                ] + \
                ['%s %s' % (sequencesets[sequenceset]['epitopesbysitefile'], sequenceset.replace('_', ' ')) for sequenceset in sequencesets]
            ))
    os.system("%s %s" % (script, infile))


    # compare overall distributions
    script = 'epitopefinder_plotdistributioncomparison.py'
    print "Comparing distribution of all sites in different sequence sets..."
    sequencesetlist = sequencesets.keys()
    for isequenceset1 in range(len(sequencesets)):
        sequenceset1 = sequencesetlist[isequenceset1]
        countlist1 = [int(line.split(',')[1]) for line in open(sequencesets[sequenceset1]['epitopesbysitefile']).readlines()[1 : ] if not line.isspace()]
        for sequenceset2 in sequencesetlist[isequenceset1 + 1 : ]:
            countlist2 = [int(line.split(',')[1]) for line in open(sequencesets[sequenceset2]['epitopesbysitefile']).readlines()[1 : ] if not line.isspace()] 
            (d, p) = scipy.stats.ks_2samp(countlist1, countlist2)
            print "Kolmogorov-Smirnov test of whether we can reject null hypothesis that %s and %s are drawn from the same distribution: %.5f" % (sequenceset1, sequenceset2, p)
            infile = '%s_%s_vs_%s_infile.txt' % (os.path.splitext(script)[0], sequenceset1, sequenceset2)
            open(infile, 'w').write('\n'.join([\
                    '# input file for %s' % script,
                    'plotfile distributioncomparison_%s_vs_%s.pdf' % (sequenceset1, sequenceset2),
                    'epitopesfile1 %s' % sequencesets[sequenceset1]['epitopesbysitefile'],
                    'epitopesfile2 %s' % sequencesets[sequenceset2]['epitopesbysitefile'],
                    'set1name %s' % sequenceset1.replace('_', ' '),
                    'set2name %s' % sequenceset2.replace('_', ' '),
                    'pvalue None',
                    'title None',
                    ]))
            os.system('%s %s' % (script, infile))

    # compare all versus substituted
    script = 'epitopefinder_plotdistributioncomparison.py'
    print "Comparing distribution of all sites versus substituted sites for each sequence set..."
    for sequenceset in sequencesets:
        infile = '%s_%s_all_vs_substituted_infile.txt' % (os.path.splitext(script)[0], sequenceset)
        open(infile, 'w').write('\n'.join([\
                '# input file for %s' % script,
                'plotfile distributioncomparison_%s_all_vs_substituted.pdf' % sequenceset,
                'epitopesfile1 %s' % sequencesets[sequenceset]['epitopesbysitefile'],
                'epitopesfile2 %s' % sequencesets[sequenceset]['epitopecounts_substitutedsites'],
                'set1name all',
                'set2name substituted',
                'pvalue 100000',
                'pvaluewithreplacement False',
                'title %s influenza NP' % sequenceset.replace('_', ' '),
                'ymax 0.68',
                ]))
        os.system('%s %s' % (script, infile))

    # compare substituted for the sequence sets
    print "Comparing substituted sites among the sequence sets..."
    for isequenceset1 in range(len(sequencesets)):
        sequenceset1 = sequencesetlist[isequenceset1]
        alllist1 = [int(line.split(',')[1]) for line in open(sequencesets[sequenceset1]['epitopesbysitefile']).readlines()[1 : ] if not line.isspace()]
        countlist1 = [int(line.split(',')[1]) for line in open(sequencesets[sequenceset1]['epitopecounts_substitutedsites']).readlines()[1 : ] if not line.isspace()]
        for sequenceset2 in sequencesetlist[isequenceset1 + 1 : ]:
            alllist2 = [int(line.split(',')[1]) for line in open(sequencesets[sequenceset2]['epitopesbysitefile']).readlines()[1 : ] if not line.isspace()]
            countlist2 = [int(line.split(',')[1]) for line in open(sequencesets[sequenceset2]['epitopecounts_substitutedsites']).readlines()[1 : ] if not line.isspace()] 
            p = ProbDiffer(alllist1, countlist1, alllist2, countlist2)
            print "The probability that the mean number of epitopes per site for %s exceeds that for %s by more than expected by chance is %g" % (sequenceset1, sequenceset2, p)


    # compare epistatic versus all and substituted
    script = 'epitopefinder_plotdistributioncomparison.py'
    print "Comparing distribution of epistatic sites to all and substituted sites..."
    for sequenceset in sequencesets:
        if 'epistaticsites' not in sequencesets[sequenceset]:
            continue
        for (comparisonset, comparisonfile) in [('all', 'epitopesbysitefile'), ('substituted', 'epitopecounts_substitutedsites')]:
            infile = '%s_%s_%s_vs_epistatic_infile.txt' % (os.path.splitext(script)[0], sequenceset, comparisonset)
            open(infile, 'w').write('\n'.join([\
                    '# input file for %s' % script,
                    'plotfile distributioncomparison_%s_%s_vs_epistatic.pdf' % (sequenceset, comparisonset),
                    'epitopesfile1 %s' % sequencesets[sequenceset][comparisonfile],
                    'epitopesfile2 %s' % sequencesets[sequenceset]['epitopecounts_epistaticsites'],
                    'set1name %s' % comparisonset,
                    'set2name epistatic',
                    'pvalue 100000',
                    'pvaluewithreplacement True',
                    'title %s influenza NP' % (sequenceset.replace('_', ' ')),
                    ]))
            os.system('%s %s' % (script, infile))

    # make JPG and eps versions of all figures
    print "Creating JPG and EPS versions of all figures..."
    pdfs = [f for f in os.listdir('.') if os.path.splitext(f)[1] == '.pdf']
    for pdf in pdfs:
        os.system('convert -density 400 %s %s.jpg' % (pdf, os.path.splitext(pdf)[0]))
        os.system('pdf2ps %s' % pdf)
        os.system('ps2eps -l -C -f %s.ps' % os.path.splitext(pdf)[0])
        os.remove('%s.ps' % os.path.splitext(pdf)[0])


if __name__ == '__main__':
    main() # run the script
