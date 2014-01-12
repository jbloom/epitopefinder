"""Module for reading and handling epitopes.

This module is for handling information about epitopes. It is
designed principally with the goal of parsing epitopes from the
Immune Epitope Database (www.iedb.org).

Written by Jesse Bloom.


Functions defined in this module
------------------------------------
* `ReadCompactIEDB` : reads compact CSV downloads from Immune Epitope Database

* `SplitCSVLine` : returns entries in CSV file line (allowing quotations)

* `ReadSuperTypeClassification` : reads in supertype classification scheme.

* `AssignAlleleInfo` : assigns allele to mhcclass, mhcgene, supertype

* `PurgeRedundantEpitopes` : purges potentially redundant epitopes.

* `PrintAlleleSummary` : prints summary of allele assignments

* `EpitopeSummaryString` : string summary of epitope in CSV format.


Classes defined in this module
--------------------------------
* `Epitope` : a class for storing epitope information.


Details of classes and functions
---------------------------------
Details of classes and functions are provided in their individual
docstrings below.

"""


import os
import re


def PurgeRedundantEpitopes(epitopes, alignmentpositions, overlap, mhc_classification):
    """Purges potentially redundant epitopes.

    Given a list of epitopes (*epitopes*) and the positions to which they align
    in a target protein (*alignmentpositions*), find all epitopes that overlap
    by >= *overlap* sites in their alignment to the target. These epitopes are
    then examined to see if they clearly come from a different MHC classification
    (according to the option *mhc_classification*) -- if they do not, then
    the redundant epitopes are removed.

    For this process, the epitopes are first sorted so that those with
    the shortest epitope sequence (best defined sequence) come first.
    Among epitopes with the same sequence lengths, those with the
    most detailed MHC classification are sorted to come first. We then
    move down this list asking which epitopes are redundant with these
    epitopes. Moving from shortest to longest ensures that distinct
    short epitopes that overlap the same longer epitope are not considered
    redundant.

    CALLING VARIABLES:

    * *epitopes* is a list of *Epitope* objects specifying the epitopes.
      At a minimum, they should all have their *ep.sequence* attribute set
      to a string.

    * *alignmentpositions* is a list of 2-tuples of the same length
      as *epitopes*. For each epitope *epitopes[i]*, the 2-tuple
      *alignmentpositions[i] = (alignstart, alignend)* gives the position
      that this epitope aligns to in the target protein. The epitope aligns
      to positions *j* where *alignstart <= j <= alignend*. These
      alignment positions are used to determine the extent to which
      two epitopes overlap in their alignment.

    * *overlap* is a string specifying the amount of overlap that two epitopes
      must possess before they are considered potentially redundant. Any two
      epitopes that overlap (as determined by *alignmentpositions*) by >=
      *overlap* residues are potentially redundant.

    * *mhc_classification* is how we determine whether two epitopes with
      >= *overlap* alignment are indeed redundant. Possible values are
      the following strings:

       - *MHCgene* two overlapping epitopes are considered redundant
         if they share their *ep.mhcgene* attributes (or if one has 
         this attribute set to *None*).

       - *MHCsupertype* two overlapping epitopes are considered redundant
         if they share their *ep.mhcgene* attribute and if they share
         their *ep.supertype* attribute. An attribute set to *None*
         is considered to match any other value of the attribute.

       - *MHCgroup* two overlapping epitopes are considered redundant
         if they share their *ep.mhcgene*, *ep.supertype*, and *ep.mhcgroup*
         attributes. An attribute set to *None* is considered to match 
         any other value of the attribute.

    RETURN VARIABLE:

    This function returns the 4-tuple 
    *(unique_epitopes, unique_positions, redundant_epitopes, redundant_positions)*.
    All three elements of this 3-tuple are lists of the same length. The list
    *unique_epitopes* contains one representative of each set of redundant
    epitopes, and it always contains the one with the shortest *ep.sequence*.
    For each epitope *unique_epitopes[i]*, the element *unique_positions[i]*
    is a 2-tuple of numbers specifying the epitope alignment position as
    in *alignmentpositions*.
    For each element *i* of *unique_epitopes*, the element *redundant_epitopes[i]*
    is a list of all other epitopes redundant with *unique_epitopes[i]*. If there
    are no such redundant epitopes, then *redundant_epitopes[i]* is an empty list.
    For each element *j* of *redundant_epitopes[i]*, the element 
    *redundant_positions[i][j]* is the 2-tuple of numbers specifying the
    epitope alignment position for redundant epitope *redundant_epitopes[i][j]*
    as in *alignmentpositions*.
    """
    # sort epitopes / indices by length
    sorted_epitopes = []
    assert len(epitopes) == len(alignmentpositions)
    for i in range(len(epitopes)):
        ep = epitopes[i] 
        if not ep.sequence:
            raise ValueError("epitope does not have sequence")
        lacksgroup = lackssupertype = lacksgene = True
        if ep.mhcgroup:
            lacksgroup = False
        if ep.supertype:
            lackssupertype = False
        if ep.mhcgene:
            lacksgene = False
        sorted_epitopes.append((len(ep.sequence), lacksgroup, lackssupertype, lacksgene, ep, alignmentpositions[i]))
    sorted_epitopes.sort()
    # now look for redundant epitopes
    unique_epitopes = []
    unique_positions = []
    redundant_epitopes = []
    redundant_positions = []
    for i in range(len(sorted_epitopes)):
        (eplength, lacksgroup, lackssupertype, lacksgene, ep, (alignstart, alignend)) = sorted_epitopes[i]
        assert alignstart <= alignend
        for i2 in range(len(unique_epitopes)):
            (ep2, (alignstart2, alignend2)) = (unique_epitopes[i2], unique_positions[i2])
            epoverlap = max(min(alignend, alignend2) - max(alignstart, alignstart2) + 1, 0)
            if epoverlap >= overlap: # ep overlaps with ep2
                if mhc_classification == 'MHCgroup':
                    if (ep.mhcgene != ep2.mhcgene and ep.mhcgene and ep2.mhcgene) or (ep.supertype != ep2.supertype and ep.supertype and ep2.supertype) or (ep.mhcgroup != ep2.mhcgroup and ep.mhcgroup and ep2.mhcgroup):
                        continue # not redundant
                elif mhc_classification == 'MHCsupertype':
                    if (ep.supertype != ep2.supertype and ep.supertype and ep2.supertype) or (ep.mhcgene != ep2.mhcgene and ep.mhcgene and ep2.mhcgene):
                        continue # not redundant
                elif mhc_classification == 'MHCgene':
                    if ep.mhcgene != ep2.mhcgene and ep.mhcgene and ep2.mhcgene:
                        continue # not redundant
                else:
                    raise ValueError("Invalid mhc_classification of %s" % mhc_classification)
                redundant_epitopes[i2].append(ep)
                redundant_positions[i2].append((alignstart, alignend))
                break
        else: # ep is unique
            unique_epitopes.append(ep)
            unique_positions.append((alignstart, alignend))
            redundant_epitopes.append([])
            redundant_positions.append([])
    return (unique_epitopes, unique_positions, redundant_epitopes, redundant_positions)


def SplitCSVLine(line):
    """Splits CSV line where entries can be in quotes.

    This function takes a single argument *line*, which is a string
    specifying a line from a CSV (comma-separated values) file. The
    line is split into entries based on the comman, with the exception
    that entries themselves can be in quotes and a quoted entry is not
    split even if it has a comma. This makes this function different
    than just *line.split(',')* which will split entries even if
    the comma is in quotes. However, this function currently does not handle
    lines with quotes as part of the line entries.

    The returned value is a list of the entries in the line, with leading
    / trailing whitespace removed.

    >>> SplitCSVLine('entry1,entry2,entry3')
    ['entry1', 'entry2', 'entry3']

    >>> SplitCSVLine('"entry1","entry2,two",entry3 ')
    ['entry1', 'entry2,two', 'entry3']

    >>> SplitCSVLine('entry1,entry"2",entry3')
    Traceback (most recent call last):
       ...
    ValueError: entry contains internal quote

    >>> SplitCSVLine(' ')
    []

    >>> SplitCSVLine(',,entry3 ')
    ['', '', 'entry3']

    """
    assert isinstance(line, str)
    entries = []
    istart = i = 0
    quotestart = quoteend = False
    while i < (len(line.strip())):
        x = line[i]
        if x == ',' and not (quotestart and not quoteend): 
            if not quotestart:
                entries.append(line[istart : i].strip())
            elif quotestart and quoteend:
                quotestart = quoteend = False
                entries.append(line[istart + 1 : i - 1].strip())
            else:
                raise ValueError("Problem with line:\n%s" % line)
            i += 1
            istart = i
        elif quotestart and quoteend:
            raise ValueError("entry contains internal quote")
        elif x == '"':
            if i == istart:
                assert not quotestart
                quotestart = True
            elif quotestart:
                quoteend = True
            else:
                raise ValueError("entry contains internal quote")
            i += 1
        else:
            i += 1
    if quotestart and not quoteend:
        raise ValueError("entry contains internal quote")
    if quotestart and quoteend:
        entry = line[istart + 1 : i - 1].strip()
    else:
        entry = line[istart : i].strip()
    if entry or ',' in line:
        entries.append(entry)
    return entries




def ReadCompactIEDB(infile):
    """Reads epitopes from compact CSV downloads from Immune Epitope Database.

    This script reads the compact form of the CSV file downloads that
    can be made of sets of epitopes from the Immune Epitope Database 
    (www.iedb.org). If the format of these CSV files is radically 
    reconfigured then this function may no longer work -- howevever, it
    is designed to hopefully raise an Exception in that case rather than
    give spurious output. It functioned well on files downloaded
    on April-15-2013, and will presumably function on downloads from other
    dates although that has not been confirmed.

    *infile* should be a string giving the name of a CSV file download.
    It is assumed that this file contains only epitopes for which the assay
    result is 'Positive'. If it contains entries for which the result is
    not 'Positive', then raises an exception.

    The function returns a list *epitopelist* where each entry is
    an *Epitope* object containing the relevant information for an
    epitope in *infile*.
    """
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile:\n%s" % infile)
    indices = { # list of required entries, values will be assigned to indices
            'Host Organism Name':None,
            'Epitope Source Organism Name':None,
            'Epitope Source Molecule Name':None,
            'Epitope Linear Sequence':None,
            'MHC Allele Name':None,
            'Qualitative Measure':None,
            'Method/Technique':None,
            'Assay Group':None,
            'PubMed ID':None,
            'Author':None,
            'Journal':None,
            'Year':None,
            'Antigen Starting Position':None,
            'Antigen Ending Position':None,
            }
    # first read the header, make sure it has the required entries,
    # and assign the indices of these entries.
    f = open(infile)
    try:
        header = f.next().strip()
    except StopIteration:
        raise IOError("infile contains no lines:\n%s" % infile)
    header = [x.strip() for x in header.split(',')]
    nentries = len(header)
    for entry in indices.keys():
        n = header.count(entry)
        if n == 0:
            raise ValueError("No column in infile header for:\n%s" % entry)
        elif n == 1:
            indices[entry] = header.index(entry)
        else:
            raise ValueError("Duplicate columns in infile header for:\n%s" % entry)
    # now iterate over the remaining lines, which should contain the epitopes
    epitopes = []
    firstline = True
    for line in f:
        entries = SplitCSVLine(line)
        if not entries:
            continue # empty line
        if nentries != len(entries):
            if firstline and header[-1] == '' and nentries == len(entries) + 1:
                # this is a hack, but it appears necessary to fix the fact that
                # the IEDB files have a trailing comma on the header lines but
                # not the others.
                nentries -= 1
            else:
                raise ValueError("Line has incorrect number of entries, expected %d but found %d. Line is:\n%s" % (nentries, len(entries), line))
        if 'POSITIVE' not in entries[indices['Qualitative Measure']].upper():
            raise ValueError('Assay result is not positive, make sure you only downloaded positive assays:\n%s' % entries[indices['Qualitative Measure']])
        if entries[indices['Host Organism Name']]:
            host = entries[indices['Host Organism Name']]
        else:
            host = None
        if entries[indices['Epitope Source Organism Name']]:
            sourceorganism = entries[indices['Epitope Source Organism Name']]
        else:
            sourceorganism = None
        if entries[indices['Epitope Source Molecule Name']]:
            sourcemolecule = entries[indices['Epitope Source Molecule Name']]
        else:
            sourcemolecule = None
        if entries[indices['Epitope Linear Sequence']]:
            sequence = entries[indices['Epitope Linear Sequence']]
        else:
            sequence = None
        if entries[indices['MHC Allele Name']]:
            mhcallele = entries[indices['MHC Allele Name']]
        else:
            mhcallele = None
        assay = []
        for key in ['Method/Technique', 'Assay Group']:
            if entries[indices[key]]:
                assay.append(entries[indices[key]])
        if len(assay) > 1:
            assay = ' / '.join(assay)
        elif len(assay) == 1:
            assay = assay[0]
        else:
            assay = None
        reference = []
        for key in ['Author', 'Journal', 'Year', 'PubMed ID']:
            x = entries[indices[key]]
            if x:
                if key == 'PubMed ID':
                    reference.append('PMID %s' % x)
                else:
                    reference.append(x)
        if len(reference) > 1:
            reference = '. '.join(reference)
        elif len(reference) == 1:
            reference = reference[0]
        else:
            reference = None
        positionstart = entries[indices['Antigen Starting Position']]
        positionend = entries[indices['Antigen Ending Position']]
        if positionstart and positionend:
            position = (int(positionstart), int(positionend))
            if position[0] > position[1]:
                raise ValueError("Antigen position start is after position end: values are %d and %d" % (position))
        else:
            position = None
        ep = Epitope(host=host, sourceorganism=sourceorganism,\
                sourcemolecule=sourcemolecule, sequence=sequence,\
                mhcallele=mhcallele, assay=assay, reference=reference,\
                position=position)
        epitopes.append(ep)
        firstline = False
    f.close()
    return epitopes



def ReadSuperTypeClassification(supertypesfile, unclassified_to_none=True):
    """Reads in supertype classification scheme.

    *supertypesfile* is the name of a file specifying the *supertype*
    classifications. This file contains lines that list the alleles
    it should contain lines listing
    at the ``mhcgene*mhcgroup:protein`` level of detail, followed by the supertype.
    Supertypes of 'Unclassified' are assigned *None*::

            A*01:01 A01
            A*01:02 Unclassified
            A*68:15 A02
            B*07:43 B07
            B*08:01 B08

    *unclassified_to_none* is an optional argument which is *True* by
    default. It specifies that if an allele has a supertype classification
    of 'Unclassified', in the returned dictionary then the value is set to
    *None* for the supertype.

    This function returns a dictionary keyed by the alleles (at
    the ``mhcgene*mhcgroup:protein`` level, and with values being the
    supertype classification for that allele.
    """
    allelematch = re.compile('^(HLA\-){0,1}(?P<gene>A|B|C|DRB1|DRB3|DRB4|DRB5|DQA1|DQB1|DPA1|DPB1)\*(?P<group>\d{2})\:(?P<protein>\d{2,3})(\:\d{2})*[LSCAQ]*$')
    # read supertype_d, keyed by allele (gene*group:protein) values supertype
    supertype_d = {}
    if not os.path.isfile(supertypesfile):
        raise IOError("Cannot find file supertypesfile of %s" % supertypesfile)
    lines = [line.strip() for line in open(supertypesfile).readlines() if not line.isspace()]
    for line in lines:
        entries = line.split()
        if len(entries) != 2:
            raise ValueError("Invalid line in supertypesfile, does not contain two entries:\n%s" % line)
        m = allelematch.search(entries[0])
        if not m:
            raise ValueError("Failed to match allele: %s" % entries[0])
        key = '%s*%s:%s' % (m.group('gene'), m.group('group'), m.group('protein'))
        if key != entries[0]:
            raise ValueError("Parsed allele doesn't match line")
        if key in supertype_d:
            raise ValueError("Duplicate supertype assignment for %s" % key)
        if unclassified_to_none and entries[1].upper() == 'UNCLASSIFIED':
            supertype_d[key] = None
        else:
            supertype_d[key] = entries[1]
    return supertype_d



def PrintAlleleSummary(eps, out):
    """Summarizes allele classifications for epitopes.

    *eps* should be a list of *Epitope* objects. Each of them should
    have been passed through *AssignAlleleInfo*.

    *out* is a writable file-like object.

    Writes a summary to *out* of the allele assignments for the
    epitopes in *eps*, as based on the *mhcclass*, *mhcgene*, *mhcgroup*,
    *supertype* attributes of the objects in this list. This summary
    contains information such as how many have each of these attributes
    assigned.
    """
    n = len(eps)
    out.write("\nSummary of allele assignments for the %d epitopes.\n" % n)
    classi = [ep for ep in eps if ep.mhcclass == 'I']
    classii = [ep for ep in eps if ep.mhcclass == 'II']
    for (mhcclass, eplist) in [('I', classi), ('II', classii)]:
        out.write('%d of %d epitopes were assigned to MHC class %s\n' % (len(eplist), n, mhcclass))
        withgene = [ep for ep in eplist if ep.mhcgene]
        withgroup = [ep for ep in eplist if ep.mhcgroup]
        genes = dict([(ep.mhcgene, True) for ep in withgene]).keys()
        genes.sort()
        out.write('  Of these, %d of %d were assigned an MHC gene (%s)\n' % (len(withgene), len(eplist), ', '.join(genes)))
        if mhcclass == 'I':
            withsupertype = [ep for ep in eplist if ep.supertype]
            supertypes = dict([(ep.supertype, True) for ep in withsupertype]).keys()
            supertypes.sort()
            out.write('  Of these, %d of %d were assigned an MHC supertype (%s)\n' % (len(withsupertype), len(eplist), ', '.join(supertypes)))
        out.write('  Of these, %d of %d were assigned an MHC group within the MHC gene\n' % (len(withgroup), len(eplist)))
    nunassigned = n - len(classi) - len(classii)
    if nunassigned:
        out.write("%d of %d epitopes were not assigned a MHC class.\n" % (nunassigned, n))



def AssignAlleleInfo(ep, supertype_d):
    """Assigns an MHC allele to its mhcclass, mhcgene, and supertype.

    Currently this only assigns alleles for *ep.host* equal
    to 'Homo sapiens' or some subset (such as 'Homo sapiens Caucasian').

    CALLING VARIABLES: 

    * *ep* is an epitope object that should have its *ep.mhcallele*
      and *ep.host* attributes assigned. Using these attributes, attempts
      to assign *ep.mhcclass*, *ep.mhcgene*, *mhc.group*, and *ep.supertype*
      attributes. On completion of this function, *ep* will have these
      attributes assigned if possible. Note, however, that some of them still
      may be *None* if for example *ep.mhcallele* specifically
      does not contain this information (for example if it is
      'HLA-Class II, allele undetermined', or if *ep.mhcallele* is *None*.

    * *supertype_d* is a dictionary classifying alleles to supertypes,
      in the format returned by *ReadSuperTypeClassifications*. It should
      assign supertypes for all class I alleles.

    RETURN VALUE

    On return, *ep* will have its *mhcclass*, *mhcgene*, *mhcgroup*,
    and *supertype* attributes updated. Note that they may still be
    *None* if they are not assignable from *ep.mhcallele*.

    This function will raise an exception if *ep.mhcallele* cannot be
    processed.
    """
    allelematch = re.compile('^(HLA\-){0,1}(?P<gene>A|B|C|DRB1|DRB3|DRB4|DRB5|DQA1|DQB1|DPA1|DPB1)\*(?P<group>\d{2})\:(?P<protein>\d{2,3})(\:\d{2})*[LSCAQ]*$')
    supertypematch = re.compile('^HLA\-(?P<gene>[AB])(?P<group>\d+)$')
    classiimatch = re.compile('^HLA\-(?P<gene>D\w+)$')
    mixedclassiimatch = re.compile('^HLA\-D[RQ]A.*/D[RQ]B\S*$')
    # now parse through epitopes
    if 'Homo sapiens' in ep.host:
        if not ep.mhcallele:
            return # no allele to parse
        m = allelematch.search(ep.mhcallele)
        if not m:
            if ep.mhcallele == 'HLA-Class I, allele undetermined':
                ep.mhcclass = 'I'
                return # this is all we can parse
            elif ep.mhcallele == 'HLA-Class II, allele undetermined':
                ep.mhcclass = 'II'
                return # this is all we can parse
            elif 'mutant' in ep.mhcallele:
                return # don't assign mutant alleles
            elif mixedclassiimatch.search(ep.mhcallele):
                ep.mhcclass = 'II'
                return
            else:
                m = supertypematch.search(ep.mhcallele)
                if m:
                    if len(m.group('group')) == 1:
                        supertype = '%s0%s' % (m.group('gene'), m.group('group'))
                    else:
                        supertype = '%s%s' % (m.group('gene'), m.group('group'))
                    if supertype in supertype_d.values():
                        ep.mhcclass = 'I'
                        ep.mhcgene = m.group('gene')
                        ep.supertype = supertype
                        return
                    elif supertype in ['A11', 'A26', 'A28', 'A68', 'A69', 'B14', 'B35', 'B37', 'B39']:
                        # not assigned a supertype, so just give mhcclass and mhcgene
                        ep.mhcclass = 'I'
                        ep.mhcgene = m.group('gene')
                    else:
                        raise ValueError("Unrecognized supertype of %s in %s" % (supertype, ep.mhcallele))
                else:
                    m = classiimatch.search(ep.mhcallele)
                    if m:
                        ep.mhcclass = 'II'
                        if m.group('gene') in ['DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1']:
                            ep.mhcgene = m.group('gene')
                        elif m.group('gene') == 'DR':
                            pass # not enough information to assign mhcgene
                        # the remaining class II mhcgene assignments are based on
                        # Marsh et al, Tissue Antigens, 75:291-455 (2010)
                        elif m.group('gene') in ['DR4', 'DR7', 'DR15', 'DR1', 'DR2', 'DR11', 'DR3']:
                            ep.mhcgene = 'DRB1'
                        elif m.group('gene') in ['DR53']:
                            ep.mhcgene = 'DRB4'
                        elif m.group('gene') in ['DQ5']:
                            ep.mhcgene = 'DQB1'
                        elif m.group('gene') in ['DQ1', 'DP', 'DQ']:
                            pass # cannot assign mhcgene
                        else:
                            raise ValueError("Unassigned class II: %s" % ep.mhcallele)
                        return
                    elif 'allele undetermined' in ep.mhcallele:
                        pass
                    else:
                        raise ValueError("Cannot parse mhcallele: %s" % ep.mhcallele)
        else:
            ep.mhcgene = m.group('gene')
            ep.mhcgroup = m.group('group')
            if ep.mhcgene in ['A', 'B', 'C']:
                ep.mhcclass = 'I'
                key = "%s*%s:%s" % (ep.mhcgene, ep.mhcgroup, m.group('protein'))
                try:
                    ep.supertype = supertype_d[key]
                    return
                except KeyError:
                    raise ValueError("No supertype specified for %s" % key)
            elif ep.mhcgene in ['DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1']:
                ep.mhcclass = 'II'
                return
            else:
                raise ValueError("Parsed mhcgene but can't assign class: %s" % ep.mhcgene)
    else: # cannot assign for this source species
        raise ValueError("Currently only handles host of 'Homo sapiens', this epitope has: %s" % ep.host)


def EpitopeSummaryString(ep, position):
    """Returns a string summary of epitope in CSV format.

    *ep* is an *Epitope* object.

    *position* is a 2-tuple of the format *(start, end)* indicating the
    starting and ending position to which *ep* aligns in the target sequence.

    Returns a string in CSV format (entries each surrounded by quotes and
    separated by commas) of epitope *ep* and its alignment position. There are
    five entries:

        1) The epitope sequence.

        2) The starting and ending alignment position separated by a dash.

        3) The MHC allele information.

        4) Information about the epitope source / host and alignment
           position as taken from the source.

        5) The reference for the epitope.

    Here is an example:

    >>> ep = Epitope(assay='51 chromium release, cytotoxicity', mhcclass='I', reference='L G Tussey; S Rowland-Jones; T S Zheng; M J Androlewicz; P Cresswell; J A Frelinger; A J McMichael. Immunity. 1995. PMID 7542549', sequence='LRSRYWAI', sourcemolecule='NP', mhcallele='HLA-B*27:02', host='Homo sapiens', sourceorganism='Influenza A virus (A/X-31(H3N2))', mhcgroup='27', supertype='B27', position=(381, 388), mhcgene='B')
    >>> position = (381, 388)
    >>> s = EpitopeSummaryString(ep, position)
    >>> print s
    "LRSRYWAI","381-388","allele HLA-B*27:02; MHC class I; MHC gene B; MHC supertype B27; MHC group 27","host = Homo sapiens; source = Influenza A virus (A/X-31(H3N2)); molecule = NP; reported position = (381, 388); assay = 51 chromium release, cytotoxicity","L G Tussey; S Rowland-Jones; T S Zheng; M J Androlewicz; P Cresswell; J A Frelinger; A J McMichael. Immunity. 1995. PMID 7542549"

    """
    s = []
    s.append('"%s"' % ep.sequence)
    s.append('"%d-%d"' % position)
    if ep.mhcclass == 'I':
        s.append('"allele %s; MHC class %s; MHC gene %s; MHC supertype %s; MHC group %s"' % (ep.mhcallele, ep.mhcclass, ep.mhcgene, ep.supertype, ep.mhcgroup))
    else:
        s.append('"allele %s; MHC class %s; MHC gene %s; MHC group %s"' % (ep.mhcallele, ep.mhcclass, ep.mhcgene, ep.mhcgroup))
    s.append('"host = %s; source = %s; molecule = %s; reported position = %s; assay = %s"' % (ep.host, ep.sourceorganism, ep.sourcemolecule, ep.position, ep.assay))
    s.append('"%s"' % ep.reference)
    s = ','.join(s)
    if s.count('"') != 10:
        raise ValueError("String contains extra quotation marks:\n%s" % s)
    return s




class Epitope(object):
    """Class for storing epitope information.

    This class defines *Epitope* objects, which can be used to store
    information about epitopes.

    Each Epitope object *ep* possesses the following attributes. If the attribute
    is not defined for an epitope on initialization, then that
    attribute is set to *None*:

        * *ep.host* : a string giving the host organism, for example 'Homo sapiens'

        * *ep.sourceorganism* : a string giving the source organism, for example
          'Influenza A virus'

        * *ep.sourcemolecule* : a string giving the source molecule, for example
          'Nucleoprotein'

        * *ep.sequence* : a string giving the sequence of the epitope, 
          for example 'GILGFVFTL'. If *sequence* is assigned upon
          initialization of *ep*, then it will be converted to all upper case.

        * *ep.mhcallele* : a string giving the MHC allele, for example 
          'HLA-B*35:01'

        * *ep.assay* : a string giving the assay used identify the epitope, 
          for example 'ELISPOT; cytokine release IFNg'

        * *ep.reference* : a string giving the reference for the epitope, for
          example: 'T Linnemann; G Jung; P Walden. J Virol. (2000). PMID 10954576'

        * *ep.position* : a 2-tuple giving the position in the protein as
          the starting and ending integer sequence positions, for example
          *(46, 54)*

        * *ep.mhcclass* : a string 'I' or 'II' specifying whether the epitope
          is MHC class I or MHC class II.

        * *ep.mhcgene* : a string specifying the MHC gene. For example,
          for humans could be 'A', 'B', 'C', 'DQA1', 'DQB1',
          'DPA1', 'DPB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5'.

        * *ep.supertype* : a string giving the supertype assigned to
           an allele. For example, 'A01' or 'B27'.

        * *ep.mhcgroup* : a string specifying the MHC group. For example,
          for 'HLA-B*35:01' this would be '35'.

    To initialize an *Epitope* object, there are no required arguments, however,
    each of the above attributes can be assigned at initialization (unassigned
    attributes are set to *None*). For instance::

        ep = Epitope(host='Homo sapiens', sequence='GILGFVFTL')

    returns an *Epitope* object *ep* with *ep.host* set to 'Homo sapiens', 
    *ep.sequence* set to 'GILGFVFTL', and all other attributes set to
    *None*.
    """

    def __init__(self, host=None, sourceorganism=None, sourcemolecule=None,\
            sequence=None, mhcallele=None, assay=None, reference=None,\
            position=None, mhcclass=None, mhcgene=None, mhcgroup=None,\
            supertype=None):
        """Initializes a new *Epitope* object with the specified attributes."""
        assert host == None or isinstance(host, str)
        self.host = host
        assert sourceorganism == None or isinstance(sourceorganism, str)
        self.sourceorganism = sourceorganism
        assert sourcemolecule == None or isinstance(sourcemolecule, str)
        self.sourcemolecule = sourcemolecule
        assert sequence == None or isinstance(sequence, str)
        self.sequence = sequence.upper()
        assert mhcallele == None or isinstance(mhcallele, str)
        self.mhcallele = mhcallele
        assert assay == None or isinstance(assay, str)
        self.assay = assay
        assert reference == None or isinstance(reference, str)
        self.reference = reference
        assert position == None or (isinstance(position, tuple) and len(position) == 2)
        self.position = position
        assert mhcclass == None or isinstance(mhcclass, str)
        self.mhcclass = mhcclass
        assert mhcgene == None or isinstance(mhcgene, str)
        self.mhcgene = mhcgene
        assert mhcgroup == None or isinstance(mhcgroup, str)
        self.mhcgroup = mhcgroup
        assert supertype == None or isinstance(supertype, str)
        self.supertype = supertype



# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()

