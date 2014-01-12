"""This file contains utilities for manipulating sequences.

Performs functions such as reading / writing FASTA files, 
translating sequences, etc.

Written by Jesse Bloom."""



import os
import re
import datetime
import cStringIO


def ClassifyBySeason(genomes, startyear, endyear, subsample, season):
    """
    """
    pass


def ParseInfluenzaGenomes(infile, genes):
    """Parses influenza genomes from FASTA file.

    This function is designed to read a FASTA file that contains
    each of the proteins or genes for an influenza genome, and
    then return all strains with full genomes.

    CALLING VARIABLES:

    * *infile* is a string giving the name of the FASTA file. The
      headers should be of the following format::

        >AAX11456 A/New York/61A/2003(H3N2) 2003/12/20 M1

        >AAB06984 A/Louisiana/4/93(H3N2) 1993// NA

        >CAD29965 A/Panama/2007/1999 1999// NA

      These headers given the sequence number, the strain name (it
      is optional whether the subtype, see the third example header),
      the year/month/date of isolation, and the gene.

    * *genes* is a list of all required gene segments, as listed
      in the header. For example, you might commonly have
      *genes = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M1', 'M2', 'NS1', 'NS2']*.

    RETURN VARIABLE

    This function returns the dictionary *genomes*. This dictionary is keyed
    by the 4-tuples *(strain, year, month, day)*. If no month is specified
    then *month* is *None*, otherwise it is the numeric month. Likewise for 
    *day*. So the keys for the examples above would be:

     * *('A/New York/61A/2003(H3N2)', 2003, 12, 20)*

     * *('A/Louisiana/4/93(H3N2)', 1993, None, None)*

     * *('A/Panama/2007/1999', 1999, None, None)*

    The value for each key is another dictionary. It is keyed by each
    string *gene* in *genes*, and the values are the sequences for those
    genes. The sequences are all in upper case.

    *genomes* only contains those strains for which each *gene* in *genes*
    has a sequence specified. If there are multiple sequences provided for
    a gene, takes the first one encountered. If there are multiple full genome
    entries for the same strain name with some providing more detailed month/day
    information, takes the one with more detailed month/day information.

    """
    headmatch = re.compile('^\S+ (?P<strain>[ABC]\/[^\/]+\/[\w\-\.]+\/[\w\(\)]+) (?P<year>\d{4})\/(?P<month>\d*)\/(?P<day>\d*) (?P<gene>[\w\-]+)$')
    seqs = Read(infile)
    genomes = {}
    for (head, seq) in seqs:
        m = headmatch.search(head)
        if not m:
            raise ValueError("Failed to match header:\n%s" % head)
        gene = m.group('gene')
        if gene not in genes:
            continue # we are not saving this gene
        seq = seq.upper()
        strain = m.group('strain')
        year = int(m.group('year'))
        if m.group('month'):
            month = int(m.group('month'))
        else:
            month = None
        if m.group('day'):
            day = int(m.group('day'))
        else:
            day = None
        key = (strain, year, month, day)
        if key not in genomes:
            genomes[key] = {gene:seq}
        else:
            genomes[key][gene] = seq
    # remove all non full-length
    for key in genomes.keys():
        if len(genomes[key].keys()) != len(genes):
            del genomes[key]
    # if duplicate keys for a strain, keep one with most detailed dates
    strain_d = {}
    for (strain, year, month, day) in genomes.keys():
        if strain in strain_d:
            strain_d[strain].append((strain, year, month, day))
        else:
            strain_d[strain] = [(strain, year, month, day)]
    for (strain, keylist) in strain_d.iteritems():
        if len(keylist) > 1:
            best = keylist[0]
            others = []
            for key in keylist[1 : ]:
                if best[2] == None:
                    if key[2] != None:
                        others.append(best)
                        best = key
                    else:
                        others.append(key)
                elif best[3] == None:
                    if key[3] != None:
                        others.append(best)
                        best = key
                    else:
                        others.append(key)
            assert best not in others
            assert len(others) == len(keylist) - 1
            for key in others:
                del genomes[key]
    return genomes





def WriteNEXUS(seqs, outfile, seqtype):
    """Write sequences to a NEXUS file.

    *seqs* is a list of 2-tuples specifying sequences as (name, sequence). 
    The sequences should all be aligned so that they are of the same
    length.

    *outfile* is the name of the output file we are writing, such as 'file.nex'

    *seqtype* is the sequence type, such as 'PROTEIN' or 'DNA'
    """
    ntax = len(seqs)
    if not ntax:
        raise ValueError("no sequences specified")
    nchar = len(seqs[0][1])
    f = open(outfile, 'w')
    f.write("#NEXUS\n\nBegin DATA;\n\tDimensions ntax=%d nchar=%d;\n\tFormat datatype=%s gap=-;\n\tMatrix\n" % (ntax, nchar, seqtype))
    for (name, seq) in seqs:
        if len(seq) != nchar:
            raise ValueError("all sequences not of the same length")
        f.write("\t%s %s\n" % (name, seq))
    f.write(';\nEnd;')



def Write(headers_seqs, filename, writable_file=False):
    """Writes sequences to a FASTA file.

    *headers_seqs* is a list of 2-tuples specifying sequences and their
    corresponding headers.  Each entry is the 2-tuple *(header, seq)*
    where *header* is a string giving the header (without the leading ">"),
    and *seq* is the corresponding sequence.

    *filename* is a string that specifies the name of the file to which the
    headers and sequences should be written.  If this file already exists,
    it is overwritten. 

    *writable_file* is a Boolean switch specifying that rather than *filename*
    giving a string specifying the name of a file to which the sequences
    should be written, it instead specifies a writable file object to which
    the sequences should be written.

    The sequences are written to the file in the same order that they are specified
    in *headers_seqs*.
    """
    assert isinstance(writable_file, bool)
    if writable_file:
        f = filename
    else:
        f = open(filename, 'w')
    for (header, seq) in headers_seqs:
        f.write(">%s\n%s\n" % (header, seq))
    if not writable_file:
        f.close()


def Read(fastafile):
    """Reads sequences from a FASTA file.

    *fastafile* should specify the name of a FASTA file.
    
    This function reads all sequences from the FASTA file.  It returns the
    list *headers_seqs*.  This list is composed of a 2-tuple '(header, seq)'
    for every sequence entry in FASTA file.  'header' is the header for
    a sequence, with the leading ">" and any trailing spaces removed. 'seq'
    is the corresponding sequence.
    """
    lines = open(fastafile).readlines()
    headers_seqs = []
    header = None
    seq = []
    for line in lines:
        if line[0] == '>':
            if (not header) and (not seq):
                pass # first sequence in file
            elif header and not seq:
                raise ValueError, "Empty sequence for %s" % header
            elif seq and not header:
                raise ValueError, "File does not begin with header."
            else:
                seq = ''.join(seq)
                seq = seq.replace(' ', '')
                headers_seqs.append((header, seq))
            header = line.strip()[1 : ]
            seq = []
        else:
            seq.append(line.strip())
    if (not header) and (not seq):
        pass # first sequence in file
    elif header and not seq:
        raise ValueError, "Empty sequence for %s" % header
    elif seq and not header:
        raise ValueError, "File does not begin with header."
    else:
        seq = ''.join(seq)
        seq = seq.replace(' ', '')
        headers_seqs.append((header, seq))
    return headers_seqs



def GetSequence(header, headers_sequences):
    """Gets a particular sequence based on its header name.

    *header* specifies the name of a sequence's FASTA header.

    *headers_sequences* is a list of tuples '(head, seq)' as would
    be returned by **Read**.

    This function searches through *headers_sequences* and returns the
    sequence corresponding to the first header found that matches
    the calling argument *header*.  If no such header is found,
    raises an exception.
    """
    for (head, seq) in headers_sequences:
        if head == header:
            return seq
    else:
        raise ValueError, "Could not find a header matching %s" % header



def Translate(headers_sequences, readthrough_n=False, readthrough_stop=False, truncate_incomplete=False, translate_gaps=False):
    """Translates a set of nucleotide sequences to amino acid sequences.

    This function takes as input a single calling argument *header_sequences*,
    which is a list of tuples *(header, seq)* as would be returned by
    `Read`.  The sequences should all specify valid coding nucleotide
    sequences.  The returned variable is a new list in which
    all of the nucleotide sequences have been translated to their
    corresponding protein sequences, given by one letter codes. Stop codons
    are translated to '*'.
    If any of the nucleotide sequences do not translate to
    valid protein sequences, an exception is raised.

    The optional argument *readthrough_n* specifies that if any nucleotides
    in the sequence are equal to to an ambiguous nt code and cannot therefore
    be unambiguously translated into an amino acid, we simply translate through these 
    nucleotides by making the corresponding amino acid equal to "X".  By
    default, this option is *False*.  Note that even when this option is *False*,
    certain ambiguous nucleotides may still be translatable if they all lead to 
    the same amino acid.

    The optional argument *readthrough_stop* specifies that if we encounter
    any stop codons, we simply translation them to '*'.  By default,
    this option is *False*, meaning that we instead raise an error
    of an incomplete stop codon.

    The optional argument *truncate_incomplete* specifies that if the sequence
    length is not a multiple of three, we simply truncate off the one or two
    final nucleotides to make the length a multiple of three prior to translation.
    By default, this option is *False*, meaning that no such truncation is done.

    The optional argument *translate_gaps* specifies that a codon containing '-' 
    is translated to '-'.

    >>> Translate([('seq1', 'ATGTAA'), ('seq2', 'gggtgc')])
    [('seq1', 'M*'), ('seq2', 'GC')]

    >>> Translate([('seq2', 'GGNTGC')])
    [('seq2', 'GC')]
    
    >>> Translate([('seq2', 'NGGTGC')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate codon NGG
    
    >>> Translate([('seq2', 'NGGTGC')], readthrough_n=True)
    [('seq2', 'XC')]

    >>> Translate([('seq2', 'TAATGC')])
    Traceback (most recent call last):
       ...
    ValueError: Premature stop codon

    >>> Translate([('seq2', 'TAATGC')], readthrough_stop=True)
    [('seq2', '*C')]

    >>> Translate([('seq2', 'TGCA')])
    Traceback (most recent call last):
       ...
    ValueError: Sequence length is not a multiple of three

    >>> Translate([('seq2', 'TGCA')], truncate_incomplete=True)
    [('seq2', 'C')]

    >>> Translate([('seq2', 'TGC---')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate gap.

    >>> Translate([('seq2', 'TGC---')], translate_gaps=True)
    [('seq2', 'C-')]
    """
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
            'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'STOP', 'TAG':'STOP',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
            'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'TGT':'C', 'TGC':'C', 'TGA':'STOP', 'TGG':'W', 'CGT':'R',
            'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
            'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    assert isinstance(headers_sequences, list)
    translated_headers_sequences = []
    for (head, seq) in headers_sequences:
        seq = seq.upper()
        if len(seq) % 3:
            if truncate_incomplete:
                seq = seq[ : -(len(seq) % 3)]
            else:
                raise ValueError, "Sequence length is not a multiple of three"
        prot_length = len(seq) // 3
        prot = []
        for i in range(prot_length):
            codon = seq[3 * i : 3 * (i + 1)]
            try:
                aa = genetic_code[codon]
            except KeyError:
                if '-' in codon:
                    if translate_gaps:
                        aa = '-'
                    else:
                        raise ValueError("Cannot translate gap.")
                else:
                    # see if we have an ambiguous nucleotide codon that doesn't matter in the translation
                    possible_nt1 = AmbiguousNTCodes(codon[0])
                    possible_nt2 = AmbiguousNTCodes(codon[1])
                    possible_nt3 = AmbiguousNTCodes(codon[2])
                    possible_codons = []
                    for nt1 in possible_nt1:
                        for nt2 in possible_nt2:
                            for nt3 in possible_nt3:
                                possible_codons.append("%s%s%s" % (nt1, nt2, nt3))
                    try:
                        aa = genetic_code[possible_codons[0]]
                    except KeyError:
                        raise KeyError("Cannot translate codon %s in %s" % (codon, head))
                    for possible_codon in possible_codons:
                        if genetic_code[possible_codon] != aa:
                            if readthrough_n:
                                aa = 'X'
                            else:
                                raise ValueError("Cannot translate codon %s" % codon)
            if aa == 'STOP' and i == prot_length - 1:
                aa = '*'
            elif aa == 'STOP':
                if readthrough_stop:
                    aa = '*'
                else:
                    raise ValueError("Premature stop codon")
            prot.append(aa)
        translated_headers_sequences.append((head, ''.join(prot)))
    return translated_headers_sequences



def UnknownsToGaps(headers_sequences):
    """Converts all unknown ambiguous amino acids to gap characters.

    This function converts all of the common codes for unknown amino acids
    to gaps.  That is, "B" (Asp or Asn), "Z" (Glu or Gln), and "X" (any
    amino acid) are all converted to the character "-", which is usually
    taken to denote a gap.  You would use this method if you are 
    subsequently processing the sequences with a program that recognizes
    the gap character but not these three other unknown characters.  Note,
    of course, that you are changing your sequence, so don't do this unless
    you need to.

    This method takes as a single calling variable the list *header_sequences*
    as would be returned by `Read`.  It returns a copy of this list, but
    with all unknown amino acids replaced by "-".
    """
    new_headers_sequences = []
    for (head, seq) in headers_sequences:
        seq = seq.replace('B', '-')
        seq = seq.replace('Z', '-')
        seq = seq.replace('X', '-')
        new_headers_sequences.append((head, seq))
    return new_headers_sequences



def PurgeDuplicates(headers_sequences):
    """Removes all duplicate sequences and those that are substrings.

    This function takes a single calling argument *headers_sequences*, which
    is a list of tuples *(header, seq)* as would be returned by `Read`.
    It returns a new list in which each sequence appears exactly once.
    If a sequence appears more than once in the original calling
    list, the sequence that is kept is the first one encountered in
    the list. Sequences that are substrings of others
    are also removed.

    The ordering of sequences is preserved except for the removal
    of duplicates.  Sequences are returned as all upper case.

    >>> PurgeDuplicates([('seq1', 'atgc'), ('seq2', 'GGCA'), ('seq3', 'ATGC')])
    [('seq1', 'ATGC'), ('seq2', 'GGCA')]

    >>> PurgeDuplicates([('seq1', 'atgcca'), ('seq2', 'TGCC')])
    [('seq1', 'ATGCCA')]
    """
    unique_headers_sequences = []
    assert isinstance(headers_sequences, list)
    for (head, seq) in headers_sequences:
        assert isinstance(seq, str)
        seq = seq.upper()
        i = 0
        for (head2, seq2) in unique_headers_sequences:
            if seq == seq2 or seq in seq2:
                break
            elif seq2 in seq:
                unique_headers_sequences[i] = (head, seq)
                break
            i += 1
        else:
            unique_headers_sequences.append((head, seq))
    return unique_headers_sequences



def AmbiguousNTCodes(nt):
    """Returns all possible nucleotides corresponding to an ambiguous code.

    This method takes as input a single nucleotide character *nt*, which is
    assumed to represent a nucleotide as one of the accepted
    codes for an ambiguous character.  Returns a list giving
    all possible codes for which a nucleotide might stand.  Raises
    an exception if *nt* is not a valid nucleotide code.

    >>> AmbiguousNTCodes('N')
    ['A', 'T', 'G', 'C']

    >>> AmbiguousNTCodes('R')
    ['A', 'G']

    >>> AmbiguousNTCodes('A')
    ['A']

    >>> AmbiguousNTCodes('-')
    ['-']

    >>> AmbiguousNTCodes('F')
    Traceback (most recent call last):
       ...
    ValueError: Invalid nt code of "F"
    """
    if nt in ['A', 'T', 'G', 'C', '-']:
        return [nt]
    elif nt == 'R':
        return ['A', 'G']
    elif nt == 'Y':
        return ['T', 'C']
    elif nt == 'K':
        return ['G', 'T']
    elif nt == 'M':
        return ['A', 'C']
    elif nt == 'S':
        return ['G', 'C']
    elif nt == 'W':
        return ['A', 'T']
    elif nt == 'B':
        return ['C', 'G', 'T']
    elif nt == 'D':
        return ['A', 'G', 'T']
    elif nt == 'H':
        return ['A', 'C', 'T']
    elif nt == 'V':
        return ['A', 'C', 'G']
    elif nt == 'N':
        return ['A', 'T', 'G', 'C']
    else: 
        raise ValueError('Invalid nt code of "%s"' % nt)



def PurgeAmbiguousDNA(headers_sequences):
    """Removes all sequences with ambiguous positions from nucleotide sequences.

    This function takes a single calling argument *headers_sequences*, which
    is a list of tuples *(header, seq)* as would be returned by `Read`.
    These sequences should specify nucleotide sequences.  It returns
    a new list which is a copy of *headers_sequences*, except that all
    sequences that contain ambiguous nucleotide entries (i.e. characters
    that are 'A', 'T', 'C', 'G', 'a', 't', 'c', or 'g') have been
    removed.
    """
    valid_bases = {'A':1, 'T':1, 'C':1, 'G':1, 'a':1, 't':1, 'c':1, 'g':1}
    assert isinstance(headers_sequences, list)
    purged_headers_sequences = []
    for (head, seq) in headers_sequences:
        for nt in seq:
            if nt not in valid_bases:
                break
        else:
            purged_headers_sequences.append((head, seq))
    return purged_headers_sequences




def GetEntries(namelist, fastafile, allow_substring=False):
    """Gets selected entries from a (potentially very large) FASTA file.

    This method is designed to extract sequences from a FASTA file.  It will work
    even if the FASTA file is very large, since it avoids reading the entire
    file into memory at once. 

    *namelist* specifies the "names" of the sequences that we want to extract from
    the FASTA file.  The "name" of a sequence is the string immediately following
    the ">" in the FASTA file header for a sequence, terminated by a space character
    (space, tab, or return).  So for example, the header::

            >E_coli_thioredoxin: the thioredoxin protein from E. coli

    would correspond to a name of "E_coli_thioredoxin".  'namelist' specifies
    a list of these names.

    *fastafile* is the name of a FASTA file that contains the sequences we are searching
    for.  For this method to be guaranteed to work properly, each sequence in the FASTA
    file must contain a unique name, where a "name" is as defined above.  Note that this
    uniqueness of names is not rigorously checked for, so if there are not unique names,
    the function may raise an exception, or it may continue along and give no hint of
    the problem.

    *allow_substring* is an optional Boolean switch that specifies that the name given in
    *namelist* need only be a substring of the first entry in the *fastafile* header.

    The function expects to find exactly one entry in *fastafile* for each name listed in
    *namelist*.  If it does not, it will raise an exception.  The returned variable
    is a list composed of 2-tuples.  Element i of this list corresponds to the name
    given by *namelist[i]*.  Each 2-tuple has the form '(header, sequence)' where
    'header' is the full FASTA header, but with the leading ">" character and any trailing
    linebreaks/spaces removed.  'sequence' is a string giving the sequence, again with the
    trailing linebreak removed.
    """
    namedict = {}
    for name in namelist:
        namedict[name] = None
    f = open(fastafile)
    line = f.readline()
    header = None
    seq = []
    while line:
        if line[0] == '>':
            if (not header) and (not seq):
                pass # first sequence in file
            elif header and not seq:
                raise ValueError, "Empty sequence for %s" % header
            elif seq and not header:
                raise ValueError, "File does not begin with header."
            else:
                name = header.split()[0]
                if name in namedict:
                    if namedict[name]:
                        raise ValueError, "Duplicate entries for name %s" % name
                    else:
                        namedict[name] = (header, ''.join(seq))
                elif allow_substring:
                    subnames = [iname for iname in namedict.iterkeys() if iname in name]
                    if len(subnames) > 1:
                        raise ValueError("Multiple subname matches for %s." % name)
                    elif subnames:
                        if namedict[subnames[0]]:
                            raise ValueError("Duplicate subname entries for name %s" % name)
                        else:
                            namedict[subnames[0]] = (header, ''.join(seq))
            header = line.strip()[1 : ]
            seq = []
        else:
            seq.append(line.strip())
        line = f.readline()
    f.close()
    if (not header) and (not seq):
        pass # first sequence in file
    elif header and not seq:
        raise ValueError, "Empty sequence for %s" % header
    elif seq and not header:
        raise ValueError, "File does not begin with header."
    else:
        name = header.split()[0]
        if name in namedict:
            if namedict[name]:
                raise ValueError, "Duplicate entries for name %s" % name
            else:
                namedict[name] = (header, ''.join(seq))
    header_seq_list = []
    for name in namelist:
        if not namedict[name]:
            raise ValueError, "No entry for %s" % name
        header_seq_list.append(namedict[name])
    return header_seq_list



def AAThreeToOne(aathree):
    """Converts a three letter amino acid code into a one letter code.

    The single input argument is the three letter amino acid code.  It
    can be of any case.

    This function returns, in upper case, the one letter amino acid code.
    It raises an exception if *aathree* is not a valid one letter code.
    'Xaa' is converted to 'X'.

    >>> AAThreeToOne('Ala')
    'A'

    >>> AAThreeToOne('cys')
    'C'

    >>> AAThreeToOne('Xaa')
    'X'

    >>> AAThreeToOne('hi')
    Traceback (most recent call last):
       ...
    ValueError: Invalid amino acid code of hi.
    """
    aa_mapping = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HIS',
                'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN',
                'R':'ARG', 'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR', 'X':'XAA'}
    aa_reverse_mapping = {}
    for (one, three) in aa_mapping.iteritems():
        aa_reverse_mapping[three] = one
    try:
        return aa_reverse_mapping[aathree.upper()]
    except KeyError:
        raise ValueError("Invalid amino acid code of %s." % aathree)



def ReverseComplement(heads_seqs):
    """Converts nucleotide sequences to their reverse complements.

    The single input argument *heads_seqs* is a list of sequences
    in the format of tuples *(head, seq)* where *head* is the
    header and *seq* is the sequence.  The sequences should all
    be nucleotide sequences composed exclusively of A, T, C, or G
    (or their lowercase equivalents).  Ambiguous nucleotide codes
    are currently not accepted.

    The returned variable is a copy of *heads_seqs* in which the headers
    are unchanged but the sequences are converted to reverse complements.

    >>> ReverseComplement([('seq1', 'ATGCAA'), ('seq2', 'atgGCA')])
    [('seq1', 'TTGCAT'), ('seq2', 'TGCcat')]

    >>> ReverseComplement([('seq1', 'ATGNAA')])
    Traceback (most recent call last):
        ...
    ValueError: Invalid nucleotide code.
    """
    mapping = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
    reversecomplements = []
    for (head, seq) in heads_seqs:
        seq = list(seq)
        seq.reverse()
        try:
            rc = [mapping[nt] for nt in seq]
        except KeyError:
            raise ValueError("Invalid nucleotide code.")
        reversecomplements.append((head, ''.join(rc)))
    return reversecomplements


def FindMotifs(seq, motif):
    """Finds occurrences of a specific motif in a nucleotide sequence.

    *seq* is a string giving a nucleotide sequence.

    *motif* is a string giving the motif that we are looking for.  It should
    be a string of valid nucleotide characters:
    * A   Adenine
    * G   Guanine
    * C   Cytosine
    * T   Thymine
    * U   Uracil
    * R   Purine (A or G)
    * Y   Pyrimidine (C or T)
    * N   Any nucleotide
    * W   Weak (A or T)
    * S   Strong (G or C)
    * M   Amino (A or C)
    * K   Keto (G or T)
    * B   Not A (G or C or T)
    * H   Not G (A or C or T)
    * D   Not C (A or G or T)
    * V   Not T (A or G or C)

    The returned variable is a list *motif_indices* of the indices that
    each occurrence of *motif* in *seq* begins with.  For example,
    if there is a motif beginning at *seq[7]*, then 7 will be
    present in *motif_indices*.  So the number of occurrences of
    the motif will be equal to len(*motif_indices*).

    This function is not case sensitive, so nucleotides can be either upper
    or lower case.  In addition, T (thymine) and U (uracil) nucleotides
    are treated identically, so the function can handle either DNA
    or RNA sequences.

    >>> FindMotifs('ATCGAA', 'WCGW')
    [1]
    """
    assert isinstance(seq, str)
    assert isinstance(motif, str)
    seq = seq.upper()
    motif = motif.upper()
    seq = seq.replace('U', 'T')
    motif = motif.replace('U', 'T')
    nt_mapping = {  # maps nucleotide codes to regular expression code
        'A' : 'A',
        'G' : 'G',
        'C' : 'C',
        'T' : 'T',
        'R' : '[A,G]',
        'Y' : '[C,T]',
        'N' : '[A,T,C,G]',
        'W' : '[A,T]',
        'S' : '[G,C]',
        'M' : '[A,C]',
        'K' : '[G,T]',
        'B' : '[G,C,T]',
        'H' : '[A,C,T]',
        'D' : '[A,G,T]',
        'V' : '[A,G,C]',
    }
    try:
        m = re.compile(''.join([nt_mapping[nt] for nt in motif]))
    except KeyError:
        raise ValueError("motif contains an invalid character: %s" % motif)
    motif_indices = [match.start() for match in m.finditer(seq)]
    return motif_indices



def CondenseSeqs(seqs, maxdiffs, exclude_positions):
    """Removes nearly identical protein sequences.

    *seqs* is a list of sequences as (head, seq) 2-tuples.  The sequences
    are assumed to be aligned.

    *maxdiffs* specifies the maximum number of differences that a sequence
    can have from another sequence in order to be removed.

    *exclude_positions* is a list of integers specifying positions at which
    sequence can NOT differ and still be removed.  These integers
    are for numbering the sequences as 1, 2, ...

    The method proceeds as follows:

        1) For the first sequence iterates through the rest of the sequences,
           and removes any that differ at <= ndiff sites, and do not differ
           at the sites specified by 'exclude_position'.

        2) Repeats this process for each of the remaining sequences.

    The returned variable is *seqs* with the removed sequences gone.

    >>> seqs = [('s1', 'ATGC'), ('s2', 'ATGA'), ('s3', 'ATC-'), ('s4', 'ATGC')]
    >>> CondenseSeqs(seqs, 0, [])
    [('s1', 'ATGC'), ('s2', 'ATGA'), ('s3', 'ATC-')]
    >>> CondenseSeqs(seqs, 1, [])
    [('s1', 'ATGC'), ('s3', 'ATC-')]
    >>> CondenseSeqs(seqs, 1, [4])
    [('s1', 'ATGC'), ('s2', 'ATGA'), ('s3', 'ATC-')]

    """
    i = 0
    while i < len(seqs):
        newseqs = list(seqs[ : i + 1])
        iseq = seqs[i]
        for jseq in seqs[i + 1 : ]:
            assert len(iseq[1]) == len(jseq[1])
            differ_at_excluded = False
            for position in exclude_positions:
                if iseq[1][position - 1] != jseq[1][position - 1]:
                    differ_at_excluded = True
                    break
            if differ_at_excluded:
                newseqs.append(jseq)
            else:
                ndiffs = len([k for k in range(len(iseq[1])) if iseq[1][k] != jseq[1][k]])
                if ndiffs > maxdiffs:
                    newseqs.append(jseq)
        seqs = newseqs
        i += 1
    return seqs


def DateToOrdinal(datestring, refyear):
    """Converts a date string to an ordinal date.

    *datestring* is a date given by a string such as '2007/2/13' (for
    Feb-13-2007), or '2007/2//' if no day is specified, or
    '2007//' if no day or month is specified. The '/' characters can
    also be '-'.

    *refdate* is an integer year from the approximate timeframe we are examining
    which is used to anchor the datestring date on the assumption
    that each year has 365.25 days.

    The returned value is a number (decimal) giving the date. If no
    day is specified, the 15th (halfway through the month) is chosen.
    If no month or day is specified, July 1 (halfway through the
    year) is chosen.

    >>> print "%.2f" % DateToOrdinal('2007/4/27', 1968)
    2007.32

    >>> print "%.2f" % DateToOrdinal('2007/4/', 1968)
    2007.29

    >>> print "%.2f" % DateToOrdinal('2007//', 1968)
    2007.50

    >>> print "%.2f" % DateToOrdinal('2007-4-27', 1968)
    2007.32

    """
    if not isinstance(refyear, int):
        raise ValueError('refyear is not an integer')
    refdate = datetime.date(1968, 1, 1).toordinal()
    try:
        if '/' in datestring:
            (year, month, day) = datestring.split('/')
        else:
            (year, month, day) = datestring.split('-')
    except ValueError:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    if year and month and day:
        (year, month, day) = (int(year), int(month), int(day))
        date = datetime.date(year, month, day)
    elif year and month:
        (year, month) = (int(year), int(month))
        date = datetime.date(year, month, 15)
    elif year:
        year = int(year)
        date = datetime.date(year, 7, 1)
    else:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    return (date.toordinal() - refdate) / 365.25 + refyear


# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
