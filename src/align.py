"""This module contains functions to run sequence alignment programs.

It can be used to run PROBCONS or MUSCLE.  

Written by Jesse Bloom."""


import os
import tempfile
import sys
import subprocess
import sequtils


def GetEpitopeAlignments(epitope, targetseqs, maxmismatches, musclepath):
    """Finds alignments of an epitope to one or more target sequences.

    *epitope* is a string giving the sequence of an epitope sequence.

    *targetseqs* is a list of one more *(header, sequence)* tuples.
    These tuples give the target protein sequences to which we 
    attempt to align *epitope*. If there are multiple sequences,
    they must be aligned (i.e. all of the same length).

    *maxmismatches* is the maximum number of mismatches that are
    allowed for an epitope to still be considered a match.

    *musclepath* is the path to a directory containing the
    MUSCLE alignment program.

    This function uses MUSCLE to align *epitope* to each of the
    sequences in *targetseq*. If it aligns with <= *maxmismatches*
    mismatches and with no gaps, then this is considered a valid
    alignment. The possible outcomes are:

        * A return value of *False* if *epitope* does not align
          to any of the sequences in *targetseqs* with the
          specified identity and no gaps.

        * A 3-tuple of the form *(alignstart, alignend, heads)*
          if *epitope* aligns to one or more of the sequences
          in *targetseqs*. In this case, *epitope* must align
          to the same position in all of the target sequences
          to which it aligns. *alignstart* is the index of
          the first position in the target sequence in which
          this alignment starts (in 1, 2, ... numbering of
          the sequence as it appears in *targetseq*), *alignend*
          is the index of the last position in the target 
          sequence (1, 2, ... numbering), and *heads* is a list
          of the headers for all sequences in *targetseqs* to
          which a successful alignment at these indices was found.

        * An Exception is raised if *epitope* aligns to multiple
          sequences in *targetseqs* but with different indices.
    """
    matchlist = []
    n = len(targetseqs[0][1])
    epitope = epitope.upper()
    for (head, seq) in targetseqs:
        seq = seq.upper()
        if len(seq) != n:
            raise ValueError('Sequences in targetseqs differ in length')
        a = Align([('epitope', epitope), (head, seq)], musclepath, 'MUSCLE')
        strippeda = StripLeadingTrailingGapsToFirstSequence(a)
        if ('-' in strippeda[0][1]) or ('-' in strippeda[1][1]):
            continue # gap in alignment
        assert len(epitope) == len(strippeda[0][1]) == len(strippeda[1][1])
        ndiffs = len([i for i in range(len(epitope)) if strippeda[0][1][i] != strippeda[1][1][i]])
        if ndiffs > maxmismatches:
            continue
        # now find alignstart and alignend. We use this procedure 
        # because MUSCLE will strip gaps from the aligned sequence
        alignstart = alignend = ialignment = -1
        alignedepitope = a[0][1]
        alignedseq = a[1][1]
        for i in range(n):
            if seq[i] != '-':
                ialignment += 1
                assert ialignment == -1 or seq[i] == alignedseq[ialignment]
            if ialignment != -1:
                if (alignedepitope[ialignment] != '-') and alignstart == -1:
                    alignstart = i
                if (alignedepitope[ialignment] == '-') and alignstart != -1:
                    alignend = i
                    break
        else:
            alignend = min(i + 1, len(seq))
        assert len(epitope) == alignend - alignstart, "Here is alignment:\n>%s\n%s\n>%s\n%s\nlen(epitope) = %d, alignend = %d, alignstart = %d" % (a[0][0], a[0][1], a[1][0], a[1][1], len(epitope), alignend, alignstart)
        extractedmatch = seq[alignstart : alignend].replace('-', '')
        assert len(extractedmatch) == len(epitope)
        diffs = len([j for j in range(len(epitope)) if epitope[j] != extractedmatch[j]])
        assert diffs <= maxmismatches
        matchlist.append((alignstart + 1, alignend, head))
    if matchlist:
        (alignstart, alignend, head) = matchlist[0]
        for (x, y, z) in matchlist[1 : ]:
            if (x, y) != (alignstart, alignend):
                raise ValueError("Epitope %s aligned to different locations in %s and %s" % (epitope, head, z))
        return matchlist
    else:
        return False



def Align(headers_seqs, progpath, program='MUSCLE'):
    """Performs a multiple sequence alignment of two or more sequences.

    By default, the protein sequences are aligned using MUSCLE. This
    program can be used to align either nucleotide or protein
    sequences. You can also use PROBCONS to align protein sequences.

    *headers_seqs* is a list specifying the names of the sequences that we
    want to align.  Each entry is a 2-tuple *(head, seq)* where *head* is
    a header giving the sequence name and other information (might be empty)
    and *seq* is a string giving the protein sequence.  The list must have
    at least 2 entries.

    *progpath* specifies a directory containing the alignment program executable,
    either PROBCONS or MUSCLE.  The PROBCONS executable is assumed to have
    the name "probcons" in this directory.  The MUSCLE executable is assumed to
    have the name "muscle" in this directory.

    *program* specifies what program to use for the alignment. By default, it is
    "MUSCLE".  If you wish to use PROBCONS instead, set it to "PROBCONS".  

    This executable is used to perform a multiple sequence alignment 
    with the default settings of PROBCONS or MUSCLE.  The returned variable is a
    new list *aligned_headers_seqs*. Each entry is a 2-tuple *(head, aligned_seq)*.  
    *head* has the same meaning as on input (the sequence header) and 
    *aligned_seq* is the aligned sequence, with gaps inserted as '-'
    as appropriate.  Therefore, all of the *aligned_seq* entries in
    *aligned_headers_seqs* are the same length. Entries in *aligned_headers_seq*
    are in the same order as in the input list *headers_seqs*.
    """
    if not (isinstance(headers_seqs, list) and len(headers_seqs) >= 2):
        raise ValueError, 'header_seqs does not specify a list with at least two entries.'
    if not os.path.isdir(progpath):
        raise ValueError, "Cannot find directory %s." % progpath
    if program == 'PROBCONS':
        exe = os.path.abspath("%s/probcons" % progpath) # the executable
    elif program == 'MUSCLE':
        exe = os.path.abspath("%s/muscle" % progpath) # the executable
    else:
        raise ValueError, "Invalid value of %s for 'program'." % (str(program))
    if not os.path.isfile(exe):
        raise IOError, "Cannot find executable at %s." % exe
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in a temporary directory
        infile = "%s/in.fasta" % tempdir # input file 
        outfile = "%s/out.fasta" % tempdir # output file 
        sequtils.Write(headers_seqs, infile) # write sequences to the input file
        if program == 'PROBCONS':
            p = subprocess.Popen("%s %s" % (exe, infile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run ProbCons
            (output, errors) = p.communicate()
            open(outfile, 'w').write(output)
        elif program == 'MUSCLE':
            p = subprocess.Popen("%s -in %s -out %s" % (exe, infile, outfile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run MUSCLE
            (output, errors) = p.communicate()
        try:
            aligned_headers_seqs = sequtils.Read(outfile)
        except:
            sys.stderr.write("Error getting alignment output, error of %s" % errors)
            raise
    finally:
        os.chdir(currdir) # return to the original directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file)) # remove files from temporary directory
        os.rmdir(tempdir) # remove temporary directory
    if len(aligned_headers_seqs) != len(headers_seqs):
        raise ValueError, "Did not return the correct number of aligned sequences."
    # put the aligned sequences in the same order as the input sequences
    n = len(aligned_headers_seqs[0][1]) # length of aligned sequences
    d = dict(aligned_headers_seqs)
    aligned_headers_seqs = []
    for (head, seq) in headers_seqs:
        try:
            alignedseq = d[head]
        except KeyError:
            raise ValueError, "After alignment, the following header is missing: %s" % head
        if len(alignedseq) != n:
            open('errors.temp', 'w').write(errors)
            raise ValueError, "Aligned sequence %s is not of length %d: if you are using MUSCLE, you may be running out of memory.  Errors have been written to errors.temp." % (alignedseq, n)
#        if len(seq) - seq.count('-') > n:
#            open('errors.temp', 'w').write(errors)
#            raise ValueError, "Unaligned seq %s is longer than aligned length of %d: if you are using MUSCLE, you many be running out of memory.  Errors have been written to errors.temp." % (seq, n)
        aligned_headers_seqs.append((head, alignedseq))
    return aligned_headers_seqs # return the aligned sequences


def PairwiseStatistics(aligned_headers_seqs):
    """Computes the number of gaps and identities in a pairwise alignment.

    This method is designed to compute statistics about an alignment of two
        sequences.  

    *aligned_headers_seqs* is a pair of aligned sequences, as would be returned by
    calling `Align` with two sequences (which should be of the same 
    length since they have been aligned).  That is, it is a list of 2-tuples::

            [(head1, alignedseq1), (head2, alignedseq2)]

    The method returns the 2-tuple *(identities, gaps)*.  *identities* is a number
    between zero and one.  It is the fraction of residues in one
    sequence that are aligned with identical residues in the other sequence, gaps
    not being included in the tally.  This is computed by dividing the number of
    identities by the total length of the aligned sequence excluding gaps.
    *gaps* is the fraction of gaps in the 
    alignment.  It is the fraction of the positions in the alignment where
    either sequence has a gap.  So it is computed by dividing the total number
    of gaps by the length of the aligned sequences.  Upper and lower case 
    nucleotides are treated equivalently.

    >>> print round(PairwiseStatistics([('seq1', 'TGCAT'), ('seq2', 'AG-AT')])[0], 2)
    0.75
    
    >>> print round(PairwiseStatistics([('seq1', 'tgcat'), ('seq2', 'AG-AT')])[1], 2)
    0.2
    """
    try:
        [(head1, alignedseq1), (head2, alignedseq2)] = aligned_headers_seqs
    except ValueError:
        raise ValueError, "Invalid input to 'PairwiseStatistics'."
    n = len(alignedseq1)
    if not n:
        raise ValueError, "Empty sequences of zero length."
    if n != len(alignedseq2):
        raise ValueError, "Aligned sequences are not of the same length."
    ngaps = nidentities = 0
    for (r1, r2) in zip(alignedseq1.upper(), alignedseq2.upper()):
        if r1 == '-' == r2:
            raise ValueError, "In a proper pairwise alignment, both residues should not be gaps."
        elif r1 == '-' or r2 == '-':
            ngaps += 1
        elif r1 == r2:
            nidentities += 1
    return (nidentities / float(n - ngaps), float(ngaps) / n)



def StripLeadingTrailingGapsToFirstSequence(aligned_headers_seqs):
    """Strips leading / trailing gaps from first sequence, and trims corresponding alignments.

    On input, *aligned_headers_seqs* is a set of two or more aligned sequences,
    as would be returned by `Align`.

    The first sequence in the alignment corresponds to the reference sequence.
    The returned variable is a list similar to *aligned_headers_seqs*, but all
    leading / trailing gaps have been stripped from the reference sequence.  A leading
    gap ('-') is one that precedes the first non-gap character; a trailing gap is
    one that follows the last non-gap character.  All other sequences have all
    positions corresponding the leading/trailing gaps of the reference sequence
    trimmed as well.  The headers and order of sequences are preserved.

    >>> StripLeadingTrailingGapsToFirstSequence([('s1', '--ATA-GC-'), ('s2', 'TGATTA-CA')])
    [('s1', 'ATA-GC'), ('s2', 'ATTA-C')]
    """
    if not (isinstance(aligned_headers_seqs, list) and len(aligned_headers_seqs) >= 2):
        raise ValueError, "aligned_headers_seqs does not specify at least two aligned sequences."
    (ref_head, ref_seq) = aligned_headers_seqs[0]
    firstnongap = 0
    while ref_seq[firstnongap] == '-':
        firstnongap += 1
    lastnongap = len(ref_seq)
    while ref_seq[lastnongap - 1] == '-':
        lastnongap -= 1
    non_strip_positions = [i for i in range(firstnongap, lastnongap)] # positions not to strip away
    stripped_headers_seqs = []
    for (ihead, iseq) in aligned_headers_seqs:
        istrippedseq = ''.join([iseq[i] for i in non_strip_positions])
        stripped_headers_seqs.append((ihead, istrippedseq))
    return stripped_headers_seqs




def StripGapsToFirstSequence(aligned_headers_seqs):
    """Strips gaps from a reference sequence, and all corresponding alignments.

    On input, *aligned_headers_seqs* should be a set of two or more aligned sequences,
    as would be returned by `Align`.

    The first sequence in this alignment is taken to correspond to the reference sequence.
    The returned variable is a list similar to aligned_headers_seqs, but with
    all positions corresponding to gaps in this reference sequence stripped away.
    All gaps ('-') characters are removed from this reference sequence.  In addition,
    in all other aligned sequences in *aligned_headers_seqs*, every character at
    the same position as a gap in the reference sequence is removed.  Therefore,
    at the end of this procedure, all of the alignments have the same length
    as the reference sequence with its gaps stripped away.  The headers are 
    unchanged.  The order of sequences in this stripped alignment is also
    unchanged.

    >>> StripGapsToFirstSequence([('s1', '-AT-A-GC'), ('s2', 'AAT-TAGC'), ('s3', '--T-A-GC')])
    [('s1', 'ATAGC'), ('s2', 'ATTGC'), ('s3', '-TAGC')]
    """        
    if not (isinstance(aligned_headers_seqs, list) and len(aligned_headers_seqs) >= 2):
        raise ValueError, "aligned_headers_seqs does not specify at least two aligned sequences."
    (ref_head, ref_seq) = aligned_headers_seqs[0]
    non_strip_positions = [] # positions not to strip away
    stripped_ref_seq = []
    for i in range(len(ref_seq)):
        r = ref_seq[i]
        if r != '-':
            non_strip_positions.append(i)
            stripped_ref_seq.append(r)
    stripped_headers_seqs = [(ref_head, ''.join(stripped_ref_seq))]
    for (ihead, iseq) in aligned_headers_seqs[1 : ]:
        istrippedseq = ''.join([iseq[i] for i in non_strip_positions])
        stripped_headers_seqs.append((ihead, istrippedseq))
    return stripped_headers_seqs


def AddDots(aligned_headers_seqs):
    """Adds dots at identities in multiple sequence alignment.

    Takes as an argument a list of two or more aligned sequences, as would be
    returned by `Align`.

    Returns a copy of this list.  The first sequence in the list is unchanged.
    In all remaining sequences, dot characters (".") have been used to
    replace any amino acids that are identical to the amino acid at the 
    same position in the first sequence, except for gap characters.
    Characters are replaced even if they are not of the same case

    >>> AddDots([('seq1', '-TGC'), ('seq2', 'AGGC'), ('seq3', '-tac')])
    [('seq1', '-TGC'), ('seq2', 'AG..'), ('seq3', '-.a.')]
    """
    if not (isinstance(aligned_headers_seqs, list) and len(aligned_headers_seqs) >= 2):
        raise ValueError, "aligned_headers_seqs does not specify at least two aligned sequences."
    (ref_head, ref_seq) = aligned_headers_seqs[0]
    dotted_headers_seqs = [(ref_head, ref_seq)]
    for (head, seq) in aligned_headers_seqs[1 : ]:
        assert len(seq) == len(ref_seq)
        seq = list(seq)
        for r in range(len(seq)):
            if ref_seq[r].upper() == seq[r].upper() and ref_seq[r] != '-':
                seq[r] = '.'
        seq = ''.join(seq)
        dotted_headers_seqs.append((head, seq))
    return dotted_headers_seqs


def RemoveDots(aligned_headers_seqs):
    """Removes dots at identities in a multiple sequence alignment.

    This function effectively undoes what can be done by 'AddDots'.

    Takes as an argument a list of two or more aligned sequences
    in the form of (header, sequence) tuples, as would be returned
    by **Align**.

    Returns a copy of this list.  Any positions where one of the sequences
    after the first sequence has a '.' character, the amino acid
    is changed to that found at the same position in the first sequence.

    >>> RemoveDots([('seq1', '-TGC'), ('seq2', 'AG..'), ('seq3', '-.a.')])
    [('seq1', '-TGC'), ('seq2', 'AGGC'), ('seq3', '-TaC')]
    """
    if not (isinstance(aligned_headers_seqs, list) and len(aligned_headers_seqs) >= 2):
        raise ValueError, "aligned_headers_seqs does not specify at least two aligned sequences."
    (ref_head, ref_seq) = aligned_headers_seqs[0]
    undotted_headers_seqs = [(ref_head, ref_seq)]
    for (head, seq) in aligned_headers_seqs[1 : ]:
        assert len(seq) == len(ref_seq)
        seq = list(seq)
        for r in range(len(seq)):
            if seq[r] == '.':
                seq[r] = ref_seq[r]
        seq = ''.join(seq)
        undotted_headers_seqs.append((head, seq))
    return undotted_headers_seqs


# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
