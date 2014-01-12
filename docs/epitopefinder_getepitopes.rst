.. _epitopefinder_getepitopes.py:

================================================
epitopefinder_getepitopes.py
================================================
This is a main script of the `epitopefinder`_ package. 

This script answers the question: How many epitopes of a specified type map to each position in a protein's sequence? A user created input file specifies parameters used to filter the epitopes, map them to a protein sequence, and determine whether different mapped epitopes are unique.

The script begins by reading collection of epitopes from a CSV file *iedbfile* downloaded from the `Immune Epitope Database`_, and then classifying the epitopes as best as possible by MHC class and MHC gene. MHC class II epitopes are also classified by supertype using the classification system in *supertypesfile*.

The epitopes are then further filtered by MHC class (according to the *mhcclass* parameter) and by epitopelength (according to the *epitopelength* parameter) to give a set of epitopes under consideration. We then attempt to align these epitopes to the target protein sequence (or collection aligned homologous target protein sequences) specified in *targetprotsfile*. Epitopes that align to at least one of these target protein sequences with no gaps and <= *maxmismatches* mismatches are considered to match the target protein(s). The alignment is performed using the MUSCLE program in *musclepath*. Finally, potentially redundant epitopes are removed according to the options set by *purgeredundant* and *purgeredundantoverlap*.

The results are written to the output CSV (comma-separated value) files *epitopelistfile* and *epitopesbysitefile*.

Running the script
---------------------
:ref:`epitopefinder_getepitopes.py` takes as input the name of a single file, the format of which is detailed below. If you have installed the package so that the scripts are the search path, you can run this script directly from the command line. For example, if you called your input file ``infile.txt`` then run::

    epitopefinder_getepitopes.py infile.txt

If the script is not executable on your platform, then run::

    python epitopefinder_getepitopes.py infile.txt

This will create the various output files described below.


Input file format
---------------------
The input file should contain a series of keys followed by their values. Lines beginning with # or empty lines are ignored. The input file must have the following keys and their values:

    * *iedbfile* is a file containing epitopes in the format of downloads from the `Immune Epitope Database`_. It is important that you download only epitopes for which the assay result is positive, otherwise an exception will be raised. These files are downloadable from this database after you do a search for some specific class of epitopes. You will be given an option to download a *full* or *compact* file with the epitopes. This script is only confirmed to work with the *compact* formats. It is possible that it will break if the format of these downloads changes; however, it works with the download format as of April-15-2013. The files themselves are CSV (comma-separated values) files with an initial header line.

    * *supertypesfile* is a file that gives the supertype classification of HLA-A and HLA-B alleles. The suggested classification scheme is that given in `BMC Immunology, 2008, 9(1)`_, although you may want to update this if the scheme is updated further. A file containing the classification scheme from this reference is contained in some of the example analyses in the ``examples`` directories in this package. This file should contains lines that list the alleles at the ``gene*group:protein`` level of detail, followed by the supertype, as in::

        A*01:01 A01
        A*01:02 Unclassified
        A*68:15 A02
        B*07:43 B07
        B*08:01 B08

    * *mhcclass* specifies the class of MHC epitopes for which we are searching. It can be either I or II depending on whether we searching for class I or class II epitopes.

    * *epitopelength* specifies the length of epitopes for which we include. It should be a list of two numbers, first the minimum length and then the maximum length. We only retain and consider epitopes with lengths >= the minimum length and <= the maximum length.

    * *targetprotsfile* specifies a FASTA file containing the target protein sequence(s) to which we are mapping the epitope. It should contain one or more sequences. If there are multiple sequences, they should be aligned. The reason that you might want to include multiple homologous sequences is if the epitopes were defined on a set of closely homologous sequences. For instance, epitopes in influenza might be determined for several closely related sequences, and you might want to find epitopes that closely match any of these sequences. If there is more than one sequence in *targetprotsfile*, then these sequences should be aligned. The sequence positions for the epitope alignments are then numbered according to consecutive numbering of the first sequence. Here is an example *targetprotsfile*::

        >A/Aichi/2/1968_NP
        MAS-GTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDY
        >A/Brisbane/10/2007_NP
        MASQGTKRSYEQMETDGDRQNATEIRASVGKMIDGIGRFYIQMCTELKLSDH
    
    In this example, the residues are numbered sequentially according to the first sequence, so that the M is 1, the A is 2, the S is 3, the G is 5 (since 4 is a gap), etc. 

    * *maxmismatches* specifies the stringency of the epitope alignment with the target protein sequences in *targetprotsfile*. We only retain epitopes that align to at least one of the sequences in *targetprotsfile* with <= *maxmismatches*.

    * *musclepath* specifies the path to a directory that contains the `MUSCLE`_ alignment program in the form of an executable named ``muscle``. This script has been tested with `MUSCLE`_ version 3.8, but should also work with other similar versions. This alignment program is used to align the epitopes to the target proteins in *targetprotsfile*.

    * *purgeredundant* specifies how we purge redundant epitopes (i.e. epitopes listed in *iedbfile* that are likely to actually refer to the same epitope). We find epitopes that overlap by >= *purgeredundantoverlap* sites in the position to which they align to the target proteins. Depending on the option set for *purgeredundant*, such overlapping epitopes are removed. The epitope that is retained among overlapping epitopes is the one with the shortest epitope sequence, and the purging moves through epitopes from shortest to longest (ensuring that distinct short epitopes that overlap with the same longer epitope are not considered redundant). Possibilities for *purgeredundant* are:

        - *None* specifies that no purging of redundant epitopes is done.

        - *MHCgene* specifies that we retain overlapping epitopes if they differ in their MHC gene. So for example two overlapping epitopes that are both HLA-A would have one purged, but two overlapping epitopes that are HLA-A and HLA-B would not have one purged.

        - *MHCsupertype* is an allowable option only if *MHCclass* was set to I. In this case, we retain overlapping epitopes if they differ in their MHC gene or if they differ in their MHC supertype.

        - *MHCgroup* specifies that we retain overlapping epitopes if they differ in their MHC gene, MHC supertype (if class I), or MHC allele group. Recall that according to `HLA nomenclature`_, the allele group is the classification after the gene, so for example HLA-B*35:01 has a group of 35.

     If two epitopes overlap and one does not have the lowest level of MHC gene / supertype / group set, then it is not retained as unique. So for example, if the same epitope is classified as HLA-B and HLA-B*35:01 and *purgeredundant* is set to *MHCgroup*, then one of the epitopes is purged because there is not enough information to ensure that they do not overlap in group.

    * *purgeredundantoverlap* specifies the amount of overlap two epitopes must possess in the positions in which they align to the target gene before they are considered potentially redundant, and candidates for removal according to the options specified by *purgeredundant*. It should be a number, and epitopes with overlap >= this number and candidates for removal.

   * *epitopelistfile* is the name of an output CSV file created by this script. It lists each of the unique epitopes identified, along with the sites where it aligns, its allele, its references, and other redundant epitopes. If this file already exists, it is overwritten.

   * *epitopesbysitefile* is the name of an output CSV file created by this script. For each site in the target protein, reports the number of unique epitopes that span this site. If this file already exits, it is overwritten.


Example input file
---------------------
Here is an example input file for MHC class I epitopes::

    # input file for epitopefinder_getepitopes.py
    iedbfile Grant2013_plus_IEDB_Influenza_Tcell_compact_2013-04-15.csv
    supertypesfile supertype_classification.txt
    mhcclass I
    epitopelength 8 12
    targetprotsfile targetprotsfile.fasta
    maxmismatches 1
    musclepath /Users/jbloom/muscle3.8/
    purgeredundant MHCsupertype
    purgeredundantoverlap 8
    epitopelistfile epitopeslist.csv
    epitopesbysitefile epitopesbysite.csv


Output files
-------------------
This script creates the two output CSV files *epitopeslistfile* and *epitopesbysitefile*. 

*epitopeslistfile* lists each unique epitope along with its redundant counterparts. Here is an example of the first two lines::

   unique_epitope_sequence,unique_epitope_position,unique_epitope_allele,unique_epitope_information,unique_epitope_reference,redundant-1_epitope_sequence,redundant-1_epitope_position,redundant-1_epitope_allele,redundant-1_epitope_information,redundant-1_epitope_reference
   "GERQNATEI","17-25","allele HLA-B*18:01; MHC class I; MHC gene B; MHC supertype B44; MHC group 18","host = Homo sapiens; source = Influenza A virus; molecule = nucleoprotein; reported position = (17, 25); assay = ELISPOT / cytokine release IFNg","Jeff Alexander; Pamuk Bilsel; Marie-France del Guercio; Aleksandra Marinkovic-Petrovic; Scott Southwood; Stephani Stewart; Glenn Ishioka; Maya F Kotturi; Jason Botten; John Sidney; Mark Newman; Alessandro Sette. Hum Immunol. 2010. PMID 20156506",,,,,
   
*epitopesbysitefile* lists all of the sites as they are numbered in sequential 1, 2, ... numbering of the aligned target proteins in *targetprotsfile*, followed by the number of unique epitopes containing that site. Here is an example of the first few lines::

    Site,NumberUniqueEpitopes
    1,0
    2,0
    3,0
    4,0
    5,0
    6,0
    7,0
    8,0
    9,0
    10,0
    11,0
    12,0
    13,0
    14,0
    15,0
    16,0
    17,3
    18,3
    19,3
    20,3

.. include:: weblinks.txt
