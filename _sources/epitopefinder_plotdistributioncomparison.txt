.. _epitopefinder_plotdistributioncomparison.py:

=======================================================
epitopefinder_plotdistributioncomparison.py
=======================================================

This script analyzes of output from :ref:`epitopefinder_getepitopes.py`. You can use this script if you want to see if some set of the sites in a protein contain more epitopes of some other set of the sites in the same or another protein.

The input data for this script is two or more files listing the number of epitopes per site. These input files are in the format of the *epitopesbysite* file created by :ref:`epitopefinder_getepitopes.py` or the *selectsitesfile* created by :ref:`epitopefinder_selectsites.py`. For example, you might use :ref:`epitopefinder_getepitopes.py` to define the number of epitopes per site for all sites in a protein, and then :ref:`epitopefinder_selectsites.py` to define the number of epitopes per site for selected subset of sites in the protein. You could then use this script to compare the distribution of epitope counts per site.

This script utilizes `matplotlib`_, and will fail if that package is not available for importation. 


Running the script
---------------------
:ref:`epitopefinder_plotdistributioncomparison.py` takes as input the name of a single file, the format of which is detailed below. If you have installed the package so that the scripts are the search path, you can run this script directly from the command line. For example, if you called your input file ``infile.txt`` then run::

    epitopefinder_plotdistributioncomparison.py infile.txt

If the script is not executable on your platform, then run::

    python epitopefinder_plotdistributioncomparison.py infile.txt

This will create the output plot file described below.


Input file format
---------------------
The input file is a text file that should contain the following key / value pairs. Each line begins with the key, and is followed by the value for that key. Empty lines or lines beginning with # are ignored:

    * *plotfile* is the name of the PDF plot file that is being created. This file must end with the extension ``.pdf``.

    * *epitopesfile1* specifies the name of a file containing the number of epitopes per site. This should be a CSV file matching the format of the *epitopesbysite* file created by :ref:`epitopefinder_getepitopes.py` or the *selectsitesfile* created by :ref:`epitopefinder_selectsites.py` . After an initial header line, each line lists the site number followed by the number of epitopes at that site. For example::

        Site,NumberUniqueEpitopes
        1,0
        2,1
        3,1
        4,3
        5,3
        6,3

    * *epitopesfile2* specifies the name of a second file containing the number of epitopes per site in the same format as *epitopesfile1*.


    * *set1name* is a string (LaTex formatting) that is used to label the distribution of the number of epitopes per site in *epitopesfile1*.

    * *set2name* is a string (LaTex formatting) that is used to label the distribution of the number of epitopes per site in *epitopesfile2*.

    * *pvalue* is an option that should only be used if *epitopesfile2* specifies a subset of the sites found in *epitopesfile1*. It allows you to compute the P-value that the subset of sites in *epitopesfile2* has a higher (or lower) average number of epitopes per site than the full set of sites in *epitopesfile1*. If you do not want to use this option, set *pvalue* to *None*.
    
      Otherwise, *pvalue* is implemented as follows. We draw random subsets of the same number of sites that are found in *epitopesfile2* from the full set of sites in *epitopesfile1*. The number of such random subsets that are drawn is the integer specified by *pvalue*, so a reasonable value of *pvalue* is 100000. We then compute the average number of epitopes per site for the random subsets and compare it to the average number of epitopes per site in the actual subset in *epitopesfile2*. If more than half of the random subsets have fewer average epitopes than the actual subset in *epitopesfile2*, then the plot reports the P-value that *epitopesfile2* has more epitopes than a random subset. If less than half of the random subsets have fewer average epitopes than the actual subset, then the plot reports the P-value that *epitopesfile2* has fewer epitopes than a random subset. So the computed value is a **one-sided** P-value for the hypothesis that the mean number of epitopes per site for the subset of sites in *epitopesfile2* is < or > than the mean number of epitopes per site expected for a random subset of that many sites from *epitopesfile1*. The hypothesis that is being test (< or >) is indicated on the plot. 

    * *pvaluewithreplacement* is an option that is only required if *pvalue* is being used. In this case, *pvaluewithreplacment* should be either *True* or *False*. If it is *True*, then the random subsets are drawn with replacement (so the same site can be drawn more than once into a random subset). If it is *False*, then the random subsets are drawn without replacement (so the same site can only be drawn once into a random subset). You will need to figure out which is more appropriate. In general, if you have created the subset in *epitopesfile2* using :ref:`epitopefinder_selectsites.py` with *retainmultiple* set to *True*, then you will want to make *pvaluewithreplacement* also *True*. If you have created the subset in *epitopesfile2* using :ref:`epitopefinder_selectsites.py` with *retainmultiple* set to *False*, then you will want to make *pvaluewithreplacement* also *False*.

      If *pvalue* is *None*, then you don't need to specify any key for *pvaluewithreplacement*. If you do specify *pvaluewithreplacement* when *pvalue* is *None*, it has no meaning.

    * *title* specifies the title placed above the plot. It can be *None* if no title is to be used. Otherwise, it should be the title (using LaTex formatting, spaces are allowed) that is placed above the plot.

    * *ymax* is an optional key that specifies that we fix the maximum value of the y-axis. You can simply leave out this option or set it to *None* if you do not want to fix the y-axis. Otherwise set it to the number that you would like to make the y-maximum.


Example input file
---------------------
Here is an example input file::

    # input file for epitopefinder_plotdistributioncomparison.py
    plotfile distributioncomparison_all_vs_subset.pdf
    epitopesfile1 epitopesbysite.csv
    epitopesfile2 selectedsites.csv
    set1name all sites
    set2name selected sites
    pvalue 100000
    pvaluewithreplacement True
    title None


Output files
-------------------
This script creates the PDF plot file *plotfile*, which shows the distribution of number of epitopes per site for *epitopesfile1* and *epitopesfile2*. This plot may also display a P-value depending on the setting for *pvalue*.

Here is an example plot.

.. figure:: distributioncomparison_example.jpg
   :width: 45%
   :align: center
   :alt: distributioncomparison_example.jpg


.. include:: weblinks.txt
