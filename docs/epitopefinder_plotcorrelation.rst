.. _epitopefinder_plotcorrelation.py:

================================================
epitopefinder_plotcorrelation.py
================================================
This script helps with analysis of output from :ref:`epitopefinder_getepitopes.py`. This script plots the correlation between two variables. For instance, you might want to use it to plot the correlation between the number of epitopes per site (as contained in the *epitopesbysite* file created by :ref:`epitopefinder_getepitopes.py` and some other per-site property.

The plot created is a scatter plot, with one data set on the x-axis and the other data set on the y-axis. The correlation coefficient can optionally also be displayed. Essentially, you provide the script with two input files giving some property (such as epitopes per site) for all or some sites in a protein. For all sites for which both input files specify a property value, a point is plotted on the correlation plot.

This script utilizes `matplotlib`_, and will fail if that package is not available for importation. If you select the option for displaying correlation coefficients, the script also utilizes `scipy`_ and so will fail if that package is not available for importation.


Running the script
---------------------
``epitopefinder_plotcorrelation.py`` takes as input the name of a single file, the format of which is detailed below. If you have installed the package so that the scripts are the search path, you can run this script directly from the command line. For example, if you called your input file ``infile.txt`` then run::

    epitopefinder_plotcorrelation.py infile.txt

If the script is not executable on your platform, then run::

    python epitopefinder_plotcorrelation.py infile.txt

This will create the output described below.


Input file format
---------------------
The input file is a text file that should contain the following key / value pairs. Each line begins with the key, and is followed by the value for that key. Empty lines or lines beginning with # are ignored:

    * *plotfile* : This line should contain the word *plotfile* followed by a string giving the name of the PDF plot file that is being created. This file must end with the extension ``.pdf``.

    * *title* specifies the title placed above the plot. It can be *False* if no title is to be used. Otherwise, it should be the title (using LaTex formatting, spaces are allowed) that is placed above the plot.

    * *correlation* specifies whether we compute and display the correlation coefficient between the two data sets. It can be *None* if you do not want to display this correlation coefficient. Otherwise, it should be either the string *Pearson* or the string *Spearman* depending on whether you want to compute Pearson's parametric correlation coefficient or Spearman's non-parametric correlation coefficient.

    * *xdatafile* specifies the name of an existing CSV file that contains one of the data values for sites in the protein. This value is plotted on the x-axis of the correlation plot. The first line of the file is assumed to be a header and is ignored, as are any lines beginning with # or empty lines. All other lines should contain two numbers: the site (residue number) and the value associated with that number. Here is an example::

        Site,NumberUniqueEpitopes
        1,0
        2,1
        3,1
        4,3
        5,3
        6,3

    * *ydatafile* is like *xdatafile*, but specifies the second value for the sites (which is plotted on the y-axis). Note that not all the same site numbers need to be present in both *xdatafile* and *ydatafile*, but only data pairs for which values are present in both files are plotted. There must be at least two data pairs or the script will raise an exception.

    * *xlabel* specifies a string (in LaTex format) that is placed on the x-axis.

    * *ylabel* specifies a string (in LaTex format) that is placed on the y-axis.

Example input file
---------------------
Here is an example input file::

    # input file for epitopefinder_plotcorrelation.py
    plotfile epitopestabilitycorrelation.pdf
    title False
    correlation Spearman
    xdatafile mutation_dtms.csv
    ydatafile epitopesbysite.csv
    xlabel $\Delta T_m$
    ylabel number of MHC class I epitopes



Output files
-------------------
This script creates the PDF plot file *plotfile*. This plot uses points to show the correlation between the data values in *xdatafile* and *ydatafile*. 

.. include:: weblinks.txt
