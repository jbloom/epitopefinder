.. _epitopefinder_plotlineardensity.py:

================================================
epitopefinder_plotlineardensity.py
================================================
This script helps with analysis of output from :ref:`epitopefinder_getepitopes.py`. Specifically, it plots the density of epitopes (or some other variable) as a function of the protein primary sequence position. Most commonly, you might use it to visualize the results contained in the *epitopesbysite* file(s) created by :ref:`epitopefinder_getepitopes.py`.

This script utilizes `matplotlib`_, and will fail if that package is not available for importation.

Running the script
---------------------
:ref:`epitopefinder_plotlineardensity.py` takes as input the name of a single file, the format of which is detailed below. If you have installed the package so that the scripts are the search path, you can run this script directly from the command line. For example, if you called your input file ``infile.txt`` then run::

    epitopefinder_plotlineardensity.py infile.txt

If the script is not executable on your platform, then run::

    python epitopefinder_plotlineardensity.py infile.txt

This will create the output described below.


Input file format
---------------------
The input file is a text file that should contain the following lines in the order indicated below. Empty lines or lines beginning with # are ignored:

    * *plotfile* : This line should contain the word *plotfile* followed by a string giving the name of the PDF plot file that is being created. This file must end with the extension ``.pdf``.

    * *title* specifies the title placed above the plot. It can be *False* if no title is to be used. Otherwise, it should be the title (using LaTex formatting, spaces are allowed) that is placed above the plot.

    * *fixymax* specifies that we fix the y-maximum to the specified value.
      This may be useful if you are making multiple plots for comparisons
      between them, and want them all to have the same y-maximum.
      Note that the value specified here is taken to be the data
      maximum -- the actually maximum of the y-axis is somewhat
      higher to provide some padding space. If you do not want to fix the
      y-maximum then set this option to the string *False*. Otherwise
      set it to the number to which you want to fix this maximum.

    * A listing of the data files. Each listing should be a line containing two entries: *datafile* and *label*. The *datafile* entry should be the name of a CSV data file; this file name cannot contain any spaces and should end in the extension ``.csv``. The *label* entry is the name of the label given to that data series in the plot. These labels use LaTex string formatting, and are allowed to contain spaces. Each *datafile* should list the *x* (usually the residue number) and *y* (usually the number of epitopes for that site) in columns separated by a comma. Both entries must be numbers. The first line of the file is always ignored as a title line, as are any additional lines beginning with a # character. Here is an example::

        Site,NumberUniqueEpitopes
        1,0
        2,1
        3,1
        4,3
        5,3
        6,3


Example input file
---------------------
Here is an example input file::

    # input file for epitopefinder_plotlineardensity.py
    plotfile epitopelineardensity.pdf
    title CTL epitopes in human H3N2 influenza NP
    fixymax False
    epitopesbysite.csv MHC class I epitopes



Output files
-------------------
This script creates the PDF plot file *plotfile*. This plot uses lines to plot the variables specified by the data file as a function of the primary sequence position. Most commonly, this would be the number of epitopes as a function of primary sequence.

Here is an example of such a plot.

.. figure:: epitopelineardensity.jpg
   :width: 90%
   :align: center
   :alt: epitopelineardensity.jpg

.. include:: weblinks.txt
