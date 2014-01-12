"""Module for performing plotting.

This module uses ``pylab`` and ``matplotlib`` to make plots. 
Before running a function in this module, you should use the *PylabAvailable*
function to determine if ``pylab`` and ``matplotlib`` are available. Otherwise,
calling any other function will raise an Exception if thise modules are
not available. The ``pdf`` backend is used for ``matplotlib`` / ``pylab``. 
This means that plots must be created as PDF files.

A few functions also utilize ``scipy`` for calculations. Before using these
functions, you should use *ScipyAvailable* to see if ``scipy`` is available.
Otherwise an exception will be raised.



List of functions
--------------------

`PylabAvailable`

`CumulativeFractionPlot`

`Base10Formatter`

`SubsetPValue`

`SplitLabel`

`PlotLinearDensity`

`CorrelationPlot`

`PlotDistributionComparison`


Details of functions
----------------------
Provided in their individual documentation strings below.
"""


import os
import time
import sys
import math
import random


# global variable _pylabavailable indicates if pylab/matplotlib present
try:
    import matplotlib
    matplotlib.use('pdf')
    import pylab
    _pylabavailable = True
except ImportError:
    _pylabavailable = False

# global variable _scipyavailable indicates if scipy is present
try:
    import scipy.stats
    _scipyavailable = True
except ImportError:
    _scipyavailable = False


def PylabAvailable():
    """Returns True if pylab/matplotlib available, False otherwise.
    
    You should call this function to test for the availability of the
    pylab/matplotlib plotting modules before using other functions in
    this module.
    """
    return _pylabavailable


def ScipyAvailable():
    """Returns *True* if scipy is available, *False* otherwise."""
    return _scipyavailable


def SubsetPValue(subset, fullset, nrandom, withreplacement):
    """Computes P-value that mean of *subset* is < or > than mean of *fullset*.

    *subset* is a list of numbers.

    *fullset* is a list of numbers with *len(fullset) > len(subset)*

    *nrandom* is the number of random draws to use to compute the P-value.

    *withreplacement* should be *True* or *False*. If *True*, the random
    draws are done with replacement (same value can be drawn multiple
    times). If *False*, the draws are done without replacement.

    Computes the mean of the numbers in *subset*. Then performs *nrandom*
    draws of *len(subset)* samples (with or without replacement
    depending on the value of *withreplacement*) of *fullset*.
    Determines if the mean of the random subsets is < or >= to the mean
    in *subset*. If it is <, computes the fraction of random subsets
    where the random subsets have a mean >= *subset*. If it is >=, computes
    the fraction where the random subsets have a mean <= *subset*. Then returns
    the 2-tuple *(gt_or_lt, fraction)* where *gt_or_lt* is either "<" (if
    the mean of *subset* is >= to the random or ">". So *fraction*
    represents the one-sided P-value for the hypothesis that *subset*
    has a mean > or < than the value of a random subset from *fullset*.
    """
    nsubset = len(subset)
    if len(fullset) <= nsubset:
        raise ValueError("Length of fullset does not exceed that of subset")
    subset_total = sum(subset)
    nge = nle = 0
    for irandom in range(int(nrandom)):
        if withreplacement:
            isubset = []
            for idraw in range(nsubset):
                isubset.append(random.choice(fullset))
        else:
            isubset = random.sample(fullset, nsubset)
        totalrandom = sum(isubset)
        if totalrandom >= subset_total:
            nge += 1
        if totalrandom <= subset_total:
            nle += 1
    if nge >= 0.5 * nrandom:
        return ("<", nle / float(nrandom))
    else:
        return (">", nge / float(nrandom))


def PlotDistributionComparison(fullset, subset, fullsetname, subsetname,\
        plotfile, xlabel, ylabel, title, nrandom, withreplacement):
    """Compares two distributions and tests if one has a greater mean.

    This function can be generally used to compare and plot two
    distributions. Specifically, this function creates a plot of the
    distributions of integers in the two distributions *fullset* and
    *subset*. For generating this plot, there is no actual requirement
    that *subset* be a true subset of *fullset*.
    
    However, if *subset* is a true subset of *fullset*, then this
    function can also calculate and display the P-value
    for the hypothesis that the mean of *subset* is greater
    than the mean of *fullset*.

    This function uses ``pylab`` / ``matplotlib``. It will raise an Exception if
    these modules cannot be imported (if *PylabAvailable() == False*).

    The calling variables use LaTex format for strings. So for
    example, '$10^5$' will print the LaTex equivalent of this 
    string. Similarly, certain raw text strings (such as those
    including underscores) will cause problems if you do not
    escape the LaTex format meaning. For instance, 'x_label'
    will cause a problem since underscore is not valid
    outside of math mode in LaTex, so you would need to use
    'x\_label' to escape the underscore.

    CALLING VARIABLES:

    * *fullset* is a list of integers giving the first data set.

    * *subset* is a list of integers giving the second data set.
      If you are using *nrandom* then *subset* should be a true
      subset of *fullset* (for the calculated P-value to make
      sense) and in this case there is a strict requirement that
      *len(subset) < len(fullset)*.

    * *fullsetname* is a string giving the name used to label
      the distribution in *fullset*.

    * *subsetname* is a string giving the name used to label
      the distribution in *subset*.

    * *plotfile* is a string giving the name of the PDF plot file
      generated by this function. It must end in the extension ``.pdf``.
      If this file already exists, it is overwritten.

    * *xlabel* is a string giving the label placed on the x-axis.

    * *ylabel* is a string giving the label placed on the y-axis.

    * *title* is a string that is used to label the plot. If it is set
      to an expression that evaluates to *False*, then no title is displayed.

    * *nrandom* specifies how we calculate the P-value that the mean of 
      *subset* is < or >= the mean of *len(subset)* random samples drawn
      from *fullset*. If *nrandom* evaluates to *False*,
      then no P-value is computed or displayed. Otherwise, *nrandom*
      should give the name of random subsets of *fullset* that we test
      to compute the P-value. For example, a reasonable number might 
      be *nrandom=1e5*. Whether the draws are done with or without 
      replacement is specified by *withreplacement*.

    * *withreplacement* specifies how we calculate the P-value. The value
      of *withreplacement* is arbitrary if *nrandom* evaluates to *False*.
      Otherwise, *withreplacement* must be a *bool* variable of either *True*
      or *False*. If it is *True*, then the draws of the random subsets are 
      done with replacement (so the same number can be drawn multiple times).
      If it is *False*, then the draws are done without replacement (so
      the same number is drawn at most once).
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if not os.path.splitext(plotfile)[1].upper() == '.PDF':
        raise ValueError("plotfile does not end in PDF extension: %s " % plotfile)
    if not (isinstance(fullset, list) and isinstance(subset, list) and fullset and subset):
        raise ValueError("subset and fullset are not both non-empty lists.")
    if (len(fullset) <= len(subset)) and nrandom:
        raise ValueError("Cannot use nrandom unless subset has fewer entries than fullset")        
    for x in fullset + subset:
        if not isinstance(x, int):
            raise ValueError("fullset or subset contains non-integer entries")
    if not ((not nrandom) or (isinstance(nrandom, int) and nrandom > 0)):
        raise ValueError("Invalid value of nrandom: %s" % str(nrandom))
    if nrandom:
        if not (isinstance(withreplacement, bool)):
            raise ValueError("withreplacement must be of type bool")
    (lmargin, rmargin, tmargin, bmargin) = (0.16, 0.02, 0.03, 0.19)
    (figwidth, figheight) = (3, 1.95)
    if title:
        tmargin = 0.1
        figheight *= 1.1
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=9)
    matplotlib.rc('legend', fontsize=10)
    figure = pylab.figure(figsize=(figwidth, figheight), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    xmin = min(subset + fullset)
    xmax = max(subset + fullset)
    xs = [x for x in range(xmin, xmax + 1)]
    y_fullset = []
    y_subset = []
    for x in xs:
        y_fullset.append(fullset.count(x) / float(len(fullset)))
        y_subset.append(subset.count(x) / float(len(subset)))
    ymin = 0
    ymax = max(y_fullset + y_subset)
    pylab.plot(xs, y_fullset, 'bo-', label=fullsetname)
    pylab.plot(xs, y_subset, 'rs--', label=subsetname)
    xrange = xmax - xmin
    yrange = ymax - ymin
    ax.set_ylim(ymin - 0.02 * yrange, ymax + 0.09 * ymax)
    ax.set_xlim(xmin - 0.02 * xrange, xmax + 0.02 * xmax)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    pylab.legend(ncol=2, loc='upper right', numpoints=1, handlelength=1.5, borderaxespad=0)
    if nrandom:
        (gt_or_lt, p) = SubsetPValue(subset, fullset, nrandom, withreplacement)
        if p == 0:
            p = 1.0 / float(nrandom) # lower bound
        p = Base10Formatter(p, 3, 0, 2)
        pylab.text(0.99 * xmax, 0.78 * ymax, '%s $%s$ %s: $P = %s$' % (subsetname, gt_or_lt, fullsetname, p), fontsize=10, verticalalignment='bottom', horizontalalignment='right')
    if title:
        pylab.title(title, size=10)
    pylab.savefig(plotfile)
    time.sleep(0.5)



def CorrelationPlot(xs, ys, plotfile, xlabel, ylabel, corr=None, title=False):
    """Plots the correlation between two variables as a scatter plot.
    
    The data is plotted as a scatter plot.

    This function uses ``pylab`` / ``matplotlib``. It will raise an Exception if
    these modules cannot be imported (if *PylabAvailable() == False*).

    The calling variables use LaTex format for strings. So for
    example, '$10^5$' will print the LaTex equivalent of this 
    string. Similarly, certain raw text strings (such as those
    including underscores) will cause problems if you do not
    escape the LaTex format meaning. For instance, 'x_label'
    will cause a problem since underscore is not valid
    outside of math mode in LaTex, so you would need to use
    'x\_label' to escape the underscore.

    CALLING VARIABLES:

    * *xs* and *ys* are lists of numbers, with the lists
      being of the same length. Entry *xs[i]* is plotted
      on the x-axis agains entrie *ys[i]* on the y-axis.

    * *plotfile* is a string giving the name of the plot PDF file
      that we create. It should end in the extension ``.pdf``.
      If this plot already exists, it is overwritten.

    * *xlabel* is a string giving the label placed on the x-axis.

    * *ylabel* is a string giving the label placed on the y-axis.

    * *corr* specifies if we calculate and include a correlation
      coefficient on the plot. If it is *None*, then no
      correlation is computed. Otherwise, the coefficient
      is calculated using ``scipy`` (so this requires
      *ScipyAvailable() == True*). In this case, *corr* should
      be set to the string *Pearson* (to calculate the 
      Pearson linear correlation coefficient) or to the string
      *Spearman* (to calculate Spearman's rho rank-order correlation).
      In both cases, the correlations are reported along with the
      two-tailed P-values. They are written on the plot.

    * *title* is a string giving the title placed above the plot. 
      It can be *False* if no title is to be used. Otherwise, it should
      be the title string (using LaTex formatting, spaces are allowed).
      Is *False* by default.

    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if corr and not _scipyavailable:
        raise ImportError("Cannot use corr option since scipy is not available")
    if not os.path.splitext(plotfile)[1].upper() == '.PDF':
        raise ValueError("plotfile does not end in PDF extension: %s " % plotfile)
    if not (len(xs) == len(ys) >= 2):
        raise ValueError("xs and ys do not specify lists of the same length with >= 2 entries")
    (bigmargin, smallmargin) = (0.15, 0.03)
    (lmargin, rmargin, bmargin, tmargin) = (bigmargin, smallmargin, bigmargin, smallmargin)
    titlemargin = 0.07
    plotmargin = 0.05 # add this much above and below the last data point
    xsize = 2.5
    if title:
        tmargin += titlemargin
    ysize = xsize * (1.0 - lmargin - rmargin) / (1.0 - tmargin - bmargin)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=9)
    matplotlib.rc('legend', fontsize=10)
    figure = pylab.figure(figsize=(xsize, ysize), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.plot(xs, ys, 'b.')
    (xmin, xmax, ymin, ymax) = (min(xs), max(xs), min(ys), max(ys))
    xmargin = plotmargin * (xmax - xmin)
    ymargin = plotmargin * (ymax - ymin)
    ax.set_xlim([xmin - xmargin, xmax + xmargin])
    ax.set_ylim([ymin - ymargin, ymax + ymargin])
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    if title:
        pylab.title(title, size=10)
    if corr:
        if corr == 'Pearson':
            (r, p) = scipy.stats.pearsonr(xs, ys)
            r = '$R = %.2f' % r
        elif corr == 'Spearman':
            (r, p) = scipy.stats.spearmanr(xs, ys)
            r = '$\\rho = %.2f' % r
        else:
            raise ValueError("Invalid value of %s for corr" % corr)
        p = Base10Formatter(p, 2, 1, 2)
        text = '%s \,\, (P = %s)$' % (r, p)
        pylab.text(0.05, 0.95, text, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, size=10)
    pylab.savefig(plotfile)
    time.sleep(0.5) # this delay seems to help for some reason


def PlotLinearDensity(datalist, plotfile, xlabel, ylabel, title=False, fixymax=False):
    """Plots linear density of variable as a function of primary sequence.

    This function is designed to plot some variable (such as the number
    of epitopes as a function of the primary sequence position). It
    creates an output PDF plot *plotfile*. 

    The data is plotted as lines. If there is more than one
    data series to be plotted, a legend is included. 

    This function uses pylab / matplotlib. It will raise an Exception if
    these modules cannot be imported (if *PylabAvailable() == False*).

    The calling variables use LaTex format for strings. So for
    example, '$10^5$' will print the LaTex equivalent of this 
    string. Similarly, certain raw text strings (such as those
    including underscores) will cause problems if you do not
    escape the LaTex format meaning. For instance, 'x_label'
    will cause a problem since underscore is not valid
    outside of math mode in LaTex, so you would need to use
    'x\_label' to escape the underscore.

    CALLING VARIABLES:

    * *datalist*  is a list specifying the data to plot. It should
      be a list of one or more 2-tuples of the form 
      *(label, data)* where *label* is a string label used in
      the legend, and data is a list of 2-tuples *(x, y)*
      specifying the points to be plotted. 

    * *plotfile* is a string giving the name of the plot PDF file
      that we create. It should end in the extension ``.pdf``.
      If this plot already exists, it is overwritten.

    * *xlabel* is a string giving the label placed on the x-axis.

    * *ylabel* is a string giving the label placed on the y-axis.

    * *title* is a string giving the title placed above the plot. 
      It can be *False* if no title is to be used. Otherwise, it should
      be the title string (using LaTex formatting, spaces are allowed).
      Is *False* by default.

    * *fixymax* means that we fix the y-maximum to the specified value.
      This may be useful if you are making multiple plots for comparisons
      between them, and want them all to have the same y-maximum.
      Note that the value specified here is taken to be the data
      maximum -- the actually maximum of the y-axis is somewhat
      higher to provide some padding space. Is *False* by default.
    
    """
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if not os.path.splitext(plotfile)[1].upper() == '.PDF':
        raise ValueError("plotfile does not end in PDF extension: %s " % plotfile)
    if not (isinstance(datalist, list) and len(datalist) >= 1):
        raise ValueError("Invalid datalist")
    linestyles = ['b-', 'r-', 'g--', 'm-', 'c:', 'y--']
    linewidths = [2.5, 1, 1, 1, 1, 1]
    if len(datalist) > len(linestyles):
        raise ValueError("Number of data entries exceeds linestyles available for plotting")
    (lmargin, rmargin, tmargin, bmargin) = (0.08, 0.02, 0.03, 0.17)
    if title:
        tmargin = 0.11
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=9)
    matplotlib.rc('legend', fontsize=10)
    figure = pylab.figure(figsize=(5, 2), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    xmin = xmax = ymin = ymax = None
    for ((label, data), style, lw) in zip(datalist, linestyles, linewidths):
        xs = [tup[0] for tup in data]
        if xmin == None:
            (xmin, xmax) = (min(xs), max(xs))
        else:
            (xmin, xmax) = (min(xmin, min(xs)), max(xmax, max(xs)))
        ys = [tup[1] for tup in data]
        if ymin == None:
            (ymin, ymax) = (min(ys), max(ys))
        else:
            (ymin, ymax) = (min(ymin, min(ys)), max(ymax, max(ys)))
        pylab.plot(xs, ys, style, label=label, lw=lw)
    if title:
        pylab.title(title, size=10)
    ax.set_xlim([xmin, xmax])
    if fixymax != False:
        ymax = fixymax
    if len(datalist) > 1:
        ax.set_ylim([0, 1.25 * ymax]) # higher value to allow space for legend
        pylab.legend(loc='upper center', ncol=3, borderaxespad=0)
    else:
        ax.set_ylim([0, 1.1 * ymax]) # higher value to allow space for legend
    pylab.savefig(plotfile)
    time.sleep(0.5) # this delay seems to help for some reason



def CumulativeFractionPlot(datalist, plotfile, title, xlabel):
    """Creates a cumulative fraction plot.

    Takes a list of numeric data. Plots a cumulative fraction
    plot giving the fraction of the data points that are <=
    the indicated value.

    *datalist* is a list of numbers giving the data for which we
    are computing the cumulative fraction plot. Raises an 
    exception if this is an empty list.

    *plotfile* is the name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be '.pdf'.

    *title* is a string placed above the plot as a title. Uses LaTex
    formatting.

    *xlabel* is the label given to the X-axis. Uses LaTex formatting.

    This function uses pylab / matplotlib. It will raise an Exception if
    these modules cannot be imported (if PylabAvailable() is False).
    """
    if len(datalist) < 1:
        raise ValueError("datalist is empty")
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    datalist.sort() # sort from smallest to largest
    (xmin, xmax) = (datalist[0], datalist[-1])
    n = len(datalist)
    cumfracs = []
    cf = 0.0
    for x in datalist:
        cf += 1. / n
        cumfracs.append(cf)
    assert len(datalist) == len(cumfracs)
    assert abs(1.0 - cf) < 1e-7
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=12)
    fig = pylab.figure(figsize=(6, 4))
    (lmargin, rmargin, bmargin, tmargin) = (0.1, 0.01, 0.15, 0.1)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    pylab.plot(datalist, cumfracs, 'r-')
    pylab.gca().set_ylim([0, 1])
    pylab.gca().set_xlim([xmin, xmax])
    pylab.ylabel('cumulative fraction')
    pylab.xlabel(xlabel)
    pylab.title(title)
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def Base10Formatter(number, exp_cutoff, exp_decimal_digits, decimal_digits):
    """Converts a number into Latex formatting with scientific notation.

    Takes a number and converts it to a string that can be shown
    in LaTex using math mode. It is converted to scientific notation
    if the criteria specified by exp_cutoff.

    *number* the number to be formatted, should be a float or integer.
    Currently only works for numbers >= 0

    *exp_cutoff* convert to scientific notation if abs(math.log10(number)) >= this.

    *exp_decimal_digits* show this many digits after the decimal if number
    is converted to scientific notation.

    *decimal_digits* show this many digits after the decimal if number
    is NOT converted to scientific notation.
    
    The returned value is the LaTex' string. If the number is zero, the
    returned string is simply '0'.

    >>> Base10Formatter(103, 3, 1, 1)
    '103.0'

    >>> Base10Formatter(103.0, 2, 1, 1)
    '1.0 \\\\times 10^{2}'

    >>> Base10Formatter(103.0, 2, 2, 1)
    '1.03 \\\\times 10^{2}'

    >>> Base10Formatter(2892.3, 3, 1, 1) 
    '2.9 \\\\times 10^{3}'

    >>> Base10Formatter(0.0, 3, 1, 1) 
    '0'

    >>> Base10Formatter(0.012, 2, 1, 1)
    '1.2 \\\\times 10^{-2}'
  
    >>> Base10Formatter(-0.1, 3, 1, 1)
    Traceback (most recent call last):
        ...
    ValueError: number must be >= 0
    """
    if number < 0:
        raise ValueError('number must be >= 0')
    if number == 0:
        return '0'
    exponent = int(math.log10(number))
    if math.log10(number) < exponent and number < 1:
        exponent -= 1
    if abs(exponent) >= exp_cutoff:
        x = number / (10.**exponent)
        formatstr = '%.' + '%d' % exp_decimal_digits + 'f \\times 10^{%d}'
        return formatstr % (x, exponent)
    else:
        formatstr = '%.' + '%d' % decimal_digits + 'f'
        return formatstr % number


def SplitLabel(label, splitlen, splitchar):
    """Splits a string with a return if it exceeds a certain length.

    *label* a string giving the label we might split.

    *splitlen* the maximum length of a label before we attempt to 
    split it.

    *splitchar* the character added when splitting a label.

    If len(*label*) > *splitlen*, we attempt to split the label in the
    middle by adding *splitchar*. The label is split as close to the
    middle as possible while splitting at a space.

    No splitting as label length less than *splitlen*

    >>> SplitLabel('WT virus 1', 10, '\\n')
    'WT virus 1'

    Splitting of this label

    >>> SplitLabel('WT plasmid 1', 10, '\\n')
    'WT\\nplasmid 1'

    Splitting of this label

    >>> SplitLabel('mutated WT plasmid 1', 10, '\\n')
    'mutated WT\\nplasmid 1'

    """
    if len(label) <= splitlen:
        return label
    else:
        j = 0
        imid = len(label) // 2
        index = None
        while 0 <= imid - j <= imid + j < len(label):
            if label[imid - j].isspace():
                return "%s%s%s" % (label[ : imid - j], splitchar, label[imid - j + 1 : ])
            elif label[imid + j].isspace():
                return "%s%s%s" % (label[ : imid + j], splitchar, label[imid + j + 1 : ])
            j += 1
        else:
            return label # no white space to split


if __name__ == '__main__':
    import doctest
    doctest.testmod()
