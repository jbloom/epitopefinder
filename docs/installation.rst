===================================
Installation
===================================

`epitopefinder`_ is written in `Python`_. It requires no external `Python`_ packages outside of the standard `Python`_ library, although it can optionally do some plotting that requires installation of `matplotlib`_. A few of the plotting functionalities also require the installation of `scipy`_.

`epitopefinder`_ has been tested on Mac OS X 10.6.8 with Python 2.6.7, `matplotlib`_ 1.0.0, and `scipy`_ 0.11. However, it should work with other similar versions as well and on other platforms as well.

Alignments are performed using `MUSCLE`_, which must be downloaded and installed separately. Scripts that utilize `MUSCLE`_ will require you to specify the location of the `MUSCLE`_ executable in the script input file.

To install `epitopefinder`_, first download the source ZIP repository `on GitHub`_. After unzipping the file, run the following commands::

    cd epitopefinder
    python setup.py build
    python setup.py install

The last command might require you to use::

    sudo python setup.py install
    
if you do not have privileges to the default installation directory, you might want to install to your local installation directory using::
    
    python setup.py install --user

These commands install the Python modules and and scripts. The scripts provide the most convenient high-level interface into the `epitopefinder`_ package. These scripts should be installed to be executable; they are present in the ``scripts`` subdirectory of the main package directory contained in the ZIP repository. The source code for the package is in the ``src`` subdirectory. There is also an ``examples`` subdirectory that contains examples. This documentation is found in the ``docs`` directory.

.. include:: weblinks.txt
