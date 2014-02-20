Polyene Dimerization Simulations
--------------------------------
by Ted Pudlik, tpudlik@bu.edu

created: Jan 2013

This software package calculates the energy of polyenes in the SSH (Su,
Schriefer, Heeger) model and finds energy-minimizing dimer configurations.
Both open- and closed-chain polyenes are supported, and both electronic and
phononic energies are taken into account.  Reasonable performance can be
achieved even for polyenes tens or hundreds of units in length.


What do I need to run it?
-------------------------

The program is implemented in Python and depends on the SciPy and NumPy
libraries for optimization and linear algebra operations, and PyLab for
plotting.  If you don't have Python on your system yet, the easiest way to
acquire a distribution that will support this software is to download Enthought
Python, currently freely provided at http://www.enthought.com/products/epd.php .


What's in the box?
------------------

The following files make up the package:

    electronic_energy.py: a module for finding the electronic contribution to
        the energy of a polyene.  This module can be used by itself, but its
        primary purpose is to be invoked by optimize_polymer.py.
    optimize_polymer.py: the main module of this package.  It contains
        functions for finding and minimizing the total energy of a polyene.
    plot_benzene_vs_cyclohexatriene.py: creates a plot of bond dimerization in
        a six-carbon ring polyene as a function of the electron-phonon coupling.
    plot_dimerization_vs_ringsize.py: creates a plot of bond dimerization as a
        function of ring size for a range of even-sized polyenes.
    plot_energy_vs_bondlength.py: creates a heat map of energy versus bond
        lengths for benzene.
    tests_electronic_energy.py: unit tests of components of the
        electronic_energy.py module.

