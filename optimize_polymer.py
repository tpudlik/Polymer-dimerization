# optimize_polymer.py
#
# author: Ted Pudlik
# created: Jan 11th, 2013
# this version (rev 4): Jan 15th, 2013

"""
This module finds the energy-minimizing bond lengths of a SSH polyene, under
the assumption that the optimal bond lengths are either all equal, or take on
two alternating lengths (as in, e.g., cyclobutadiene).
"""

from scipy.optimize import fmin
from electronic_energy import electronicEnergy

def optimalDimer(pa, initialBonds=[0.01, -0.01], **kwargs):
    """Find the energy-minimizing bond lengths of a polyene.

    Input arguments:
        pa: a dictionary containing the values of the system parameters.  These
            parameters must include,
            t0: the "bare" hopping (at equilibrium bond length), in eV
            alpha: the electron-phonon coupling, in units of eV per Angstrom
            K: the spring constant of the sigma bonds, in units of eV per
                Angstrom squared
            chainlength: the number of carbon atoms in the polyene
            boundary: the type of boundary conditions ('periodic' or 'open')
        initialBonds: the starting guess for the bond lengths, relative to the
            undimerized equilibrium lengths
        keyword arguments: any options accepted by scipy.optimize.fmin.

    Returns: numpy array containing the lengths of the energy-minimizing
        short and long bonds.

    This function uses a simplex algorithm to find the configuration (both
    phononic and electronic)  which minimizes the total energy of the polyene.
    It is assumed that there are at most two bond lengths in the polymer, and
    these bond lengths alternate.  (Note that the degenerate case of two equal
    bond lengths is included.)"""
    
    if not set(['t0','alpha','chainlength','boundary','K']).issubset(set(pa.keys())):
        raise MissingParameters(set(pa.keys()))
    if pa['chainlength'] < 3: raise Exception
    def dimerEnergy(twobonds):
        if pa['boundary'] == 'open':
            bonds = (twobonds[0],twobonds[1])*((pa['chainlength'] - 1)//2) + (twobonds[0],)*((pa['chainlength'] - 1)%2) + (0,)
        elif pa['boundary'] == 'periodic':
            bonds = (twobonds[0],twobonds[1])*(pa['chainlength']//2) + (twobonds[0],)*(pa['chainlength']%2)
        return totalEnergy(bonds, pa)

    return fmin(dimerEnergy, initialBonds, **kwargs)

def totalEnergy(bonds, pa):
    """Compute the total ground state energy of a polyene, given the bond
    lengths and system parameters."""
    return electronicEnergy(bonds, pa) + pa['K']*sum(bond**2 for bond in bonds)/2    

class MissingParameters(Exception):
    def __init__(self, parameters):
        self.parameters = parameters
    def __str__(self):
        return 'Passed parameters ' + repr(self.parameters) + ', missing ' + \
               repr(set(['t0','alpha','chainlength','boundary','K']) - self.parameters)
