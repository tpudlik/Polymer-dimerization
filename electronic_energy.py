# electronic_energy.py
#
# author: Ted Pudlik
# created: Jan 11th, 2013
# this version (rev 4): Jan 15th, 2013

"""
Module for computing the energy of a polyene's electronic ground state.

The module finds the ground state energies of both cyclic and open polyenes,
using a simple noninteracting electron tight-binding model.

The main function is electronicEnergy.
"""

from scipy.linalg import eigvalsh
import numpy as np

def electronicEnergy(bonds, pa):
    """Compute the electronic contribution to the ground state energy of a
    polyene.

    Input arguments:
        bonds: a tuple of length equal to the number of carbon atoms in the
            polyene, containing the lengths of the bonds measured relative
            to their equilibrium lengths (i.e., bond[i] = u_{i+2} - u_{i+1}).
        pa: a dictionary containing the values of the system parameters.  These
        parameters must include,
            t0: the "bare" hopping (at equilibrium bond length), in eV
            alpha: the electron-phonon coupling, in units of eV per Angstrom
            chainlength: the number of carbon atoms in the polyene
            boundary: the type of boundary conditions ('periodic' or 'open')

    Returns: the ground state electronic energy of the specified polyene.
            
    Note that the length of the bonds tuple should be equal to the number of
    atoms even in the 'open' boundary condition case, in which the number of
    bonds in the molecule is one less than the number of atoms.  In this last
    case, the last entry of the bonds tuple should be set to 0."""
    
    if not set(['t0','alpha','chainlength','boundary']).issubset(set(pa.keys())):
        raise MissingParameters(set(pa.keys()))
    if len(bonds) != pa['chainlength']:
        raise BondNumberNotEqualChainLength(len(bonds), pa['chainlength'])
    if pa['boundary'] is 'open' and bonds[-1] != 0:
        raise OpenBCImplyLastBondAbsent(bonds[-1])
    energies = eigvalsh(hamiltonianMatrix(bonds, pa))
    return groundStateEnergy(energies)

def hamiltonianMatrix(bonds, pa):
    """Construct the hamiltonian matrix.

    The input arguments are those of the parent function electronicEnergy"""
    
    basis = electronicBasis(pa['chainlength'])
    hamiltonianMatrix = np.zeros(shape=(len(basis),len(basis)))
    for i in range(len(basis)):
        for j in range(len(basis)):
            hamEm = hamiltonianElement(basis[i], basis[j], bonds, pa)
            if hamEm is not 0: hamiltonianMatrix[i, j] = hamEm
    return hamiltonianMatrix

def electronicBasis(chain):
    """Return the basis for the electronic states of a polyene of length chain.

    The basis states are tight-binding orbitals, one for each CH unit, each with
    the capacity to hold two electrons.  They are represented as integers
    corresponding to the CH unit number."""
    
    return tuple(range(1, chain+1))

def hamiltonianElement(bra, ket, bonds, pa):
    """Return the matrix element of the Hamiltonian.

    The matrix element is taken between the basis vectors bra and ket.  If the
    matrix element vanishes identically, return integer 0."""
    
    hoppings = [inner(bra, wrap(ket+1, pa['chainlength'])), inner(bra, wrap(ket-1, pa['chainlength']))]
    if hoppings == [0, 0]:
        return 0
    else:
        return -t(ket, pa, bonds)*hoppings[0] - t(wrap(ket-1, pa['chainlength']), pa, bonds)*hoppings[1]

def t(ket, pa, bonds):
    """Hopping between sites ket and ket+1, given the parameters.

    If the chain is not periodic, the hopping between the last and first sites
    is identically zero."""
    
    if ket == pa['chainlength'] and pa['boundary'] == 'open':
        return 0
    else:
        return pa['t0'] - pa['alpha']*bonds[ket-1]
    
def inner(bra, ket):
    """Inner product of basis states bra and ket"""
    return 0 if bra == 0 or ket == 0 else bra == ket

def wrap(site_index, chain):
    """A modified modulo function that maps integers to the range (1, chain).

    The usual modulo sign "%" maps integers to the range (0, chain - 1)."""
    
    if site_index == chain or site_index == 0:
        return chain
    else:
        return site_index % chain

def groundStateEnergy(energies):
    """Return the energy of the half-filled ground state.

    The list of electronic state energies is assumed to be sorted in ascending
    order."""
    
    chain = len(energies)
    if chain % 2 == 0:
        return 2*sum(energies[0:(chain/2)])
    else:
        return 2*sum(energies[0:(chain/2)]) + \
                     energies[chain/2]

class BondNumberNotEqualChainLength(Exception):
    def __init__(self, bonds, chainlength):
        self.bonds = bonds
        self.chainlength = chainlength
    def __str__(self):
        return 'bonds: ' + repr(self.bonds) + \
               ', chainlength: ' + repr(self.chainlength)

class OpenBCImplyLastBondAbsent(Exception):
    def __init__(self, lastBond):
        self.bond = lastBond
    def __str__(self):
        return 'Bond between first and last unit has length ' + repr(self.bond)

class MissingParameters(Exception):
    def __init__(self, parameters):
        self.parameters = parameters
    def __str__(self):
        return 'Passed parameters ' + repr(self.parameters) + ', missing ' + \
               repr(set(['t0','alpha','chainlength','boundary']) - self.parameters)
