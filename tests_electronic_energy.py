# tests_electronic_energy: tests for the electronic_energy module, rev 3+.
#
# Ted Pudlik
# created: Jan 14th, 2013
# this version (rev 4): Jan 15th, 2013

import electronic_energy as ds
from scipy import array
import numpy as np

# ground state energy from sorted list of state energies
assert ds.groundStateEnergy(array([-5, -3, 0, 1])) == -16
assert ds.groundStateEnergy(array([-7, -2, -1, 3, 4])) == -19
assert ds.groundStateEnergy(array([-3,3])) == -6

# Modified modulo function
assert ds.wrap(1, 3) == 1
assert ds.wrap(2, 3) == 2
assert ds.wrap(3, 3) == 3
assert ds.wrap(4, 3) == 1
assert ds.wrap(5, 3) == 2
assert ds.wrap(0, 3) == 3

# Inner product
assert ds.inner(1, 2) == 0
assert ds.inner(2, 2) == 1
assert ds.inner(0, 1) == 0
assert ds.inner(0, 0) == 0 # note that 0 is not the ket |0>, but just zero!

# hopping t
pa = {'t0': 5, 'alpha': 2, 'K': 10, 'chainlength': 5, 'boundary': 'open'}
bonds = (1,1,1,1,0)
assert ds.t(1, pa, bonds) == pa['t0'] - pa['alpha']
assert ds.t(5, pa, bonds) == 0
assert ds.t(4, pa, bonds) == pa['t0'] - pa['alpha']
pa = {'t0': 5, 'alpha': 2, 'K': 10, 'chainlength': 5, 'boundary': 'periodic'}
bonds = (1,1,2,1,1)
assert ds.t(5, pa, bonds) == pa['t0'] - pa['alpha']
assert ds.t(3, pa, bonds) == pa['t0'] - 2*pa['alpha']

# Hamiltonian matrix elements
pa = {'t0': 5, 'alpha': 2, 'K': 10, 'chainlength': 5, 'boundary': 'open'}
bonds = (1,1,1,1,0)
assert ds.hamiltonianElement(1, 3, bonds, pa) == 0
assert ds.hamiltonianElement(5, 1, bonds, pa) == 0
assert ds.hamiltonianElement(1, 2, bonds, pa) == -(pa['t0'] - pa['alpha'])
assert ds.hamiltonianElement(2, 1, bonds, pa) == -(pa['t0'] - pa['alpha'])
assert ds.hamiltonianElement(3, 4, bonds, pa) == -(pa['t0'] - pa['alpha'])
pa = {'t0': 5, 'alpha': 2, 'K': 10, 'chainlength': 5, 'boundary': 'periodic'}
bonds = (1,2,1,4,1)
assert ds.hamiltonianElement(2, 3, bonds, pa) == -(pa['t0'] - 2*pa['alpha'])
assert ds.hamiltonianElement(4, 5, bonds, pa) == -(pa['t0'] - 4*pa['alpha'])
assert ds.hamiltonianElement(1, 5, bonds, pa) == -(pa['t0'] - pa['alpha'])
assert ds.hamiltonianElement(5, 1, bonds, pa) == -(pa['t0'] - pa['alpha'])

# Electronic ground state energy: full function tests
## Ethene
pa = {'t0': 5, 'alpha': 2, 'K': 10, 'chainlength': 2, 'boundary': 'open'}
assert round(ds.electronicEnergy((1,0), pa) + 2*(pa['t0']-pa['alpha']), 7) == 0
pa = {'t0': 8, 'alpha': 0.23, 'K': 7, 'chainlength': 2, 'boundary': 'open'}
assert round(ds.electronicEnergy((1,0), pa) + 2*(pa['t0']-pa['alpha']), 7) == 0
## Propene
pa = {'t0': 3, 'alpha': 1, 'K': 7, 'chainlength': 3, 'boundary': 'open'}
bonds = (1,1,0)
assert round(ds.electronicEnergy(bonds, pa) + 2*np.sqrt((pa['t0'] - pa['alpha']*bonds[0])**2 + (pa['t0'] - pa['alpha']*bonds[1])**2), 7) == 0
pa = {'t0': 3, 'alpha': 2, 'K': 7, 'chainlength': 3, 'boundary': 'open'}
bonds = (2,3,0)
assert round(ds.electronicEnergy(bonds, pa) + 2*np.sqrt((pa['t0'] - pa['alpha']*bonds[0])**2 + (pa['t0'] - pa['alpha']*bonds[1])**2), 7) == 0
pa = {'t0': 3, 'alpha': 2.1, 'K': 7, 'chainlength': 3, 'boundary': 'open'}
bonds = (1.2,1.3,0)
assert round(ds.electronicEnergy(bonds, pa) + 2*np.sqrt((pa['t0'] - pa['alpha']*bonds[0])**2 + (pa['t0'] - pa['alpha']*bonds[1])**2), 7) == 0
## Cyclopropene
pa = {'t0': 2, 'alpha': 1, 'K': 7, 'chainlength': 3, 'boundary': 'periodic'}
bonds = (1,1,1)
assert round(ds.electronicEnergy(bonds, pa) + 3, 7) == 0
pa = {'t0': 3, 'alpha': 1, 'K': 7, 'chainlength': 3, 'boundary': 'periodic'}
bonds = (1, 2, 1)
assert round(ds.electronicEnergy(bonds, pa) + np.sqrt(33), 7) == 0

print "tests pass"
