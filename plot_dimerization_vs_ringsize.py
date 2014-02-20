import optimize_polymer as op
import numpy as np
import pylab

def get_dimerization(ring_sizes, pa):
    bonds = {}
    for size in ring_sizes:
        pa['chainlength'] = size
        bonds[size] = op.optimalDimer(pa)
    return np.array([abs(bonds[i][0]-bonds[i][1]) for i in ring_sizes])

pa = {'t0': 2.5, 'alpha': 1.66*2.5, 'K': 21}
ring_sizes = range(4, 60, 2)
odd_ring_sizes = range(3, 60, 2)
markersize = 40

pa['boundary'] = 'periodic'
dimerizationPeriodic = get_dimerization(ring_sizes, pa)
dimerizationPeriodicOdd = get_dimerization(odd_ring_sizes, pa)
pa['boundary'] = 'open'
dimerizationOpen = get_dimerization(ring_sizes, pa)
dimerizationOpenOdd = get_dimerization(odd_ring_sizes, pa)

fig = pylab.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(ring_sizes, dimerizationPeriodic, c='b', marker='o',
            label='rings: 2n CH units', s=markersize)
ax1.scatter(odd_ring_sizes, dimerizationPeriodicOdd, c='#99ccff', marker='o',
            label='rings: 2n+1 CH units', s=markersize)
ax1.scatter(ring_sizes, dimerizationOpen, c='r', marker='s',
            label='open chains: 2n CH units', s=markersize)
ax1.scatter(odd_ring_sizes, dimerizationOpenOdd, c='#ff9999', marker='s',
            label='open chains: 2n+1 CH units', s=markersize)
ax1.plot([ring_sizes[0]-2,ring_sizes[-1]+2],[0.084,0.084], c='#ffa500',
            label='polyacetylene')
pylab.legend()
pylab.grid(True)
pylab.xlim([ring_sizes[0]-2, ring_sizes[-1]+2])
pylab.xlabel('chain length')
pylab.ylabel(ur"bond length difference (\u00c5)")
pylab.title('t0 = ' + str(pa['t0']) + ' eV, alpha = ' + str(pa['alpha']) + \
            ur' eV/\u00c5, K = ' + str(pa['K']) + ur' eV/\u00c5^2')
pylab.xticks(ring_sizes[::2])

pylab.savefig('dimerization_vs_ringsize')
pylab.show()

