import optimize_polymer as op
import numpy as np
import pylab

pa = {'t0': 2.5, 'K': 21, 'boundary': 'periodic', 'chainlength': 6}
alphas = np.arange(8, 11, 0.1)

bonds = {}
for alpha in alphas:
    pa['alpha'] = alpha
    bonds[alpha] = op.optimalDimer(pa)
dimerization = np.array([abs(bonds[a][0]-bonds[a][1]) for a in alphas])

pylab.plot(alphas, dimerization, 'o')
pylab.grid(True)
pylab.xlabel(ur'electron-phonon coupling (eV/\u00c5)')
pylab.ylabel(ur'bond length difference (\u00c5)')
pylab.title('Benzene: bond dimerization vs alpha')
pylab.savefig('benzene_dimerization')

pylab.show()
