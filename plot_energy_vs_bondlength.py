import optimize_polymer as op
import pylab
import numpy as np

pa = {'t0': 2.5, 'alpha': 1.66*2.5, 'K': 21, 'boundary': 'periodic', 'chainlength': 6}

def bondchain(twobonds, pa):
    if pa['boundary'] == 'periodic':
        return (twobonds[0],twobonds[1])*(pa['chainlength']//2) + \
               (twobonds[0],)*(pa['chainlength']%2)
    elif pa['boundary'] == 'open':
        return (twobonds[0],twobonds[1])*((pa['chainlength'] - 1)//2) + \
               (twobonds[0],)*((pa['chainlength'] - 1)%2) + (0,)

x = y = np.arange(-0.4, 0, 0.005)
z = np.zeros(shape=(len(x), len(y)))
for i in range(len(x)):
    for j in range(len(y)):
        z[i,j] = op.totalEnergy(bondchain((x[i], y[j]), pa), pa)

pylab.imshow(z, origin='lower', extent=[-0.4,0,-0.4,0])
pylab.xlabel(ur"bond 1 length change (\u00c5)")
pylab.ylabel(ur"bond 2 length change (\u00c5)")
pylab.title('Benzene: energy (eV) vs bondlengths')
pylab.colorbar()

pylab.savefig('energy_vs_bondlength')
pylab.show()
