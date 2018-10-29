from orang._orang import sample

import dimod

bqm = dimod.BinaryQuadraticModel.from_ising({}, {'ab': -1, 'bc': -1})

ss = sample(bqm, 100)

for sample, energy in ss.data(['sample', 'energy']):
    print(bqm.energy(sample), energy)
