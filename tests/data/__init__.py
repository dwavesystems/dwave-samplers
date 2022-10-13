from os.path import abspath, dirname, join

import dimod


fp39 = join(dirname(abspath(__file__)), 'FrustTriangleL39.tsv')

bqm_L39 = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
with open(fp39, 'r') as fo:
    for line in fo:
        u, v, bias = line.split()
        bqm_L39.add_interaction(int(u), int(v), float(bias))


def pos_L39(v):
    y = v % 39
    x = (v - y) // 39
    return x, y
