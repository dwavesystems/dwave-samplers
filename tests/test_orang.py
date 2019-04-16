import unittest

import dimod
import numpy as np

from orang.solve import solve_coo


class TestSolveCoo(unittest.TestCase):
    def test_single_interaction_spin(self):
        ldata = [.1, 0]
        irow, icol = [0], [1]
        qdata = [-1]
        offset = 0
        vartype = dimod.SPIN
        kwargs = dict()
        complexity = 5

        quadratic = (irow, icol, qdata)
        samples, energies = solve_coo(ldata, quadratic, offset, vartype, complexity, **kwargs)

        np.testing.assert_array_equal(samples, [[-1, -1]])
        np.testing.assert_array_equal(energies, [-1.1])

    def test_empty(self):
        ldata = []
        irow, icol = [], []
        qdata = []
        offset = 0
        vartype = dimod.SPIN
        kwargs = dict()
        complexity = 5

        quadratic = (irow, icol, qdata)
        with self.assertRaises(ValueError):
            solve_coo(ldata, quadratic, offset, vartype, complexity, **kwargs)
