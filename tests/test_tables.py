import unittest

# import dimod
import numpy as np

from orang.tables import Tables, solve


class TestSolve(unittest.TestCase):
    def test_single_variable_binary_down(self):
        ldata = [1]
        irow = []
        icol = []
        qdata = []
        offset = 0
        vartype = 'BINARY'
        order = [0]
        complexity = 2
        num_reads = 1

        tables = Tables.from_coo(ldata, (irow, icol, qdata), offset, vartype)
        samples, energies = solve(tables, order, complexity, num_reads)

        self.assertEqual(samples.shape, (1, 1))
        self.assertEqual(energies.shape, (1,))

        np.testing.assert_array_equal(samples, [[0]])
        np.testing.assert_array_equal(energies, [0])

    def test_single_variable_binary_up(self):
        ldata = [-1]
        irow = []
        icol = []
        qdata = []
        offset = 0
        vartype = 'BINARY'
        order = [0]
        complexity = 2
        num_reads = 1

        tables = Tables.from_coo(ldata, (irow, icol, qdata), offset, vartype)
        samples, energies = solve(tables, order, complexity, num_reads)

        self.assertEqual(samples.shape, (1, 1))
        self.assertEqual(energies.shape, (1,))

        np.testing.assert_array_equal(samples, [[1]])
        np.testing.assert_array_equal(energies, [-1])
