# Copyright 2019 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
# =============================================================================
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
