# Copyright 2019 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import List, Optional

import dimod
import numpy as np

from dimod.typing import Variable

from dwave.samplers.tree.sample import sample_bqm_wrapper
from dwave.samplers.tree.solve import solve_bqm_wrapper, samples_dtype, energies_dtype
from dwave.samplers.tree.utilities import elimination_order_width, min_fill_heuristic

__all__ = ['TreeDecompositionSolver', 'TreeDecompositionSampler']


class TreeDecompositionSolver(dimod.Sampler):
    """Tree decomposition-based solver for binary quadratic models.

    The tree decomposition solver uses `tree decomposition`_ to find ground
    states of the given binary quadratic model.

    .. Important:: The performance of this sampler is highly dependant on the
        treewidth of the problem.

    Examples:
        Create a solver:

        >>> from dwave.samplers import TreeDecompositionSolver
        >>> solver = TreeDecompositionSolver()

        Create a simple Ising problem:

        >>> h = {'a': .1, 'b': 0}
        >>> J = {('a', 'b') : -1}

        Find the ground state.

        >>> sampleset = solver.sample_ising(h, J)
        >>> sampleset.first.sample
        {'a': -1, 'b': -1}

        Take multiple reads to find states of increasing energy.

        >>> sampleset = solver.sample_ising(h, J, num_reads=3)
        >>> print(sampleset)
           a  b energy num_oc.
        0 -1 -1   -1.1       1
        1 +1 +1   -0.9       1
        2 -1 +1    0.9       1
        ['SPIN', 3 rows, 3 samples, 2 variables]

    .. _tree decomposition: https://en.wikipedia.org/wiki/Tree_decomposition

    """
    parameters = {'num_reads': [],
                  'elimination_order': ['max_treewidth']}
    """Keyword arguments accepted by the sampling methods.

    Examples:

        >>> from dwave.samplers import TreeDecompositionSolver
        >>> solver = TreeDecompositionSolver()
        >>> solver.parameters.keys()
        dict_keys(['num_reads', 'elimination_order'])

    See :meth:`.sample` for descriptions.

    """

    properties = {'max_treewidth': 25}
    """Information about the solver.

    Examples:

        >>> from dwave.samplers import TreeDecompositionSolver
        >>> solver = TreeDecompositionSolver()
        >>> solver.properties
        {'max_treewidth': 25}

    """

    def __init__(self):
        # make these object properties rather than class properties
        self.parameters = dict(self.parameters)
        self.properties = dict(self.properties)

    def sample(self, bqm: dimod.BinaryQuadraticModel, num_reads: Optional[int] = 1,
               elimination_order: Optional[List[Variable]] = None) -> dimod.SampleSet:
        """Find ground states of a binary quadratic model.

        Args:
            bqm:
                Binary quadratic model (BQM).

            num_reads:
                Total number of samples to draw. The samples are drawn in
                order of energy so if ``num_reads=1``, only the ground state is
                returned. If ``num_reads=2``, the ground state and the first
                excited state are returned. If ``num_reads >= len(bqm)**2``,
                samples are duplicated.

            elimination_order:
                Variable elimination order. Should be a list of the
                variables in the binary quadratic model. If None, the min-fill
                heuristic [#gd]_ is used to generate one.

        Raises:
            ValueError:
                The treewidth_ of the given BQM and elimination order cannot
                exceed the value provided in :attr:`.properties`.

        .. _treewidth: https://en.wikipedia.org/wiki/Treewidth

        .. [#gd] Gogate & Dechter, "A Complete Anytime Algorithm for Treewidth",
           https://arxiv.org/abs/1207.4109

        """
        if not bqm:
            samples = np.empty((num_reads, 0), dtype=samples_dtype)
            energies = bqm.energies(samples, dtype=energies_dtype)
            return dimod.SampleSet.from_samples(samples, bqm.vartype, energy=energies)

        bqm = dimod.as_bqm(bqm, copy=True, dtype=float)

        max_samples = min(num_reads, 2**len(bqm))

        if elimination_order is None:
            tree_width, elimination_order = min_fill_heuristic(bqm)
        else:
            tree_width = elimination_order_width(bqm, elimination_order)

        # developer note: we start getting bad_alloc errors above tree_width 25, this
        # should be fixed in the future
        if tree_width > self.properties['max_treewidth']:
            raise ValueError(
                f"maximum treewidth of {self.properties['max_treewidth']} exceeded. "
                "To see the bqm's treewidth:\n"
                ">>> import dwave_networkx as dnx\n"
                f">>> dnx.elimination_order_width(bqm.adj, {elimination_order})"
                )

        max_complexity = tree_width + 1

        # relabel bqm variables so that we only work with linear indices
        bqm_copy, int_to_var = bqm.relabel_variables_as_integers(inplace=False)
        var_to_int = {v: k for k, v in int_to_var.items()}

        # relabel variables in the elimination order as well
        elimination_order = [var_to_int.get(var, var) for var in elimination_order]

        samples, energies = solve_bqm_wrapper(bqm=bqm_copy,
                                              order=elimination_order,
                                              max_complexity=max_complexity,
                                              max_solutions=max_samples
                                              )

        # if we asked for more than the total number of distinct samples, we
        # just resample again starting from the beginning
        num_occurrences = np.ones(max_samples, dtype=np.intc)
        if num_reads > max_samples:
            q, r = divmod(num_reads, max_samples)
            num_occurrences *= q
            num_occurrences[:r] += 1

        return dimod.SampleSet.from_samples((samples, bqm.variables),
                                            bqm.vartype,
                                            energy=energies,
                                            num_occurrences=num_occurrences)


class TreeDecompositionSampler(dimod.Sampler):
    """Tree decomposition-based sampler for binary quadratic models.

    The tree decomposition sampler uses `tree decomposition`_ to sample from a
    `Boltzmann distribution`_ defined by the given binary quadratic model.

    .. Important:: The performance of this sampler is highly dependant on the
        treewidth of the problem.

    Examples:
        Create a sampler:

        >>> from dwave.samplers import TreeDecompositionSampler
        >>> sampler = TreeDecompositionSampler()

        Create a simple Ising problem:

        >>> h = {'a': .1, 'b': 0}
        >>> J = {('a', 'b') : -1}

        Sample from the given problem.

        >>> sampleset = sampler.sample_ising(h, J, num_reads=100,
        ...                                  elimination_order=['a', 'b'])
        >>> sampleset.first.sample
        {'a': -1, 'b': -1}

        View information about the distribution.

        >>> variable_marginals = sampleset.info['variable_marginals']
        >>> round(variable_marginals['a'], 3)  # prob(a == 1)
        0.354
        >>> round(1 - variable_marginals['b'], 3)  # prob(b == -1)
        0.645

        >>> pair_marg = sampleset.info['interaction_marginals']
        >>> round(pair_marg[('a', 'b')][(1, -1)], 3)  # prob(a == 1 & b == -1)
        0.001

    .. _tree decomposition: https://en.wikipedia.org/wiki/Tree_decomposition

    .. _Boltzmann distribution: https://en.wikipedia.org/wiki/Boltzmann_distribution

    """
    parameters = {'num_reads': [],
                  'elimination_order': ['max_treewidth'],
                  'beta': [],
                  'marginals': [],
                  'seed': []}
    """Keyword arguments accepted by the sampling methods.

    Examples:

        >>> from dwave.samplers import TreeDecompositionSampler
        >>> sampler = TreeDecompositionSampler()
        >>> sampler.parameters.keys()
        dict_keys(['num_reads', 'elimination_order', 'beta', 'marginals', 'seed'])

    See :meth:`.sample` for descriptions.

    """

    properties = {'max_treewidth': 25}
    """Information about the sampler.

    Examples:

        >>> from dwave.samplers import TreeDecompositionSampler
        >>> sampler = TreeDecompositionSampler()
        >>> sampler.properties
        {'max_treewidth': 25}

    """

    def __init__(self):
        # make these object properties rather than class properties
        self.parameters = dict(self.parameters)
        self.properties = dict(self.properties)

    def sample(self, bqm: dimod.BinaryQuadraticModel, num_reads: Optional[int] = 1,
               elimination_order: Optional[List[Variable]] = None, beta: Optional[float] = 3.0,
               marginals: Optional[bool] = True, seed: Optional[int] = None) -> dimod.SampleSet:
        """Draw samples and compute marginals of a binary quadratic model.

        Args:
            bqm:
                Binary quadratic model.

            num_reads:
                Number of samples to draw.

            elimination_order:
                Variable elimination order. Should be a list of the
                variables in the binary quadratic model. If None, the min-fill
                heuristic [#gd]_ is used to generate one.

            beta:
                `Boltzmann distribution`_ inverse temperature parameter.

            marginals:
                Whether or not to compute the marginals. If True, they are
                included in the return :obj:`~dimod.SampleSet`'s ``info`` field.
                See example in :class:`.TreeDecompositionSampler`.

            seed:
                Random number generator seed. Negative values cause a
                time-based seed to be used.

        Returns:
            Returned :attr:`dimod.SampleSet.info` contains:

                * ``'log_partition_function'``: The log partition function.


            If ``marginals=True``, also contains:

                * ``'variable_marginals'``: Dict of the form ``{v: p, ...}``, where
                  ``v`` is a variable in the binary quadratic model and
                  ``p = prob(v == 1)``.
                * ``'interaction_marginals'``: Dict of the form
                  ``{(u, v): {(s, t): p, ...}, ...}``, where ``(u, v)`` is an
                  interaction in the binary quadratic model and
                  ``p = prob(u == s & v == t)``.

        Raises:
            ValueError:
                The treewidth_ of the given bqm and elimination order cannot
                exceed the value provided in :attr:`.properties`.

        .. _treewidth: https://en.wikipedia.org/wiki/Treewidth
        .. _Boltzmann distribution: https://en.wikipedia.org/wiki/Boltzmann_distribution
        .. [#gd] Gogate & Dechter, "A Complete Anytime Algorithm for Treewidth",
           https://arxiv.org/abs/1207.4109

        """
        if not bqm:
            info = {'log_partition_function': 0.0}
            if marginals:
                info['variable_marginals'] = {}
                info['interaction_marginals'] = {}
            samples = np.empty((num_reads, 0), dtype=samples_dtype)
            energies = bqm.energies(samples, dtype=energies_dtype)
            return dimod.SampleSet.from_samples(samples, bqm.vartype,
                                                energy=energies, info=info)

        bqm = dimod.as_bqm(bqm, copy=True, dtype=float)

        if elimination_order is None:
            # note that this does not respect the given seed
            tree_width, elimination_order = min_fill_heuristic(bqm)
        else:
            # this also checks the order against the bqm
            tree_width = elimination_order_width(bqm, elimination_order)

        # developer note: we start getting bad_alloc errors above tree_width 25, this
        # should be fixed in the future
        if tree_width > self.properties['max_treewidth']:
            raise ValueError(
                f"maximum treewidth of {self.properties['max_treewidth']} exceeded. "
                "To see the bqm's treewidth:\n"
                ">>> import dwave_networkx as dnx\n"
                f">>> dnx.elimination_order_width(bqm.adj, {elimination_order})"
                )

        # relabel bqm variables so that we only work with linear indices
        bqm_copy, int_to_var = bqm.relabel_variables_as_integers(inplace=False)
        var_to_int = {v: k for k, v in int_to_var.items()}

        # relabel variables in the elimination order as well
        elimination_order = [var_to_int.get(var, var) for var in elimination_order]

        max_complexity = tree_width + 1

        samples, data = sample_bqm_wrapper(bqm=bqm_copy,
                                           beta=beta,
                                           max_complexity=max_complexity,
                                           order=elimination_order,
                                           marginals=marginals,
                                           num_reads=num_reads,
                                           seed=seed)

        info = {'log_partition_function': data['log_partition_function']}

        if marginals:
            info['variable_marginals'] = {}
            for i in elimination_order:
                info['variable_marginals'][int_to_var.get(i, i)] = data['variable_marginals'][i]

            info['interaction_marginals'] = {}
            low = -1 if bqm.vartype is dimod.SPIN else 0
            configs = (low, low), (1, low), (low, 1), (1, 1)
            for (i, j), probs in zip(data['interactions'],
                                     data['interaction_marginals']):
                u = int_to_var.get(i, i)
                v = int_to_var.get(j, j)
                info['interaction_marginals'][(u, v)] = dict(zip(configs, probs))

        energies = bqm.energies((samples, bqm.variables),
                                dtype=energies_dtype)

        return dimod.SampleSet.from_samples((samples, bqm.variables),
                                            bqm.vartype,
                                            energy=energies,
                                            info=info)
