from __future__ import absolute_import

import random

import dimod
import dwave_networkx as dnx
import numpy as np

from orang.conversions import bias_dtype, index_dtype
from orang.sample import sample_coo
from orang.solve import solve_coo, samples_dtype, energies_dtype

__all__ = ['OrangSolver', 'OrangSampler']


class OrangSolver(dimod.Sampler):
    """Tree decomposition-based solver for binary quadratic models.

    The Orang solver uses `tree decomposition`_ to find ground states of the
    given binary quadratic model.

    Examples:

        Create a solver

        >>> solver = orang.OrangSolver()

        Create a simple Ising problem

        >>> h = {'a': .1, 'b': 0}
        >>> J = {('a', 'b') : -1}

        We can use Orang to find the ground state

        >>> sampleset = solver.sample_ising(h, J)
        >>> sampleset.first.sample
        {'a': -1, 'b': -1}

        We can also take multiple reads to find states of increasing energy

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

    Accepted kwargs:

        * `num_reads`
        * `elimination_order`

    See :meth:`.sample` for descriptions.

    """

    properties = {'max_treewidth': 25}
    """Information about the solver.

    Properties:

        * `max_treewidth`: 25. The maximum treewidth_ allowed by the solver.

    .. _treewidth: https://en.wikipedia.org/wiki/Treewidth

    """

    def __init__(self):
        # make these object properties rathar than class properties
        self.parameters = dict(OrangSolver.parameters)
        self.properties = dict(OrangSolver.properties)

    def sample(self, bqm, num_reads=1, elimination_order=None):
        """Find ground states of a binary quadratic model.

        Args:
            bqm (:obj:`dimod.BinaryQuadraticModel`):
                A binary quadratic model.

            num_reads (int, optional, default=1):
                The total number of samples to draw. The samples are drawn in
                order of energy so if `num_reads=1` only the ground state will
                be returned. If `num_reads=2` the ground state and the first
                excited state. If `num_reads >= len(bqm)**2` then samples
                are duplicated.

            elimination_order (list, optional):
                The variable elimination order. Should be a list of the
                variables in the binary quadratic model. If None is provided,
                the min-fill heuristic [#gd]_ is used to generate one.

        Returns:
            :obj:`dimod.SampleSet`

        Raises:
            ValueError:
                The treewidth_ of the given bqm and elimination order cannot
                exceed the value provided in :attr:`.properties`.

        .. _treewidth: https://en.wikipedia.org/wiki/Treewidth

        .. [#gd] Gogate & Dechter, "A Complete Anytime Algorithm for Treewidth",
           https://arxiv.org/abs/1207.4109

        """
        max_samples = min(num_reads, 2**len(bqm))

        if not bqm:
            samples = np.empty((num_reads, 0), dtype=samples_dtype)
            energies = np.zeros((num_reads), dtype=energies_dtype)
            return dimod.SampleSet.from_samples(samples, bqm.vartype, energy=energies)

        if elimination_order is None:
            tw, elimination_order = dnx.min_fill_heuristic(bqm.adj)
        else:
            tw = dnx.elimination_order_width(bqm.adj, elimination_order)

        # developer note: we start getting bad_alloc errors above tw 25, this
        # should be fixed in the future
        if tw > self.properties['max_treewidth']:
            msg = ("maximum treewidth of {} exceeded. To see bqm's treewidth:\n"
                   ">>> import dwave_networkx as dnx\n"
                   ">>> dnx.elimination_order(bqm.adj, {})".format(self.properties['max_treewidth'], elimination_order))
            raise ValueError(msg)

        linear, quadratic, offset = bqm.to_numpy_vectors(
            variable_order=elimination_order,
            dtype=bias_dtype, index_dtype=index_dtype)  # to avoid copies

        # we used the elimination order as the variable order from the bqm, so
        # we use the default elimination order
        samples, energies = solve_coo(linear, quadratic, offset, bqm.vartype,
                                      max_complexity=tw+1,
                                      max_solutions=max_samples)

        # if we asked for more than the total number of distinct samples, we
        # just resample again starting from the beginning
        num_occurrences = np.ones(max_samples, dtype=np.intc)
        if num_reads > max_samples:
            q, r = divmod(num_reads, max_samples)
            num_occurrences *= q
            num_occurrences[:r] += 1

        return dimod.SampleSet.from_samples((samples, elimination_order), bqm.vartype,
                                            energy=energies + bqm.offset,
                                            num_occurrences=num_occurrences)


class OrangSampler(dimod.Sampler):
    """Tree decomposition-based sampler for binary quadratic models.

    The orang sampler uses tree-decomposition to sample from a boltzmann
    distribution defined by the given binary quadratic model.

    Examples:
        >>> sampler = orang.OrangSampler()
        >>> samples = sampler.sample_ising({}, {(0, 1): -1})

    """
    parameters = None
    properties = None

    def __init__(self):
        self.parameters = {'num_reads': [],
                           'elimination_order': [],
                           'beta': [],
                           'marginals': [],
                           'seed': []}
        self.properties = {'max_treewidth': 25}

    def sample(self, bqm, num_reads=1, elimination_order=None,
               beta=3,
               marginals=False,
               seed=None):
        """Draw samples and compute marginals of a binary quadratic model.

        Args:
            bqm (:obj:`dimod.BinaryQuadraticModel`):
                A binary quadratic model.

            num_reads (int, optional, default=1):
                The maximum number of solutions to return.

            elimination_order (list, optional):
                The variable elimination order. If None is provided, the min-fill heuristic is
                used to generate one.

            beta (float, optional, default=3.0):
                Boltzmann distribution inverse temperature parameter.

            marginals (bool, optional, default=False):
                Whether or not to compute the marginals. If True, they will be included in the
                return :obj:`~dimod.SampleSet`'s `info` field.

            seed (int, optional):
                random number generator seed.  Negative values will cause a time-based seed to be
                used.

        Returns:
            :obj:`dimod.SampleSet`

        Raises:
            ValueError: The treewidth_ of the given bqm and elimination order cannot exceed 25.

        .. _treewidth: https://en.wikipedia.org/wiki/Treewidth

        """

        if not bqm:
            info = dict(log_partition_function=0.0)
            if marginals:
                info['variable_marginals'] = {}
                info['interaction_marginals'] = {}
            samples = np.empty((num_reads, 0), dtype=samples_dtype)
            energies = bqm.energies(samples, dtype=energies_dtype)
            return dimod.SampleSet.from_samples(samples, bqm.vartype,
                                                energy=energies, info=info)

        if elimination_order is None:
            tw, elimination_order = dnx.min_fill_heuristic(bqm.adj)
        else:
            # this also checks the order against the bqm
            tw = dnx.elimination_order_width(bqm.adj, elimination_order)

        # developer note: we start getting bad_alloc errors above tw 25, this
        # should be fixed in the future
        if tw > self.properties['max_treewidth']:
            msg = ("maximum treewidth of {} exceeded. To see bqm's treewidth:\n"
                   ">>> import dwave_networkx as dnx\n"
                   ">>> dnx.elimination_order(bqm.adj, {})".format(self.properties['max_treewidth'], elimination_order))
            raise ValueError(msg)

        linear, quadratic, offset = bqm.to_numpy_vectors(
            variable_order=elimination_order,
            dtype=np.double, index_dtype=np.uintc)  # to avoid copies

        samples, data = sample_coo(linear, quadratic, offset, bqm.vartype,
                                   beta,
                                   max_complexity=tw+1,
                                   num_reads=num_reads,
                                   marginals=marginals,
                                   seed=seed)

        info = dict(log_partition_function=data['log_partition_function'])
        if marginals:
            info['variable_marginals'] = dict(zip(elimination_order,
                                                  data['variable_marginals']))

            low = -1 if bqm.vartype is dimod.SPIN else 0

            info['interaction_marginals'] = interaction_marginals = {}
            configs = (low, low), (1, low), (low, 1), (1, 1)
            for (i, j), probs in zip(data['interactions'],
                                     data['interaction_marginals']):
                u = elimination_order[i]
                v = elimination_order[j]
                interaction_marginals[(u, v)] = dict(zip(configs, probs))

        energies = bqm.energies((samples, elimination_order),
                                dtype=energies_dtype)

        return dimod.SampleSet.from_samples((samples, elimination_order),
                                            bqm.vartype,
                                            energy=energies,
                                            info=info)
