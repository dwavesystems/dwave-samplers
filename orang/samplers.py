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
    """Generic tree decomposition-based solver.

    Examples:
        >>> samples = orang.OrangSolver().sample_ising({0: 1}, {(0, 1): -1})
        >>> samples.first.sample
        {0: -1, 1: -1}

    """
    parameters = None
    properties = None

    def __init__(self):
        self.parameters = {'max_samples': [],
                           'elimination_order': []}
        self.properties = {'max_treewidth': 25}

    def sample(self, bqm, max_samples=1, elimination_order=None):
        """Find ground states of a binary quadratic model.

        Args:
            bqm (:obj:`dimod.BinaryQuadraticModel`):
                A binary quadratic model.

            max_samples (int, optional, default=1):
                The maximum number of solutions to return.

            elimination_order (list, optional):
                The variable elimination order. If None is provided, the
                min-fill heuristic is used to generate one.

        Returns:
            :obj:`dimod.SampleSet`

        Raises:
            ValueError: The treewidth_ of the given bqm and elimination order
            cannot exceed 25.

        .. _treewidth: https://en.wikipedia.org/wiki/Treewidth

        """

        if not bqm:
            samples = np.empty((max_samples, 0), dtype=samples_dtype)
            energies = np.zeros((max_samples), dtype=energies_dtype)
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

        return dimod.SampleSet.from_samples((samples, elimination_order), bqm.vartype,
                                            energy=energies + bqm.offset)


class OrangSampler(dimod.Sampler):
    """Generic tree decomposition-based sampler.

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
