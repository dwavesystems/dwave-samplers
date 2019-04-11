from __future__ import absolute_import

import random

import dimod
import dwave_networkx as dnx
import numpy as np

from orang.orang import solve_coo

__all__ = ['OrangSolver', 'OrangSampler']


class OrangSolver(dimod.Sampler):
    """Generic tree decomposition-based solver.

    Examples:
        >>> sampler = orang.OrangSolver()
        >>> samples = sampler.sample_ising({}, {(0, 1): -1})

    """
    parameters = None
    properties = None

    def __init__(self):
        self.parameters = {'num_reads': [],
                           'elimination_order': []}
        self.properties = {}

    def sample(self, bqm, num_reads=1, elimination_order=None):
        """Find ground states of a binary quadratic model.

        Args:
            bqm (:obj:`dimod.BinaryQuadraticModel`):
                A binary quadratic model.

            num_reads (int, optional, default=1):
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
            return dimod.SampleSet.from_samples([], bqm.vartype, energy=[])

        if elimination_order is None:
            tw, elimination_order = dnx.min_fill_heuristic(bqm.adj)
        else:
            # this also checks the order against the bqm
            tw = dnx.elimination_order_width(bqm.adj, elimination_order)

        # developer note: we start getting bad_alloc errors above tw 25, this
        # should be fixed in the future
        if tw > 25:
            msg = ("maximum treewidth of 25 exceeded. To see your bqm's treewidth run\n"
                   ">>> import dwave_networkx as dnx\n"
                   ">>> dnx.elimination_order(bqm.adj, {})".format(elimination_order))
            raise ValueError(msg)

        num_variables = len(bqm)

        linear, quadratic, offset = bqm.to_numpy_vectors(variable_order=elimination_order)

        # we used the elimination order as the variable order from the bqm, so 
        # don't need to pass in an elimination order
        samples, energies = solve_coo(linear, quadratic, offset, bqm.vartype,
                                      complexity=tw+1.,
                                      max_solutions=num_reads)

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
        self.properties = {}

    def sample(self, bqm, num_reads=1, elimination_order=None, beta=3.0, marginals=True, seed=None):
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

            marginals (bool, optional, default=True):
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
            info = {'log_partition_function': 0.0}
            if marginals:
                info['variable_marginals'] = {}

                info['interaction_marginals'] = {}

            return dimod.SampleSet.from_samples([], bqm.vartype, energy=[], info=info)

        if elimination_order is None:
            tw, elimination_order = dnx.min_fill_heuristic(bqm.adj)
        else:
            # this also checks the order against the bqm
            tw = dnx.elimination_order_width(bqm.adj, elimination_order)

        # developer note: we start getting bad_alloc errors above tw 25, this should be fixed in the
        # future
        if tw > 25:
            msg = ("maximum treewidth of 25 exceeded. To see your bqm's treewidth run\n"
                   ">>> import dwave_networkx as dnx\n"
                   ">>> dnx.elimination_order(bqm.adj, {})".format(elimination_order))
            raise ValueError(msg)

        num_variables = len(bqm)

        Q = bqm.binary.to_numpy_matrix(variable_order=elimination_order)
        qdata = np.asarray(Q.flatten(order='C'), dtype=np.double)  # make sure C ordered an correct dtype

        variable_order = np.arange(num_variables, dtype=np.intc)  # to match the given elimination_order

        if bqm.vartype is dimod.SPIN:
            low = -1
        else:
            low = 0

        if seed is None:
            # pick a random seed
            seed = random.randint(0, (1 << 16 - 1))

        Q = bqm.binary.to_numpy_matrix(variable_order=elimination_order)
        qdata = np.asarray(Q.flatten(order='C'), dtype=np.double)  # make sure C ordered an correct dtype

        variable_order = np.arange(num_variables, dtype=np.intc)  # to match the given elimination_order

        samples, logpf, linear_marginals, quadratic_marginals, pairs \
            = sample(qdata, num_variables, low,
                     num_reads,
                     variable_order, tw+1.,
                     marginals,
                     beta,
                     seed)

        # this should be sped up in the future
        energies = [bqm.energy(dict(zip(elimination_order, spl))) for spl in samples]

        # dev note: this does not appear to be correct, more investigation needed
        # info = {'log_partition_function': logpf}
        info = {}

        if marginals:
            info['variable_marginals'] = dict(zip(elimination_order, linear_marginals))

            info['interaction_marginals'] = interaction_marginals = {}
            configs = (low, low), (1, low), (low, 1), (1, 1)
            for (i, j), probs in zip(pairs, quadratic_marginals):
                u = elimination_order[i]
                v = elimination_order[j]
                interaction_marginals[(u, v)] = dict(zip(configs, probs))

        return dimod.SampleSet.from_samples((samples, elimination_order), bqm.vartype,
                                            energy=energies, info=info)
