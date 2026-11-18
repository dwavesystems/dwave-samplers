# Copyright 2022 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS F ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import copy
import typing

import dimod
import networkx as nx
import numpy as np

from dimod.typing import Variable

from dwave.samplers.planar.partition import log_partition_bqm
from dwave.samplers.planar.planar import rotation_from_coordinates, plane_triangulate, expanded_dual
from dwave.samplers.planar.util import bqm_to_multigraph

__all__ = ["PlanarGraphSolver", "PlanarGraphSampler"]


class PlanarGraphSolver(dimod.Sampler):
    """An exact solver for planar Ising problems with no linear biases."""
    parameters: typing.Dict[str, typing.Sequence[str]] = dict(pos=tuple())
    """Keyword arguments accepted by the sampling methods."""
    properties: typing.Dict[str, typing.Any] = dict()
    """Values for parameters accepted by the sampling methods."""

    def __init__(self):
        self.parameters = copy.deepcopy(self.parameters)
        self.properties = copy.deepcopy(self.properties)

    def sample(self,
               bqm: dimod.BinaryQuadraticModel,
               pos: typing.Optional[typing.Mapping[Variable, typing.Tuple[float, float]]] = None,
               **kwargs) -> dimod.SampleSet:
        """Find a ground state of a planar Ising problem without linear biases.

        Args:
            bqm: Binary quadratic model representing a planar Ising model with no linear biases.
            pos: Position for each node as a dict of form ``{node: (x-coordinate, y-coordinate),...}``.

        Examples:
            >>> import dimod
            >>> from dwave.samplers.planar import PlanarGraphSolver
            >>> bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
            >>> bqm.add_interaction('a', 'b', +1.0)
            >>> bqm.add_interaction('b', 'c', +1.0)
            >>> bqm.add_interaction('c', 'a', +1.0)
            >>> pos = {'a': (0, 0), 'b': (1, 0), 'c': (0, 1)}
            >>> sampleset = PlanarGraphSolver().sample(bqm, pos)
            >>> print(sampleset)
               a  b  c energy num_oc.
            0 +1 -1 -1   -1.0       1
            ['SPIN', 1 rows, 1 samples, 3 variables]

        """

        if len(bqm) < 3:
            raise ValueError("The provided BQM must have at least three variables")

        G, off = bqm_to_multigraph(bqm)

        if pos is None:
            pos = _determine_pos(G)

        # apply the rotation system
        r = rotation_from_coordinates(G, pos)
        nx.set_node_attributes(G, name='rotation', values=r)

        # triangulation
        plane_triangulate(G)

        # create an edge indexing scheme
        indices = {edge: idx for idx, edge in enumerate(G.edges(keys=True))}
        nx.set_edge_attributes(G, name='index', values=indices)

        dual = expanded_dual(G)
        matching = nx.max_weight_matching(dual, maxcardinality=True, weight='weight')

        assert nx.is_perfect_matching(dual, matching)

        cut = _dual_matching_to_cut(G, matching)
        state = _cut_to_state(G, cut)

        if bqm.vartype is not dimod.BINARY:
            state = {v: 2 * b - 1 for v, b in state.items()}

        return dimod.SampleSet.from_samples_bqm(state, bqm)


class PlanarGraphSampler(dimod.Sampler):
    """A fair sampler for planar Ising problems with no linear biases.

    Draws exact (fair) samples from the Boltzmann distribution of a planar Ising
    model using a contraction-tree scheme: an edge is contracted under the two
    relative-orientation hypotheses (agree/disagree), and the resulting exact
    partition functions (computed with the Kasteleyn method) give the
    conditional probability of that orientation. Iterating over a spanning tree
    of each connected component fixes a state up to the per-component global
    spin-flip symmetry, which is then broken uniformly at random.

    The log partition function is returned in the sample set's
    :attr:`~dimod.SampleSet.info` field.
    """
    parameters: typing.Dict[str, typing.Sequence[str]] = dict(
        num_reads=tuple(), beta=tuple(), marginals=tuple(), seed=tuple())
    """Keyword arguments accepted by the sampling methods."""
    properties: typing.Dict[str, typing.Any] = dict()
    """Values for parameters accepted by the sampling methods."""

    def __init__(self):
        self.parameters = copy.deepcopy(self.parameters)
        self.properties = copy.deepcopy(self.properties)

    def sample(self,
               bqm: dimod.BinaryQuadraticModel,
               num_reads: int = 1,
               beta: float = 1.0,
               marginals: bool = True,
               seed: typing.Optional[int] = None,
               **kwargs) -> dimod.SampleSet:
        """Draw fair samples from a planar Ising problem without linear biases.

        Args:
            bqm: Binary quadratic model representing a planar Ising model with no linear biases.
            num_reads: Number of samples to draw.
            beta: `Boltzmann distribution`_ inverse temperature parameter.
            marginals: Whether to compute the marginals. If True, they are
                included in the returned :class:`~dimod.SampleSet`'s ``info`` field.
            seed: Seed for the pseudo-random number generator.

        Returns:
            A :class:`~dimod.SampleSet` with ``num_reads`` samples. The ``info``
            field contains:

                * ``'log_partition_function'``: The log partition function at
                  inverse temperature ``beta``.

            If ``marginals=True``, also contains:

                * ``'correlations'``: Dict of the form ``{(u, v): c, ...}``, where
                  ``(u, v)`` is an interaction in the binary quadratic model and
                  ``c = sum_{s, t} s * t * prob(u == s & v == t)`` is the
                  two-point correlation, computed by the edge-expectation method.

        Examples:
            >>> import dimod
            >>> from dwave.samplers.planar import PlanarGraphSampler
            >>> bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
            >>> bqm.add_interaction('a', 'b', +1.0)
            >>> bqm.add_interaction('b', 'c', +1.0)
            >>> bqm.add_interaction('c', 'a', +1.0)
            >>> sampleset = PlanarGraphSampler().sample(bqm, num_reads=10)
            >>> sampleset.info['log_partition_function']  # doctest: +SKIP
            2.0825402296228025

        .. _Boltzmann distribution: https://en.wikipedia.org/wiki/Boltzmann_distribution

        """
        if len(bqm) < 3:
            raise ValueError("The provided BQM must have at least three variables")

        rng = np.random.default_rng(seed)

        # work in spin space (the planar method requires zero linear biases);
        # scaling the whole model by beta realizes the inverse temperature
        work = bqm.spin.change_vartype(dimod.SPIN, inplace=False)
        if beta != 1.0:
            work.scale(beta)

        logZ = _log_partition(work)

        info: typing.Dict[str, typing.Any] = {'log_partition_function': logZ}
        if marginals:
            info['correlations'] = _correlations(
                work, bqm.quadratic, bqm.vartype)

        variables = list(bqm.variables)
        to_binary = bqm.vartype is dimod.BINARY

        samples = np.empty((num_reads, len(variables)), dtype=np.int8)
        for i in range(num_reads):
            state = _sample_state(work, rng)
            for j, v in enumerate(variables):
                s = state[v]
                samples[i, j] = (s + 1) // 2 if to_binary else s

        energies = bqm.energies((samples, variables))

        return dimod.SampleSet.from_samples((samples, variables),
                                            bqm.vartype, energy=energies,
                                            info=info)


def _determine_pos(G: nx.MultiGraph) -> dict:
    is_planar, P = nx.check_planarity(G)
    if not is_planar:
        raise ValueError("The provided BQM does not yield a planar embedding")

    return nx.planar_layout(P)


def _dual_matching_to_cut(G: nx.MultiGraph, matching: set) -> set:
    """
    Then for each edge u,v in G, the expanded dual has two nodes (u,v), (v,u) with
    an edge between them. The matching is defined on that graph. There are
    additional edges, but we don't care about those
    """
    cut = set(G.edges)

    # need to get the cut from the matching
    for u, v in matching:
        if u[0] == v[1] and u[1] == v[0] and u[2] == v[2]:
            cut.discard(u)
            cut.discard(v)

    return cut


def _cut_to_state(G: nx.MultiGraph, cut: set, node=None, val=0) -> dict:
    if node is None:
        node = next(iter(G))  # get any node and assign it to val

    state = {}

    def _cut_to_state0(v, s):
        if v in state:
            assert state[v] == s, "Inconsistent state"
        else:
            state[v] = s
            for u in G[v]:
                for key in G[u][v]:
                    if (u, v, key) in cut or (v, u, key) in cut:
                        _cut_to_state0(u, 1 - s)
                    else:
                        _cut_to_state0(u, s)

    _cut_to_state0(node, val)

    return state


def _interaction_graph(bqm: dimod.BinaryQuadraticModel) -> nx.Graph:
    # exactly-zero couplings (e.g. from exact cancellation during contraction)
    # carry no interaction and must not connect otherwise-independent spins
    g = nx.Graph()
    g.add_nodes_from(bqm.variables)
    for (u, v), bias in bqm.quadratic.items():
        if bias != 0.0:
            g.add_edge(u, v)
    return g


def _restrict(bqm: dimod.BinaryQuadraticModel, component: typing.Iterable) -> dimod.BinaryQuadraticModel:
    """Build the sub-model on ``component`` using only its nonzero couplings."""
    cset = set(component)
    sub = dimod.BinaryQuadraticModel(dimod.SPIN)
    sub.add_linear_from((v, 0.0) for v in component)
    for (u, v), bias in bqm.quadratic.items():
        if bias != 0.0 and u in cset and v in cset:
            sub.set_quadratic(u, v, bias)
    return sub


def _log_partition(bqm: dimod.BinaryQuadraticModel) -> float:
    """Log partition function of a spin model with zero linear biases.

    The partition function factorizes over the connected components of the
    interaction graph, and within each component over its biconnected
    components (blocks): the Kasteleyn method requires a biconnected planar
    graph, which contraction does not preserve. A single bridge (two-spin
    block) is evaluated in closed form, and larger blocks use the planar
    Kasteleyn method. Each block carries its own global spin-flip symmetry, so
    joining ``B`` blocks of a component at shared cut vertices over-counts that
    symmetry by a factor ``2 ** (B - 1)``.
    """
    logZ = -bqm.offset

    graph = _interaction_graph(bqm)
    for component in nx.connected_components(graph):
        if len(component) == 1:
            logZ += np.log(2)
            continue

        blocks = list(nx.biconnected_component_edges(graph.subgraph(component)))
        for block_edges in blocks:
            nodes = {v for edge in block_edges for v in edge}
            if len(nodes) == 2:
                u, v = nodes
                j = bqm.get_quadratic(u, v)
                logZ += np.log(2) + np.logaddexp(j, -j)
            else:
                block = dimod.BinaryQuadraticModel(dimod.SPIN)
                block.add_linear_from((v, 0.0) for v in nodes)
                for u, v in block_edges:
                    block.set_quadratic(u, v, bqm.get_quadratic(u, v))
                G, _ = bqm_to_multigraph(block)
                logZ += log_partition_bqm(block, _determine_pos(G))

        logZ -= (len(blocks) - 1) * np.log(2)

    return float(logZ)


def _correlations(bqm: dimod.BinaryQuadraticModel, interactions,
                  vartype: dimod.Vartype) -> dict:
    """Edge-expectation two-point correlations for a zero-field spin model.

    ``bqm`` is a spin model whose couplings already include the
    inverse-temperature factor. For each interaction the probability that the
    two spins agree is ``Z_agree / (Z_agree + Z_disagree)``. The correlation is
    ``sum_{s, t} s * t * prob(u == s & v == t)``; by spin-flip symmetry the two
    agreeing configurations share the product ``low * low`` and the two
    disagreeing ones share ``1 * low``, so the sum collapses to the closed form
    below.
    """
    low = -1 if vartype is dimod.SPIN else 0

    correlations = {}
    for u, v in interactions:
        agree = bqm.copy()
        agree.contract_variables(u, v)

        disagree = bqm.copy()
        disagree.flip_variable(v)
        disagree.contract_variables(u, v)

        logZ_agree = _log_partition(agree)
        logZ_disagree = _log_partition(disagree)

        p_agree = np.exp(logZ_agree - np.logaddexp(logZ_agree, logZ_disagree))
        p_disagree = 1.0 - p_agree

        correlations[(u, v)] = p_agree * (low ** 2 + 1) / 2 + p_disagree * low

    return correlations


def _sample_state(bqm: dimod.BinaryQuadraticModel, rng: np.random.Generator) -> dict:
    """Draw one fair spin sample via the contraction-tree scheme.

    ``bqm`` is a spin model with zero linear biases whose couplings already
    include the inverse-temperature factor.
    """
    state = {}

    graph = _interaction_graph(bqm)
    for component in nx.connected_components(graph):
        component = list(component)

        if len(component) == 1:
            state[component[0]] = int(rng.choice((-1, 1)))
            continue

        work = _restrict(bqm, component)
        root = component[0]
        signs = {root: 1}  # spin of each variable relative to the root

        # contract along a spanning tree so that every contraction merges an
        # actual edge into the root; merging non-adjacent vertices would break
        # the planarity the Kasteleyn method relies on
        order = list(nx.bfs_tree(graph.subgraph(component), root))

        for v in order[1:]:
            agree = work.copy()
            agree.contract_variables(root, v)

            disagree = work.copy()
            disagree.flip_variable(v)
            disagree.contract_variables(root, v)

            logZ_agree = _log_partition(agree)
            logZ_disagree = _log_partition(disagree)

            p_agree = np.exp(logZ_agree - np.logaddexp(logZ_agree, logZ_disagree))

            if rng.random() < p_agree:
                work, signs[v] = agree, 1
            else:
                work, signs[v] = disagree, -1

        # break the component's global spin-flip symmetry uniformly
        flip = int(rng.choice((-1, 1)))
        for v in component:
            state[v] = flip * signs[v]

    return state

