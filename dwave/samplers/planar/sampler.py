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

from dwave.samplers.planar.planar import rotation_from_coordinates, plane_triangulate, expanded_dual
from dwave.samplers.planar.util import bqm_to_multigraph

__all__ = ["PlanarGraphSolver"]


class PlanarGraphSolver(dimod.Sampler, dimod.Initialized):
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
               pos: typing.Optional[typing.Mapping[dimod.typing.Variable, typing.Tuple[float, float]]] = None,
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
            >>> sample = PlanarGraphSolver().sample(bqm, pos)
            >>> sample.first
            Sample(sample={'a': 1, 'b': -1, 'c': -1}, energy=-1.0, num_occurrences=1)

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
