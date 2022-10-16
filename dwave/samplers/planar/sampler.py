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
import networkx
from dimod import *
from networkx import is_perfect_matching

from dwave.samplers.planar import *

__all__ = ["PlanarGraphSampler"]


class PlanarGraphSampler(dimod.Sampler, dimod.Initialized):
    """A (weighted) max-cut sampler, specifically for planar-graphs

    In "planar graphs, the Maximum-Cut Problem is dual to the route inspection
    problem (the problem of finding a shortest tour that visits each edge of a
    graph at least once), in the sense that the edges that do not belong to a
    maximum cut-set of a graph G are the duals of the edges that are doubled in
    an optimal inspection tour of the dual graph of G. The optimal inspection
    tour forms a self-intersecting curve that separates the plane into two
    subsets, the subset of points for which the winding number of the curve is
    even and the subset for which the winding number is odd; these two subsets
    form a cut that includes all of the edges whose duals appear an odd number
    of times in the tour. The route inspection problem may be solved in polynomial
    time, and this duality allows the maximum cut problem to also be solved in
    polynomial time for planar graphs."

    Ref: https://en.wikipedia.org/wiki/Maximum_cut

    """
    parameters = None
    """Keyword arguments accepted by the sampling methods."""
    properties = None
    """Values for parameters accepted by the sampling methods."""

    def __init__(self):
        self.parameters = {
            'pos': []
        }
        self.properties = {}

    def sample(self,
               bqm: BinaryQuadraticModel,
               pos: dict = None,
               **kwargs) -> SampleSet:
        """Sample from a binary quadratic model.

        Args:
            bqm: Binary quadratic model to be sampled.
            pos: Position for each node

        Examples:
            >>> import dimod
            >>> from dwave.samplers.planar import PlanarGraphSampler
            >>> bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
            >>> bqm.add_interaction('a', 'b', +1.0)
            >>> bqm.add_interaction('b', 'c', +1.0)
            >>> bqm.add_interaction('c', 'a', +1.0)
            >>> pos = {'a': (0, 0), 'b': (1, 0), 'c': (0, 1)}
            >>> sample = PlanarGraphSampler().sample(bqm, pos)
            >>> sample.first
            Sample(sample={'a': 1, 'b': -1, 'c': -1}, energy=0, num_occurrences=1)

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

        assert is_perfect_matching(dual, matching)

        cut = _dual_matching_to_cut(G, matching)
        state = _cut_to_state(G, cut)

        if bqm.vartype is not dimod.BINARY:
            state = {v: 2 * b - 1 for v, b in state.items()}

        ret = SampleSet.from_samples(state, bqm.vartype, 0)
        return ret


def _determine_pos(G: MultiGraph) -> dict:
    is_planar, P = networkx.check_planarity(G)
    if not is_planar:
        raise ValueError("The provided BQM does not yield a planar embedding")

    layout = nx.planar_layout(P)
    return {k: tuple(v) for k, v in layout.items()}


def _dual_matching_to_cut(G: MultiGraph, matching: set) -> set:
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


def _cut_to_state(G: MultiGraph, cut: set, node=None, val=0) -> dict:
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
