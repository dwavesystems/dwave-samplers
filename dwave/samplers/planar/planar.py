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

import math
from collections import OrderedDict, namedtuple

import networkx as nx

Edge = namedtuple('Edge', ['head', 'tail', 'key'])
"""Represents a (directed) edge in a Multigraph"""


def rotation_from_coordinates(G: nx.MultiGraph, pos: dict) -> dict:
    """Compute the rotation system for a planar G from the node positions.

    Args:
        G (:obj:`networkx.MultiGraph`):
            A planar MultiGraph.

        pos (callable/dict):
            The position for each node.

    Returns:
        dict: Dictionary of nodes with a rotation dict as the value.

    """
    if not isinstance(G, nx.MultiGraph):
        raise TypeError("expected G to be a MultiGraph")

    rotation = {}
    for u in G.nodes:
        x0, y0 = pos[u]

        def angle(edge):
            if len(edge) == 3:  # generic for Graph or MultiGraph
                _, v, _ = edge
            else:
                _, v = edge
            x, y = pos[v]
            return math.atan2(y - y0, x - x0)

        if isinstance(G, nx.MultiGraph):
            circle = sorted(G.edges(u, keys=True), key=angle)
        else:
            circle = sorted(G.edges(u), key=angle)

        circle = [Edge(*edge) for edge in circle]

        rotation[u] = OrderedDict((circle[i - 1], edge) for i, edge in enumerate(circle))

    return rotation


def plane_triangulate(G: nx.MultiGraph):
    """Add edges to planar graph G to make it plane triangulated.

    An embedded graph is plane triangulated iff it is biconnected and
    each of its faces is a triangle.

    Args:
        G (:obj:`nx.MultiGraph`):
            A planar graph. Must have a rotation system. Note that edges are added in-place.

    """

    if len(G) < 3:
        raise ValueError("only defined for graphs with 3 or more nodes")

    rotation_system = {v: G.nodes[v]['rotation'] for v in G}

    # following the notation from the paper
    for i in G.nodes:
        for edge in list(G.edges(i, keys=True)):
            i, j, _ = ij = Edge(*edge)

            j, k, _ = jk = rotation_system[j][(j, i, ij.key)]
            k, l, _ = kl = rotation_system[k][(k, j, jk.key)]

            assert ij.tail == jk.head

            while l != i:
                m, n, _ = rotation_system[l][(l, k, kl.key)]
                if (m, n) == (l, j):
                    break

                if i == k:
                    # avoid self-loop
                    i, j, _ = ij = jk
                    j, k, _ = jk = kl
                    k, l, _ = kl = rotation_system[k][(k, j, jk.key)]

                    assert ij.tail == jk.head

                _insert_chord(ij, jk, G, rotation_system)

                i, j, _ = ij = kl
                j, k, _ = jk = rotation_system[j][(j, i, ij.key)]
                k, l, _ = kl = rotation_system[k][(k, j, jk.key)]

                assert ij.tail == jk.head

    assert is_plane_triangulated(G), "Something went wrong, G is not plane triangulated"

    return


def is_plane_triangulated(G: nx.MultiGraph) -> bool:
    """
    Args:
        G: The graph to be tested; note that this function expects G to have the rotation system as node attributes
    Returns:
        True iff G is bi-connected and each face is triangular
    """
    if not nx.is_biconnected(G):
        return False

    # x
    # |\
    # | \
    # y--z
    for x in G.nodes:
        for xz in G.edges(x, keys=True):

            _ = x, y, xykey = G.nodes[x]['rotation'][xz]
            _ = y, z, yzkey = G.nodes[y]['rotation'][(y, x, xykey)]
            _ = z, x, xzkey = G.nodes[z]['rotation'][(z, y, yzkey)]

            if xz != (x, z, xzkey):
                return False
    return True


def odd_in_degree_orientation(H: nx.MultiGraph) -> dict:
    """
    Args:
        H:
    Returns:
    """
    G = H.copy()

    orientation = set()

    for (u, v) in reversed(list(nx.dfs_edges(H))):
        uv = u, v, uvkey = u, v, min(G[u][v])

        uv_odd = True  # for now assume that we'll mark uv as odd

        # TODO: Should "keys=True" be "data=True"? "keys" seems to be unsupported?
        for vw in G.edges(v, keys=True):
            v, w, vwkey = vw

            if vw == (v, u, uvkey):
                # we'll do this one last
                pass
            elif (w, v, vwkey) in orientation:
                # we've already done this one and it's heading in, so we need to toggle the edge
                # we walked in on
                uv_odd = not uv_odd
                continue
            else:
                # ok, it's not the edge we came in on, and it's not already been marked or it's already
                # going out, so let's just set it going out so it doesn't change our degree
                orientation.add(vw)

        if uv_odd:
            # we want uv oriented towards v
            orientation.add(uv)
        else:
            orientation.add((v, u, uvkey))

    return {(u, v, key): v for (u, v, key) in orientation}


def expanded_dual(G: nx.MultiGraph) -> nx.Graph:
    """
    Args:
        G: should be multigraph, triangulated, oriented, edges indexed
    """

    dual = nx.Graph()

    # first we add the edges of the dual that cross the edges of G
    # for an edge (u, v, key) oriented towards v, we adopt the convention
    # that the right-hand node is labelled (u, v, key) and the left-hand
    # node is (v, u, key).
    for edge in G.edges(keys=True):
        u = edge
        v = (edge[1], edge[0], edge[2])
        dual.add_edge(u, v, weight=G.edges[edge].get('weight', 0.0))

    # next we add the edges within each triangular face
    for n in G.nodes:
        # iterate through the edges around n
        for left in G.edges(n, keys=True):
            u, v, _ = left = Edge(*left)
            assert u == n
            s, t, _ = right = Edge(*G.nodes[u]['rotation'][left])
            assert s == u

            # we want to connect the node left (from n) or right
            dual.add_edge(tuple(left), (t, s, right.key), weight=0.0)

    return dual


def _inverse_rotation_system(rotation_system: dict, v: str, edge: Edge) -> Edge:
    for e1, e2 in rotation_system[v].items():
        if e2 == edge:
            return e1

    raise RuntimeError


def _insert_chord(ij: Edge, jk: Edge, G: nx.MultiGraph, rotation_system: dict):
    """Insert a chord between i and k."""
    assert ij.tail == jk.head
    i, j, _ = ij
    j, k, _ = jk

    # because G is a Multigraph, G.add_edge returns the key
    ik = Edge(i, k, G.add_edge(i, k))

    rotation_system[k][(k, i, ik.key)] = rotation_system[k][(k, j, jk.key)]
    rotation_system[k][(k, j, jk.key)] = Edge(k, i, ik.key)

    rotation_system[i][_inverse_rotation_system(rotation_system, i, ij)] = ik
    rotation_system[i][ik] = ij
