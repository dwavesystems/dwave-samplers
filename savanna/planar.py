import math
from collections import OrderedDict, Mapping, namedtuple

import networkx as nx

Edge = namedtuple('Edge', ['head', 'tail', 'key'])
"""Represents a (directed) edge in a Multigraph"""


def rotation_from_coordinates(G, pos):
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

    if isinstance(pos, Mapping):
        _pos = pos

        def pos(v):
            return _pos[v]

    rotation = {}
    for u in G.nodes:
        x0, y0 = pos(u)

        def angle(edge):
            __, v, *__ = edge  # generic for Graph or MultiGraph
            x, y = pos(v)
            return math.atan2(y - y0, x - x0)

        if isinstance(G, nx.MultiGraph):
            circle = sorted(G.edges(u, keys=True), key=angle)
        else:
            circle = sorted(G.edges(u), key=angle)

        circle = [Edge(*edge) for edge in circle]

        rotation[u] = OrderedDict((circle[i - 1], edge) for i, edge in enumerate(circle))

    return rotation


def _inverse_rotation_system(rotation_system, v, edge):
    for key, val in rotation_system[v].items():
        if val == edge:
            return key

    raise RuntimeError


def _insert_chord(ij, jk, G, rotation_system):
    """Insert a chord between i and k."""
    assert ij.tail == jk.head
    i, j, _ = ij
    j, k, _ = jk

    # because G is a Multigraph, G.add_edge returns the key
    ik = Edge(i, k, G.add_edge(i, k))  # do not add weight for added edge

    rotation_system[k][(k, i, ik.key)] = rotation_system[k][(k, j, jk.key)]
    rotation_system[k][(k, j, jk.key)] = Edge(k, i, ik.key)

    rotation_system[i][_inverse_rotation_system(rotation_system, i, ij)] = ik
    rotation_system[i][ik] = ij


def plane_triangulate(G):
    """Takes a multigraph"""

    if len(G) < 3:
        raise ValueError("only defined for graphs with 3 or more nodes")

    rotation_system = {v: G.node[v]['rotation'] for v in G}

    # following the notation from the paper
    for i in G.nodes:
        for edge in list(G.edges(i, keys=True)):
            i, j, _ = ij = Edge(*edge)

            j, k, _ = jk = rotation_system[j][(j, i, ij.key)]
            k, l, _ = kl = rotation_system[k][(k, j, jk.key)]

            # assert j in G[l] and len(G[l][j]) == 1
            while l != i:
                m, n, _ = rotation_system[l][(l, k, kl.key)]
                if (m, n) == (l, j):
                    break

                if i == k:
                    raise NotImplementedError
                _insert_chord(ij, jk, G, rotation_system)

                i, j, _ = ij = kl
                j, k, _ = jk = rotation_system[j][(j, i, ij.key)]
                k, l, _ = kl = rotation_system[k][(k, j, jk.key)]

    return


def _make_odd(uv, G, visited, orientation):

    G.remove_edge(*uv)

    u, v, _ = uv

    if v in visited:
        return True

    visited.add(v)

    odd = False

    while G.adj[v]:
        vw = Edge(*next(iter(G.edges(v, keys=True))))
        if _make_odd(vw, G, visited, orientation):
            orientation.add(Edge(vw.tail, vw.head, vw.key))
            odd = not odd
        else:
            orientation.add(vw)
    return odd


def odd_edge_orientation(G):
    """requires triangular multigraph"""

    G = G.copy()  # we will be modifying G in-place

    visited = set()

    rs = Edge(*next(iter(G.edges)))

    orientation = {rs}
    _make_odd(rs, G, visited, orientation)

    edge_attr = {}

    for u, v, key in orientation:
        edge_attr[(u, v, key)] = v

    return edge_attr
