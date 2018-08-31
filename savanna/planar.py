import itertools
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
            if len(edge) == 3:  # generic for Graph or MultiGraph
                _, v, _ = edge
            else:
                _, v = edge
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
    ik = Edge(i, k, G.add_edge(i, k))

    rotation_system[k][(k, i, ik.key)] = rotation_system[k][(k, j, jk.key)]
    rotation_system[k][(k, j, jk.key)] = Edge(k, i, ik.key)

    rotation_system[i][_inverse_rotation_system(rotation_system, i, ij)] = ik
    rotation_system[i][ik] = ij


def plane_triangulate(G):
    """Add edges to planar graph G to make it plane triangulated.

    An embedded graph is plane triangulated iff it is biconnected and
    each of its faces is a triangle.

    Args:
        G (:obj:`nx.MultiGraph`):
            A planar graph. Must have a rotation system. Note that edges are added in-place.

    """

    if len(G) < 3:
        raise ValueError("only defined for graphs with 3 or more nodes")

    rotation_system = {v: G.node[v]['rotation'] for v in G}

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


def is_plane_triangulated(G):
    # biconnected and each of its faces is a triangle
    # expects G to have the rotation system as node attributes, and is multigraph
    if not nx.is_biconnected(G):
        return False

    # x
    # |\
    # | \
    # y--z
    for x in G.nodes:
        for xz in G.edges(x, keys=True):

            xy = x, y, xykey = G.node[x]['rotation'][xz]
            yz = y, z, yzkey = G.node[y]['rotation'][(y, x, xykey)]
            zx = z, x, xzkey = G.node[z]['rotation'][(z, y, yzkey)]

            if xz != (x, z, xzkey):
                return False
    return True


def odd_in_degree_orientation(H):
    G = H.copy()

    visited = set()
    orientation = set()

    for (u, v) in reversed(list(nx.dfs_edges(H))):
        uv = u, v, uvkey = u, v, min(G[u][v])

        uv_odd = True  # for now assume that we'll mark uv as odd

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


def expanded_dual(G):
    """G should be multigraph, triangulated, oriented, edges indexed"""

    dual = nx.Graph()

    # first we add the edges of the dual that cross the edges of G
    # for an edge (u, v, key) oriented towards v, we adopt the convention that the right-hand node
    # is labelled (u, v, key) and the left-hand node is (v, u, key).
    for edge in G.edges(keys=True):
        u = edge
        v = (edge[1], edge[0], edge[2])
        dual.add_edge(u, v, weight=G.edges[edge].get('weight', 0.0))

    # next we add the edges within each triangular face
    for n in G.nodes:

        # iterate through the wedges around n
        for left in G.edges(n, keys=True):

            u, v, _ = left = Edge(*left)
            assert u == n
            s, t, _ = right = Edge(*G.node[u]['rotation'][left])
            assert s == u

            # we want to connect the node left (from n) or right
            dual.add_edge(tuple(left), (t, s, right.key), weight=0.0)

    return dual
