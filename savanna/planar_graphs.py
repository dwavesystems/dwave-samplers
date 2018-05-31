import math
from collections import OrderedDict, Mapping, namedtuple

import networkx as nx


@nx.utils.decorators.not_implemented_for('directed')
def rotation_system_from_coordinates(G, pos):
    """

    G must be planar
    G can be multigraph, graph

    r[u][v] := the next node clockwise around u after v

    """
    if isinstance(pos, Mapping):
        _pos = pos

        def pos(v):
            return _pos[v]

    rotation_system = {}
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

        rotation_system[u] = OrderedDict((circle[i - 1], edge) for i, edge in enumerate(circle))

    return rotation_system


Edge = namedtuple('Edge', ['head', 'tail', 'key'])


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
    ik = Edge(i, k, G.add_edge(i, k, weight=0.0))

    rotation_system[k][(k, i, ik.key)] = rotation_system[k][(k, j, jk.key)]
    rotation_system[k][(k, j, jk.key)] = Edge(k, i, ik.key)

    rotation_system[i][_inverse_rotation_system(rotation_system, i, ij)] = ik
    rotation_system[i][ik] = ij


@nx.utils.decorators.not_implemented_for('directed', 'multigraph')
def plane_triangulation(G, rotation_system):
    """Returns a multigraph"""

    if len(G) < 3:
        raise ValueError("only defined for graphs with 3 or more nodes")

    # plane triangulation might create multiple edges between two nodes so we need a multigraph
    triG = nx.MultiGraph()  # triangulated graph

    # add all of the nodes and edges (with weights) from G
    triG.add_nodes_from(G.nodes(data=True))
    triG.add_edges_from(G.edges(data=True))

    G = triG

    rotation_system = {v: {Edge(i, j, 0): Edge(s, t, 0) for (i, j), (s, t) in rotation_system[v].items()}
                       for v in rotation_system}

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

    return G, rotation_system


# log = logging.getLogger(__name__)


# def planar_faces(planar_G, rotation_system):
#     """Determine the faces of a planar graph.

#     Returns:
#         list[set]: Each face is a set containing the edges that define the face.

#     """

#     faces = []

#     # get the set of all edges. If we walk around each face, we pass through each edge twice,
#     # once moving clockwise, once counterclockwise
#     edges = set((u, v) for v in planar_G for u in planar_G[v])

#     while edges:
#         # grab an edge, and remove it from the list of edges we need to visit
#         root_tail, root_head = tail, head = edges.pop()

#         log.debug('starting a face walk from (%s, %s)', root_tail, root_head)

#         if tail == head:
#             raise ValueError("no self-loops allowed")

#         # start the face
#         face = set()

#         # get the next step
#         next_ = _clockwise_step(tail, head, rotation_system)
#         log.debug('walked from %s -> %s, now heading  to %s', tail, head, next_)
#         face.add((head, next_))
#         edges.discard((head, next_))

#         tail, head = head, next_

#         # while we haven't returned to our original edge
#         while (head != root_head) or (tail != root_tail):
#             next_ = _clockwise_step(tail, head, rotation_system)
#             log.debug('walked from %s -> %s, now heading to %s', tail, head, next_)
#             face.add((head, next_))
#             edges.discard((head, next_))

#             tail, head = head, next_

#         faces.append(face)

#     log.debug('faces: %s', faces)

#     return faces


# def _clockwise_step(tail, head, rotation_system):
#     circle = rotation_system[head]
#     idx = (circle.index(tail) - 1) % len(circle)
#     return circle[idx]  # the next step


# def expanded_dual(planar_G, rotation_system):

#     edual = nx.Graph()

#     faces = planar_faces(planar_G, rotation_system)

#     for face in faces:
#         # add each edge as a node, keeping track of the face
#         edual.add_nodes_from(face, face=face)

#         edual.add_edges_from(itertools.combinations(face, 2), weight=0)

#     # now add the between-face edges
#     for u, v in planar_G.edges:
#         edual.add_edge((u, v), (v, u), weight=planar_G[u][v].get('weight', 0))

#     return edual





# def is_perfect_matching(G, matching):
#     """Decides whether the given set represents a valid perfect matching in
#     ``G``.

#     A *perfect matching* in a graph is a matching in which exactly one edge
#     is incident upon each vertex.

#     Parameters
#     ----------
#     G : NetworkX graph

#     matching : dict or set
#         A dictionary or set representing a matching. If a dictionary, it
#         must have ``matching[u] == v`` and ``matching[v] == u`` for each
#         edge ``(u, v)`` in the matching. If a set, it must have elements
#         of the form ``(u, v)``, where ``(u, v)`` is an edge in the
#         matching.

#     Returns
#     -------
#     bool
#         Whether the given set or dictionary represents a valid perfect
#         matching in the graph.
#     """
#     if isinstance(matching, dict):
#         matching = nx.algorithms.matching.matching_dict_to_set(matching)

#     if not nx.is_matching(G, matching):
#         return False

#     count = Counter(sum(matching, ()))

#     return all(count[v] == 1 for v in G)
