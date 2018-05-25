import itertools
import logging
import math

from collections import Counter

import networkx as nx


log = logging.getLogger(__name__)


def planar_faces(planar_G, rotation_system):
    """Determine the faces of a planar graph.

    Returns:
        list[set]: Each face is a set containing the edges that define the face.

    """

    faces = []

    # get the set of all edges. If we walk around each face, we pass through each edge twice,
    # once moving clockwise, once counterclockwise
    edges = set((u, v) for v in planar_G for u in planar_G[v])

    while edges:
        # grab an edge, and remove it from the list of edges we need to visit
        root_tail, root_head = tail, head = edges.pop()

        log.debug('starting a face walk from (%s, %s)', root_tail, root_head)

        if tail == head:
            raise ValueError("no self-loops allowed")

        # start the face
        face = set()

        # get the next step
        next_ = _clockwise_step(tail, head, rotation_system)
        log.debug('walked from %s -> %s, now heading  to %s', tail, head, next_)
        face.add((head, next_))
        edges.discard((head, next_))

        tail, head = head, next_

        # while we haven't returned to our original edge
        while (head != root_head) or (tail != root_tail):
            next_ = _clockwise_step(tail, head, rotation_system)
            log.debug('walked from %s -> %s, now heading to %s', tail, head, next_)
            face.add((head, next_))
            edges.discard((head, next_))

            tail, head = head, next_

        faces.append(face)

    log.debug('faces: %s', faces)

    return faces


def _clockwise_step(tail, head, rotation_system):
    circle = rotation_system[head]
    idx = (circle.index(tail) - 1) % len(circle)
    return circle[idx]  # the next step


def expanded_dual(planar_G, rotation_system):

    edual = nx.Graph()

    faces = planar_faces(planar_G, rotation_system)

    for face in faces:
        # add each edge as a node, keeping track of the face
        edual.add_nodes_from(face, face=face)

        edual.add_edges_from(itertools.combinations(face, 2), weight=0)

    # now add the between-face edges
    for u, v in planar_G.edges:
        edual.add_edge((u, v), (v, u), weight=planar_G[u][v].get('weight', 0))

    return edual


def rotation_system_from_coordinates(G, pos):
    """assume G is planar"""
    if callable(pos):
        pos = {v: pos(v) for v in G}

    # dev note: we use G.adj so that G can be a bqm

    rotation_system = {}
    for u in G.adj:
        x0, y0 = pos[u]

        def angle(v):
            x, y = pos[v]
            return math.atan2(y - y0, x - x0)

        rotation_system[u] = sorted(G.adj[u], key=angle)

    return rotation_system


def is_perfect_matching(G, matching):
    """Decides whether the given set represents a valid perfect matching in
    ``G``.

    A *perfect matching* in a graph is a matching in which exactly one edge
    is incident upon each vertex.

    Parameters
    ----------
    G : NetworkX graph

    matching : dict or set
        A dictionary or set representing a matching. If a dictionary, it
        must have ``matching[u] == v`` and ``matching[v] == u`` for each
        edge ``(u, v)`` in the matching. If a set, it must have elements
        of the form ``(u, v)``, where ``(u, v)`` is an edge in the
        matching.

    Returns
    -------
    bool
        Whether the given set or dictionary represents a valid perfect
        matching in the graph.
    """
    if isinstance(matching, dict):
        matching = nx.algorithms.matching.matching_dict_to_set(matching)

    if not nx.is_matching(G, matching):
        return False

    count = Counter(sum(matching, ()))

    return all(count[v] == 1 for v in G)
