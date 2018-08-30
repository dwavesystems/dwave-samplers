import dimod
import networkx as nx

from savanna.io.dimod import bqm_to_multigraph
from savanna.planar import rotation_from_coordinates, plane_triangulate, odd_edge_orientation
from savanna.planar import expanded_dual
from savanna.transforms import cut_to_state, dual_matching_to_cut


def ground_state_bqm(bqm, pos):
    """todo"""

    if len(bqm) < 3:
        raise ValueError("bqm must have at least three variables")

    G, off = bqm_to_multigraph(bqm)

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

    cut = dual_matching_to_cut(G, matching)

    state = cut_to_state(G, cut)

    if bqm.vartype is dimod.BINARY:
        return state
    else:
        return {v: 2 * b - 1 for v, b in state.items()}


# this will be present in NetworkX 2.1, but for now use it
def is_perfect_matching(G, matching):
    from collections import Counter

    if isinstance(matching, dict):
        matching = matching_dict_to_set(matching)

    if not nx.is_matching(G, matching):
        return False

    counts = Counter(sum(matching, ()))

    return all(counts[v] == 1 for v in G)
