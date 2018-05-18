import dimod
import networkx as nx

from savanna.agreement_ising import bqm_to_agreement_graph
from savanna.planar_graphs import expanded_dual
from savanna.transforms import cut_to_state


def ground_state_sample(bqm, rotation_system):
    """todo"""

    G, offset = bqm_to_agreement_graph(bqm)

    dual = expanded_dual(G, rotation_system)

    matching = nx.max_weight_matching(dual)

    cut = set(G.edges)

    # need to get the cut from the matching
    for u, v in matching:
        cut.discard(u)
        cut.discard(v)

    state = cut_to_state(G, cut)

    if bqm.vartype is dimod.BINARY:
        return state
    else:
        return {v: 2 * b - 1 for v, b in state.items()}
