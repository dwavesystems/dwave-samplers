import networkx as nx


def bqm_to_multigraph(bqm):
    """todo"""
    if any(bqm.spin.linear.values()):
        raise NotImplementedError("not yet implemented for non-zero linear biases")

    offset = bqm.spin.offset
    G = nx.OrderedMultiGraph()
    for (u, v), bias in bqm.spin.quadratic.items():
        G.add_edge(u, v, weight=-2 * bias)

        offset += bias

    return G, offset


def agreement_energy(state, G, offset=0.0):
    """todo"""
    en = offset
    for u, v, data in G.edges(data=True):
        if state[u] != state[v]:
            en += data['weight']
    return en
