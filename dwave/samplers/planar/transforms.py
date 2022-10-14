# noinspection PyPep8Naming
def dual_matching_to_cut(G, matching):
    # the for each edge u,v in G, the expanded dual has two nodes (u,v), (v,u) with an edge between
    # them. The matching is defined on that graph. There are additional edges, but we don't care
    # about those
    cut = set(G.edges)

    # need to get the cut from the matching
    for u, v in matching:
        if u[0] == v[1] and u[1] == v[0] and u[2] == v[2]:
            cut.discard(u)
            cut.discard(v)

    return cut


# noinspection PyPep8Naming
def cut_to_state(G, cut, node=None, val=0):
    if node is None:
        node = next(iter(G))  # get any node and assign it to val

    state = {}

    def _cut_to_state(v, s):
        if v in state:
            assert state[v] == s, "Inconsistent state"
        else:
            state[v] = s
            for u in G[v]:
                for key in G[u][v]:
                    if (u, v, key) in cut or (v, u, key) in cut:
                        _cut_to_state(u, 1 - s)
                    else:
                        _cut_to_state(u, s)

    _cut_to_state(node, val)

    return state


# noinspection PyPep8Naming
def state_to_cut(G, state):
    return set((u, v) for u, v in G.edges if state[u] != state[v])
