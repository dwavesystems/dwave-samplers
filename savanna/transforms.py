

def cut_to_state(G, cut, node=None, val=0):
    if node is None:
        node = next(iter(G))  # get any node and assign it to val

    state = {}

    def _cut_to_state(v, s):
        if v in state:
            assert state[v] == s
        else:
            state[v] = s
            for u in G[v]:
                if (u, v) in cut or (v, u) in cut:
                    _cut_to_state(u, 1 - s)
                else:
                    _cut_to_state(u, s)

    _cut_to_state(node, val)

    return state


def state_to_cut(G, state):
    return set((u, v) for u, v in G.edges if state[u] != state[v])
