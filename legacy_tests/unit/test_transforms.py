import unittest
import random
import itertools

import networkx as nx

import savanna


class TestCutState(unittest.TestCase):
    def test_inversion_functional_random(self):
        for __ in range(5):
            G = nx.gnp_random_graph(10, .1)

            # connect if disconnected
            comp_nodes = [next(iter(H)) for H in nx.connected_components(G)]
            if len(comp_nodes) > 1:
                for u, v in itertools.combinations(comp_nodes, 2):
                    G.add_edge(u, v)

            assert nx.is_connected(G)

            for __ in range(5):
                sample = {v: random.choice((0, 1)) for v in G}

                # seed it with at least one node to handle the flips
                new_sample = savanna.cut_to_state(G, savanna.state_to_cut(G, sample), node=0, val=sample[0])

                self.assertEqual(sample, new_sample)
