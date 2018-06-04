import unittest
from collections import OrderedDict

import networkx as nx

import savanna


class TestHalfKasteleyn(unittest.TestCase):
    def test_three_path(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        G[0][1][0]['weight'] = 1.0
        G[1][2][0]['weight'] = .5
        # edge 0, 2 was added in the triangulation step

        for u, v, key in [(0, 1, 0), (1, 2, 0), (2, 0, 0)]:
            G[u][v][key]['oriented'] = v

        G.node[0]['pos'] = (0, 0)
        G.node[1]['pos'] = (1, 0)
        G.node[2]['pos'] = (0, 1)

        G.node[0]['rotation'] = OrderedDict([((0, 2, 0), (0, 1, 0)),
                                             ((0, 1, 0), (0, 2, 0))])
        G.node[1]['rotation'] = OrderedDict([((1, 0, 0), (1, 2, 0)),
                                             ((1, 2, 0), (1, 0, 0))])
        G.node[2]['rotation'] = OrderedDict([((2, 1, 0), (2, 0, 0)),
                                             ((2, 0, 0), (2, 1, 0))])

        # index the edges
        for idx, (u, v, k) in enumerate(G.edges(keys=True)):
            G[u][v][k]['index'] = idx

        H = savanna.half_kasteleyn(G)

        # todo


class TestKasteleyn(unittest.TestCase):
    def test_three_path(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        G[0][1][0]['weight'] = 1.0
        G[1][2][0]['weight'] = .5
        # edge 0, 2 was added in the triangulation step

        for u, v, key in [(0, 1, 0), (1, 2, 0), (2, 0, 0)]:
            G[u][v][key]['oriented'] = v

        G.node[0]['pos'] = (0, 0)
        G.node[1]['pos'] = (1, 0)
        G.node[2]['pos'] = (0, 1)

        G.node[0]['rotation'] = OrderedDict([((0, 2, 0), (0, 1, 0)),
                                             ((0, 1, 0), (0, 2, 0))])
        G.node[1]['rotation'] = OrderedDict([((1, 0, 0), (1, 2, 0)),
                                             ((1, 2, 0), (1, 0, 0))])
        G.node[2]['rotation'] = OrderedDict([((2, 1, 0), (2, 0, 0)),
                                             ((2, 0, 0), (2, 1, 0))])

        # index the edges
        for idx, (u, v, k) in enumerate(G.edges(keys=True)):
            G[u][v][k]['index'] = idx

        K = savanna.kasteleyn(G)

        # todo
