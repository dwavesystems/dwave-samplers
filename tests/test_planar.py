import unittest

from collections import OrderedDict

import networkx as nx

import savanna
from savanna import Edge


class TestSetRotationFromCoords(unittest.TestCase):
    def test_multigraph_star(self):
        G = nx.star_graph(4, nx.MultiGraph())
        G.add_edge(0, 1)  # add another edge going north

        pos = {0: (0, 0),
               1: (1, 0),  # north
               2: (0, 1),  # east
               3: (-1, 0),  # south
               4: (0, -1)}  # west

        r = savanna.rotation_from_coordinates(G, pos)

        # the order of the two edges (0, 1) are not fixed
        self.assertIn(r, [{0: {(0, 1, 1): (0, 2, 0),
                               (0, 1, 0): (0, 1, 1),
                               (0, 2, 0): (0, 3, 0),
                               (0, 3, 0): (0, 4, 0),
                               (0, 4, 0): (0, 1, 0)},
                           1: {(1, 0, 0): (1, 0, 1),
                               (1, 0, 1): (1, 0, 0)},
                           2: {(2, 0, 0): (2, 0, 0)},
                           3: {(3, 0, 0): (3, 0, 0)},
                           4: {(4, 0, 0): (4, 0, 0)}
                           },
                          {0: {(0, 1, 0): (0, 1, 1),
                               (0, 1, 1): (0, 2, 0),
                               (0, 2, 0): (0, 3, 0),
                               (0, 3, 0): (0, 4, 0),
                               (0, 4, 0): (0, 1, 0)},
                           1: {(1, 0, 0): (1, 0, 1),
                               (1, 0, 1): (1, 0, 0)},
                           2: {(2, 0, 0): (2, 0, 0)},
                           3: {(3, 0, 0): (3, 0, 0)},
                           4: {(4, 0, 0): (4, 0, 0)}
                           }])

    def test_empty(self):
        G = nx.MultiGraph()
        pos = {}
        r = savanna.rotation_from_coordinates(G, pos)
        self.assertTrue(len(G) == 0)
        self.assertTrue(len(r) == 0)


class TestPlaneTriangulation(unittest.TestCase):
    def test_triangle(self):
        # should do nothing
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        G.node[0]['rotation'] = OrderedDict([(Edge(0, 2, 0), Edge(0, 1, 0)),
                                             (Edge(0, 1, 0), Edge(0, 2, 0))])
        G.node[1]['rotation'] = OrderedDict([(Edge(1, 0, 0), Edge(1, 2, 0)),
                                             (Edge(1, 2, 0), Edge(1, 0, 0))])
        G.node[2]['rotation'] = OrderedDict([(Edge(2, 1, 0), Edge(2, 0, 0)),
                                             (Edge(2, 0, 0), Edge(2, 1, 0))])

        savanna.plane_triangulate(G)

        # G should not have changed
        self.assertEqual(set(G.edges()), set(nx.cycle_graph(3).edges()))
        self.assertEqual(set(G.nodes()), set(nx.cycle_graph(3).nodes()))

    def test_three_path(self):
        # should add an edge between 0, 2
        G = nx.path_graph(3, create_using=nx.MultiGraph())

        G.node[0]['rotation'] = OrderedDict([(Edge(0, 1, 0), Edge(0, 1, 0))])
        G.node[1]['rotation'] = OrderedDict([(Edge(1, 0, 0), Edge(1, 2, 0)),
                                             (Edge(1, 2, 0), Edge(1, 0, 0))])
        G.node[2]['rotation'] = OrderedDict([(Edge(2, 1, 0), Edge(2, 1, 0))])

        savanna.plane_triangulate(G)

        # MG should match now be a three cycle
        self.assertEqual(set(nx.cycle_graph(3).edges()), set(G.edges()))
        self.assertEqual(set(nx.cycle_graph(3).nodes()), set(G.nodes()))

    def test_paper_example(self):
        G = nx.cycle_graph([1, 2, 3, 4], create_using=nx.MultiGraph())
        G.add_edge(4, 5)
        G.add_edge(4, 2)

        pos = {1: (-1, -1), 2: (+1, -1), 3: (+1, +1), 4: (-1, +1), 5: (-.5, -.5)}

        r = savanna.rotation_from_coordinates(G, pos)
        nx.set_node_attributes(G, name='rotation', values=r)

        savanna.plane_triangulate(G)

        # we could test that it is plane triangular here but we would need to write
        # a test. For now let's just use this as a smoke test
        self.assertTrue(nx.is_biconnected(G))  # we do know planar triangular are biconnected

    def test_four_path(self):
        G = nx.path_graph([1, 2, 3, 4], create_using=nx.MultiGraph())

        pos = {1: (-1, -1), 2: (+1, -1), 3: (+1, +1), 4: (-1, +1), 5: (-.5, -.5)}

        r = savanna.rotation_from_coordinates(G, pos)
        nx.set_node_attributes(G, name='rotation', values=r)

        savanna.plane_triangulate(G)

        # we could test that it is plane triangular here but we would need to write
        # a test. For now let's just use this as a smoke test
        self.assertTrue(nx.is_biconnected(G))  # we do know planar triangular are biconnected


class TestOddEdgeOrientation(unittest.TestCase):
    def test_triangle(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        orientation = savanna.odd_edge_orientation(G)

        for u, v, key in orientation:
            self.assertIn(u, G.adj)
            self.assertIn(v, G.adj[u])
            self.assertIn(key, G[u][v])

    def test_square(self):
        G = nx.MultiGraph()

        G.add_edges_from([(0, 1, 0), (0, 3, 0), (0, 2, 0), (0, 2, 1), (1, 2, 0), (2, 3, 0)])

        orientation = savanna.odd_edge_orientation(G)

        self.assertEqual(len(orientation), len(G.edges))
        for (u, v, key), n in orientation.items():
            self.assertEqual(n, v)
            self.assertIn(u, G.adj)
            self.assertIn(v, G.adj[u])
            self.assertIn(key, G[u][v])

            self.assertNotIn((v, u, key), orientation)
