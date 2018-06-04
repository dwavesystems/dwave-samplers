import unittest

import networkx as nx

import savanna


class TestRotationSystemFromCoords(unittest.TestCase):
    def test_graph_star(self):
        G = nx.star_graph(4)

        pos = {0: (0, 0),
               1: (1, 0),  # north
               2: (0, 1),  # east
               3: (-1, 0),  # south
               4: (0, -1)}  # west

        r = savanna.rotation_system_from_coordinates(G, pos)

        self.assertEqual(r, {0: {(0, 1): (0, 2), (0, 2): (0, 3), (0, 3): (0, 4), (0, 4): (0, 1)},
                             1: {(1, 0): (1, 0)},
                             2: {(2, 0): (2, 0)},
                             3: {(3, 0): (3, 0)},
                             4: {(4, 0): (4, 0)}
                             })

    def test_multigraph_star(self):
        G = nx.star_graph(4, nx.MultiGraph())
        G.add_edge(0, 1)  # add another edge going north

        pos = {0: (0, 0),
               1: (1, 0),  # north
               2: (0, 1),  # east
               3: (-1, 0),  # south
               4: (0, -1)}  # west

        r = savanna.rotation_system_from_coordinates(G, pos)

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

#     def test_bqm_star(self):
#         pass
#         # G = nx.star_graph(4, nx.MultiGraph())
#         # G.add_edge(0, 1)  # add another edge going north

#         # pos = {0: (0, 0),
#         #        1: (1, 0),  # north
#         #        2: (0, 1),  # east
#         #        3: (-1, 0),  # south
#         #        4: (0, -1)}  # west

#         # r = savanna.rotation_system_from_coordinates(G, pos)

    # def test_single_node(self):
    #     G = nx.Graph()
    #     G.add_node(0)

    #     r = savanna.rotation_system_from_coordinates(G, {0: (0, 0)})

    #     self.assertEqual(r, {0: {}})


class TestPlaneTriangulation(unittest.TestCase):
    def test_triangle(self):
        # should do nothing
        G = nx.cycle_graph(3)

        pos = {0: (0, 0), 1: (1, 0), 2: (0, 1)}

        rotation_system = savanna.rotation_system_from_coordinates(G, pos)

        MG, rot = savanna.plane_triangulation(G, rotation_system)

        # G should not have changed
        self.assertEqual(set(G.edges()), set(nx.cycle_graph(3).edges()))
        self.assertEqual(set(G.nodes()), set(nx.cycle_graph(3).nodes()))

        # MG should match G
        self.assertEqual(set(G.edges()), set(MG.edges()))
        self.assertEqual(set(G.nodes()), set(MG.nodes()))

        self.assertEqual(rot, {w: {(u, v, 0): (s, t, 0) for (u, v), (s, t) in rotation_system[w].items()}
                               for w in G})

    def test_three_path(self):
        # should add an edge between 0, 2
        G = nx.path_graph(3)

        pos = {0: (0, 0), 1: (1, 0), 2: (0, 1)}

        rotation_system = savanna.rotation_system_from_coordinates(G, pos)

        MG, rot = savanna.plane_triangulation(G, rotation_system)

        # G should not have changed
        self.assertEqual(set(G.edges()), set(nx.path_graph(3).edges()))
        self.assertEqual(set(G.nodes()), set(nx.path_graph(3).nodes()))

        # MG should match now be a three cycke
        self.assertEqual(set(nx.cycle_graph(3).edges()), set(MG.edges()))
        self.assertEqual(set(nx.cycle_graph(3).nodes()), set(MG.nodes()))

    def test_paper_example(self):
        G = nx.cycle_graph([1, 2, 3, 4])  # , create_using=nx.MultiGraph())
        G.add_edge(4, 5)
        G.add_edge(4, 2)

        pos = {1: (-1, -1), 2: (+1, -1), 3: (+1, +1), 4: (-1, +1), 5: (-.5, -.5)}

        rotation_system = savanna.rotation_system_from_coordinates(G, pos)

        MG, rot = savanna.plane_triangulation(G, rotation_system)

        # we could test that it is planar_triangular here but we would need to write
        # a test. For now let's just use this as a smoke test

    def test__insert_chord(self):
        G = nx.path_graph(3, create_using=nx.MultiGraph())

        pos = {0: (0, 0), 1: (1, 0), 2: (0, 1)}

        rotation_system = savanna.rotation_system_from_coordinates(G, pos)

        ij = savanna.planar.Edge(0, 1, 0)
        jk = savanna.planar.Edge(1, 2, 0)

        savanna.planar._insert_chord(ij, jk, G, rotation_system)

        # G should now have an edge between i, k
        self.assertEqual(len(G.edges), 3)
        self.assertIn(2, G.adj[0])


class TestOddEdgeOrientation(unittest.TestCase):
    def test_triangle(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        oriented = savanna.odd_edge_orientation(G)

        for u, v, key in oriented.edges(keys=True):
            self.assertIn(u, G.adj)
            self.assertIn(v, G.adj[u])
            self.assertIn(key, G[u][v])
            self.assertNotIn((v, u, key), oriented.edges(keys=True))

        self.assertEqual(len(oriented.edges), len(G.edges))
