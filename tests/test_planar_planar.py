# Copyright 2022 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS F ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest
from collections import OrderedDict

import networkx as nx

from dwave.samplers.planar.planar import rotation_from_coordinates, Edge, plane_triangulate, odd_in_degree_orientation,\
    expanded_dual


class TestSetRotationFromCoords(unittest.TestCase):
    def test_multigraph_star(self):
        G = nx.star_graph(4, nx.MultiGraph())
        G.add_edge(0, 1)  # add another edge going north

        pos = {0: (0, 0),
               1: (1, 0),  # north
               2: (0, 1),  # east
               3: (-1, 0),  # south
               4: (0, -1)}  # west

        r = rotation_from_coordinates(G, pos)

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
        r = rotation_from_coordinates(G, pos)
        self.assertTrue(len(G) == 0)
        self.assertTrue(len(r) == 0)


class TestPlaneTriangulation(unittest.TestCase):
    def test_triangle(self):
        # should do nothing
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        G.nodes[0]['rotation'] = OrderedDict([(Edge(0, 2, 0), Edge(0, 1, 0)),
                                             (Edge(0, 1, 0), Edge(0, 2, 0))])
        G.nodes[1]['rotation'] = OrderedDict([(Edge(1, 0, 0), Edge(1, 2, 0)),
                                             (Edge(1, 2, 0), Edge(1, 0, 0))])
        G.nodes[2]['rotation'] = OrderedDict([(Edge(2, 1, 0), Edge(2, 0, 0)),
                                             (Edge(2, 0, 0), Edge(2, 1, 0))])

        plane_triangulate(G)

        # G should not have changed
        self.assertEqual(set(G.edges()), set(nx.cycle_graph(3).edges()))
        self.assertEqual(set(G.nodes()), set(nx.cycle_graph(3).nodes()))

    def test_three_path(self):
        # should add an edge between 0, 2
        G = nx.path_graph(3, create_using=nx.MultiGraph())

        G.nodes[0]['rotation'] = OrderedDict([(Edge(0, 1, 0), Edge(0, 1, 0))])
        G.nodes[1]['rotation'] = OrderedDict([(Edge(1, 0, 0), Edge(1, 2, 0)),
                                             (Edge(1, 2, 0), Edge(1, 0, 0))])
        G.nodes[2]['rotation'] = OrderedDict([(Edge(2, 1, 0), Edge(2, 1, 0))])

        plane_triangulate(G)

        # MG should match now be a three cycle
        self.assertEqual(set(nx.cycle_graph(3).edges()), set(G.edges()))
        self.assertEqual(set(nx.cycle_graph(3).nodes()), set(G.nodes()))

    def test_paper_example(self):
        G = nx.cycle_graph([1, 2, 3, 4], create_using=nx.MultiGraph())
        G.add_edge(4, 5)
        G.add_edge(4, 2)

        pos = {1: (-1, -1), 2: (+1, -1), 3: (+1, +1), 4: (-1, +1), 5: (-.5, -.5)}

        r = rotation_from_coordinates(G, pos)
        nx.set_node_attributes(G, name='rotation', values=r)

        plane_triangulate(G)

        # we could test that it is plane triangular here but we would need to write
        # a test. For now let's just use this as a smoke test
        self.assertTrue(nx.is_biconnected(G))  # we do know planar triangular are biconnected

    def test_four_path(self):
        G = nx.path_graph([1, 2, 3, 4], create_using=nx.MultiGraph())

        pos = {1: (-1, -1), 2: (+1, -1), 3: (+1, +1), 4: (-1, +1), 5: (-.5, -.5)}

        r = rotation_from_coordinates(G, pos)
        nx.set_node_attributes(G, name='rotation', values=r)

        plane_triangulate(G)

        # we could test that it is plane triangular here but we would need to write
        # a test. For now let's just use this as a smoke test
        self.assertTrue(nx.is_biconnected(G))  # we do know planar triangular are biconnected


class TestOddEdgeOrientation(unittest.TestCase):
    def test_triangle(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        orientation = odd_in_degree_orientation(G)

        # check that it's an orientation
        self.assertEqual(len(orientation), len(G.edges), "not every edge present in orientation")
        for u, v, key in orientation:
            self.assertIn(u, G.adj)
            self.assertIn(v, G.adj[u])
            self.assertIn(key, G[u][v])

        # apply the orientation and check that it's odd
        nx.set_edge_attributes(G, name='oriented', values=orientation)

        self.assertTrue(sum(bool(sum(data['oriented'] == v for _, _, data in G.edges(v, data=True)) % 2)
                        for v in G) >= len(G) - 1)

    def test_square(self):
        G = nx.MultiGraph()

        G.add_edges_from([(0, 1), (0, 3), (0, 2), (0, 2), (1, 2), (2, 3)])

        orientation = odd_in_degree_orientation(G)

        # check that it's an orientation
        self.assertEqual(len(orientation), len(G.edges), "not every edge present in orientation")
        for u, v, key in orientation:
            self.assertIn(u, G.adj)
            self.assertIn(v, G.adj[u])
            self.assertIn(key, G[u][v])

        # apply the orientation and check that it's odd
        nx.set_edge_attributes(G, name='oriented', values=orientation)

        self.assertTrue(sum(bool(sum(data['oriented'] == v for _, _, data in G.edges(v, data=True)) % 2)
                        for v in G) >= len(G) - 1)

    def test_grid_100x100(self):
        G = nx.MultiGraph()

        for x in range(100):
            for y in range(100):
                G.add_edge((x, y), (x + 1, y))
                G.add_edge((x, y), (x + 1, y + 1))
                G.add_edge((x, y), (x, y + 1))

        orientation = odd_in_degree_orientation(G)

        # check that it's an orientation
        self.assertEqual(len(orientation), len(G.edges), "not every edge present in orientation")
        for u, v, key in orientation:
            self.assertIn(u, G.adj)
            self.assertIn(v, G.adj[u])
            self.assertIn(key, G[u][v])

        # apply the orientation and check that it's odd
        nx.set_edge_attributes(G, name='oriented', values=orientation)

        self.assertTrue(sum(bool(sum(data['oriented'] == v for _, _, data in G.edges(v, data=True)) % 2)
                        for v in G) >= len(G) - 1)


class TestExpandedDual(unittest.TestCase):
    def test_three_cycle(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        G[0][1][0]['weight'] = 1.0
        G[1][2][0]['weight'] = .5
        # edge 0, 2 was added in the triangulation step, therfore has no weight

        for u, v, key in [(0, 1, 0), (1, 2, 0), (2, 0, 0)]:
            G[u][v][key]['oriented'] = v

        G.nodes[0]['rotation'] = OrderedDict([((0, 2, 0), (0, 1, 0)),
                                             ((0, 1, 0), (0, 2, 0))])
        G.nodes[1]['rotation'] = OrderedDict([((1, 0, 0), (1, 2, 0)),
                                             ((1, 2, 0), (1, 0, 0))])
        G.nodes[2]['rotation'] = OrderedDict([((2, 1, 0), (2, 0, 0)),
                                             ((2, 0, 0), (2, 1, 0))])

        dual = expanded_dual(G)

        self.assertEqual(len(dual.edges), 9)
        self.assertEqual(len(dual.nodes), 6)

        edges = [((1, 0, 0), (0, 1, 0)), ((1, 0, 0), (0, 2, 0)), ((1, 0, 0), (2, 1, 0)),
                 ((0, 1, 0), (2, 0, 0)),
                 ((0, 1, 0), (1, 2, 0)), ((0, 2, 0), (2, 0, 0)), ((0, 2, 0), (2, 1, 0)),
                 ((2, 0, 0), (1, 2, 0)),
                 ((2, 1, 0), (1, 2, 0))]

        for u, v in edges:
            self.assertIn(u, dual.adj)
            self.assertIn(v, dual.adj[u])

    def test_square(self):
        G = nx.cycle_graph(4, create_using=nx.MultiGraph())
        nx.set_edge_attributes(G, name='weight', values=1.0)
        G.add_edge(0, 2)
        G.add_edge(1, 3)

        # pos = {0: (0, 0), 1: (1, 0), 2: (1, 1), 3: (0, 1)}

        r = {0: OrderedDict([(Edge(head=0, tail=3, key=0), Edge(head=0, tail=2, key=0)),
                             (Edge(head=0, tail=1, key=0), Edge(head=0, tail=3, key=0)),
                             (Edge(head=0, tail=2, key=0), Edge(head=0, tail=1, key=0))]),
             1: OrderedDict([(Edge(head=1, tail=0, key=0), Edge(head=1, tail=2, key=0)),
                             (Edge(head=1, tail=2, key=0), Edge(head=1, tail=3, key=0)),
                             (Edge(head=1, tail=3, key=0), Edge(head=1, tail=0, key=0))]),
             2: OrderedDict([(Edge(head=2, tail=3, key=0), Edge(head=2, tail=1, key=0)),
                             (Edge(head=2, tail=1, key=0), Edge(head=2, tail=0, key=0)),
                             (Edge(head=2, tail=0, key=0), Edge(head=2, tail=3, key=0))]),
             3: OrderedDict([(Edge(head=3, tail=2, key=0), Edge(head=3, tail=0, key=0)),
                             (Edge(head=3, tail=0, key=0), Edge(head=3, tail=1, key=0)),
                             (Edge(head=3, tail=1, key=0), Edge(head=3, tail=2, key=0))])}
        nx.set_node_attributes(G, name='rotation', values=r)

        # orientation = {(1, 2, 0): 2, (0, 3, 0): 3, (1, 3, 0): 3, (2, 3, 0): 3,
        #                (0, 1, 0): 1, (2, 0, 0): 0}
        # nx.set_edge_attributes(G, name='oriented', values=orientation)

        dual = expanded_dual(G)

        self.assertEqual(len(dual.edges), 4 * 3 + len(G.edges))
        self.assertEqual(len(dual.nodes), 2 * len(G.edges))

        edges = [((0, 1), (1, 0)), ((0, 1), (3, 0)), ((0, 1), (1, 3)), ((1, 0), (0, 2)),
                 ((1, 0), (2, 1)), ((0, 3), (3, 0)), ((0, 3), (2, 0)), ((0, 3), (3, 2)),
                 ((3, 0), (1, 3)), ((0, 2), (2, 0)), ((0, 2), (2, 1)), ((2, 0), (3, 2)),
                 ((1, 2), (2, 1)), ((1, 2), (3, 1)), ((1, 2), (2, 3)), ((1, 3), (3, 1)),
                 ((3, 1), (2, 3)), ((2, 3), (3, 2))]

        for u, v in edges:
            u = (u[0], u[1], 0)
            v = (v[0], v[1], 0)
            self.assertIn(u, dual.adj)
            self.assertIn(v, dual.adj[u])
