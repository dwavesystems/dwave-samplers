import unittest
import logging

import networkx as nx

import savanna


# logging.basicConfig(level=logging.DEBUG)


class TestPlanarFaces(unittest.TestCase):
    def test_paper_example(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 3), (1, 2), (2, 4), (3, 4), (0, 3), (0, 2)])

        rotation_system = {0: [1, 2, 3],
                           1: [0, 3, 2],
                           2: [1, 4, 0],
                           3: [1, 0, 4],
                           4: [2, 3]}

        faces = savanna.planar_faces(G, rotation_system)

        self.assertEqual(len(faces), 4)
        self.assertIn({(0, 3), (3, 1), (1, 0)}, faces)
        self.assertIn({(0, 1), (1, 2), (2, 0)}, faces)
        self.assertIn({(1, 3), (3, 4), (4, 2), (2, 1)}, faces)
        self.assertIn({(3, 0), (0, 2), (2, 4), (4, 3)}, faces)

    def test_square_with_leg(self):
        G = nx.cycle_graph(['north', 'east', 'south', 'west'])
        G.add_edge('west', 'superwest')

        rotation_system = {'north': ['east', 'west'],
                           'east': ['north', 'south'],
                           'south': ['west', 'east'],
                           'west': ['south', 'north', 'superwest'],
                           'superwest': ['west']}

        faces = savanna.planar_faces(G, rotation_system)

        # self.assertEqual({frozenset(face) for face in faces},
        #                  {frozenset({'north', 'east', 'south', 'west'}),
        #                   frozenset({'north', 'east', 'south', 'west', 'superwest'})})

    def test_one_edge(self):
        G = nx.path_graph(2)

        rotation_system = {0: [1], 1: [0]}

        faces = savanna.planar_faces(G, rotation_system)

        self.assertEqual(faces, [{(0, 1), (1, 0)}])  # only one face

    def test__clockwise_step_one_edge(self):
        rotation_system = {0: [1], 1: [0]}  # path_graph(2)

        next_ = savanna.planar_graphs._clockwise_step(0, 1, rotation_system)
        self.assertEqual(next_, 0)

        next_ = savanna.planar_graphs._clockwise_step(1, 0, rotation_system)
        self.assertEqual(next_, 1)


class TestExpandedDual(unittest.TestCase):
    def test_paper_example(self):
        G = nx.Graph()
        G.add_edges_from([('a', 'b'), ('b', 'd'), ('b', 'c'), ('c', 'e'),
                          ('d', 'e'), ('a', 'd'), ('a', 'c')], weight=1)

        rotation_system = {'a': ['b', 'c', 'd'],
                           'b': ['a', 'd', 'c'],
                           'c': ['b', 'e', 'a'],
                           'd': ['b', 'a', 'e'],
                           'e': ['c', 'd']}

        H = savanna.expanded_dual(G, rotation_system)

        self.assertEqual(len(H), len(G.edges) * 2)

        self.assertIn(('a', 'b'), H)


class TestRotationSystemFromCoords(unittest.TestCase):
    def test_simple(self):
        G = nx.star_graph(4)

        pos = {0: (0, 0),
               1: (1, 0),  # north
               2: (0, 1),  # east
               3: (-1, 0),  # south
               4: (0, -1)}  # west

        rotation_system = savanna.rotation_system_from_coordinates(G, pos)

        self.assertEqual(rotation_system,
                         {0: [4, 1, 2, 3], 1: [0], 2: [0], 3: [0], 4: [0]})
