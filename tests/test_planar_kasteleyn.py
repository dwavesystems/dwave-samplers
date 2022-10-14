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
import numpy as np


def _half_kasteleyn(G):
    """G should be triangulated, oriented"""

    _ = sorted(G.edges(keys=True), key=lambda e: G.edges[e]['index'])
    num_edges = len(G.edges)

    # Create the empty half-kasteleyn
    H = np.zeros((2 * num_edges, 2 * num_edges), dtype=bool)

    for v in G.nodes:
        edge_s = next(iter(G.edges(v, keys=True)))
        s = G.edges[edge_s]['index']

        if G.edges[edge_s]['oriented'] == v:
            alpha = 2 * s
        else:
            alpha = 2 * s - 1

        edge_i = G.nodes[v]['rotation'][edge_s]
        i = G.edges[edge_i]['index']

        while True:
            if G.edges[edge_i]['oriented'] == v:
                H[2 * i - 1, alpha] = 1
                alpha = 2 * i

                if 'weight' not in G.edges[edge_i]:
                    # created by triangulation
                    H[2 * i - 1, 2 * i] = 1
            else:
                H[2 * i, alpha] = 1
                alpha = 2 * i - 1

            edge_i = G.nodes[v]['rotation'][edge_i]
            i = G.edges[edge_i]['index']

            if edge_i == G.nodes[v]['rotation'][edge_s]:
                break

    return H


def _half_kasteleyn_to_kasteleyn(G, H):
    K = H.astype(float)

    for edge in G.edges(keys=True):
        k = G.edges[edge]['index']
        weight = G.edges[edge].get('weight', 0.0)
        if weight:
            K[2 * k - 1, 2 * k] += np.exp(weight)

    K = K - np.transpose(K)

    return K


def _kasteleyn(G):
    return _half_kasteleyn_to_kasteleyn(G, _half_kasteleyn(G))


class TestHalfKasteleyn(unittest.TestCase):
    def test_three_path(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        G[0][1][0]['weight'] = 1.0
        G[1][2][0]['weight'] = .5
        # edge 0, 2 was added in the triangulation step

        for u, v, key in [(0, 1, 0), (1, 2, 0), (2, 0, 0)]:
            G[u][v][key]['oriented'] = v

        G.nodes[0]['pos'] = (0, 0)
        G.nodes[1]['pos'] = (1, 0)
        G.nodes[2]['pos'] = (0, 1)

        G.nodes[0]['rotation'] = OrderedDict([((0, 2, 0), (0, 1, 0)),
                                             ((0, 1, 0), (0, 2, 0))])
        G.nodes[1]['rotation'] = OrderedDict([((1, 0, 0), (1, 2, 0)),
                                             ((1, 2, 0), (1, 0, 0))])
        G.nodes[2]['rotation'] = OrderedDict([((2, 1, 0), (2, 0, 0)),
                                             ((2, 0, 0), (2, 1, 0))])

        # index the edges
        for idx, (u, v, k) in enumerate(G.edges(keys=True)):
            G[u][v][k]['index'] = idx

        _ = _half_kasteleyn(G)

        # todo


# noinspection DuplicatedCode,PyPep8Naming
class TestKasteleyn(unittest.TestCase):
    def test_three_path(self):
        G = nx.cycle_graph(3, create_using=nx.MultiGraph())

        G[0][1][0]['weight'] = 1.0
        G[1][2][0]['weight'] = .5
        # edge 0, 2 was added in the triangulation step

        for u, v, key in [(0, 1, 0), (1, 2, 0), (2, 0, 0)]:
            G[u][v][key]['oriented'] = v

        G.nodes[0]['pos'] = (0, 0)
        G.nodes[1]['pos'] = (1, 0)
        G.nodes[2]['pos'] = (0, 1)

        G.nodes[0]['rotation'] = OrderedDict([((0, 2, 0), (0, 1, 0)),
                                             ((0, 1, 0), (0, 2, 0))])
        G.nodes[1]['rotation'] = OrderedDict([((1, 0, 0), (1, 2, 0)),
                                             ((1, 2, 0), (1, 0, 0))])
        G.nodes[2]['rotation'] = OrderedDict([((2, 1, 0), (2, 0, 0)),
                                             ((2, 0, 0), (2, 1, 0))])

        # index the edges
        for idx, (u, v, k) in enumerate(G.edges(keys=True)):
            G[u][v][k]['index'] = idx

        _ = _kasteleyn(G)

        # todo
