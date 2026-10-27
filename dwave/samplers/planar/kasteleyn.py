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
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np

__all__ = ['kasteleyn']


def half_kasteleyn(G):
    """G should be triangulated, oriented"""

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


def half_kasteleyn_to_kasteleyn(G, H):
    K = H.astype(float)

    for edge in G.edges(keys=True):
        k = G.edges[edge]['index']
        weight = G.edges[edge].get('weight', 0.0)
        if weight:
            K[2 * k - 1, 2 * k] += np.exp(weight)

    K = K - np.transpose(K)

    return K


def kasteleyn(G):
    return half_kasteleyn_to_kasteleyn(G, half_kasteleyn(G))
