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

import networkx as nx
from dimod import BinaryQuadraticModel


def bqm_to_multigraph(bqm: BinaryQuadraticModel):
    """
    Args:
        bqm: Binary quadratic model to be sampled.
    Returns:
        G - the multi-graph representing the provided BQM, where there is a
            weighted-edge ``uv`` for each quadratic-term, where the weight is
            a linear-function of its bias
        offset - the sum over the biases in the BQM plus its spin-offset
    """
    if any(bqm.spin.linear.values()):
        raise NotImplementedError("not yet implemented for non-zero linear biases")

    offset = bqm.spin.offset
    G = nx.MultiGraph()
    for (u, v), bias in bqm.spin.quadratic.items():
        G.add_edge(u, v, weight=-2 * bias)
        offset += bias

    return G, offset
