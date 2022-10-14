import dimod
from dimod import *
from networkx import is_perfect_matching

from dwave.samplers.planar import *
from dwave.samplers.planar.transforms import *


__all__ = ["PlanarGraphSampler"]


# noinspection PyPep8Naming,DuplicatedCode
class PlanarGraphSampler(dimod.Sampler, dimod.Initialized):
    properties = None
    parameters = None

    def __init__(self):
        self.parameters = {
            'pos': []
        }
        self.properties = {}

    def sample(self,
               bqm: BinaryQuadraticModel,
               pos: dict = None,
               **parameters) -> SampleSet:
        """

        :param bqm:
        :param pos:
        :param parameters:
        :return:
        """

        if pos is None:
            pos = {}

        if len(bqm) < 3:
            raise ValueError("bqm must have at least three variables")

        G, off = bqm_to_multigraph(bqm)

        # apply the rotation system
        r = rotation_from_coordinates(G, pos)
        nx.set_node_attributes(G, name='rotation', values=r)

        # triangulation
        plane_triangulate(G)

        # create an edge indexing scheme
        indices = {edge: idx for idx, edge in enumerate(G.edges(keys=True))}
        nx.set_edge_attributes(G, name='index', values=indices)

        dual = expanded_dual(G)

        matching = nx.max_weight_matching(dual, maxcardinality=True, weight='weight')

        assert is_perfect_matching(dual, matching)

        cut = dual_matching_to_cut(G, matching)

        state = cut_to_state(G, cut)

        if bqm.vartype is not dimod.BINARY:
            state = {v: 2 * b - 1 for v, b in state.items()}

        ret = dimod.SampleSet.from_samples(state, bqm.vartype, 0)
        return ret
