import itertools

try:
    from collections.abc import Mapping, Sequence
except ImportError:
    from collections import Mapping, Sequence

import dimod
import dwave_networkx as dnx
import networkx as nx

from orang.poly import solve

__all__ = 'solve_polynomial',


@dimod.decorators.vartype_argument('vartype')
def solve_polynomial(poly, vartype, num_reads=1, elimination_order=None):
    """Find states that minimize the energy of the given polynomial.

    Args:
        poly (dict):
            Polynomial as a dict of form {term: bias, ...}, where `term` is a tuple of
            variables and `bias` the associated bias.

        vartype (:class:`.Vartype`/str/set):
            Variable type for the binary quadratic model. Accepted input values:

            * :class:`.Vartype.SPIN`, ``'SPIN'``, ``{-1, 1}``
            * :class:`.Vartype.BINARY`, ``'BINARY'``, ``{0, 1}``

        num_reads (int, optional, default=1):
            Maximum number of unique solutions to return.

        elimination_order (list, optional):
            The variable elimination order. If not provided, the min-fill heuristic is used.

    Returns:
        :class:`dimod.SampleSet`

    """

    # need to determine the elimination order or get the complexity. All of the interactions
    # become cliques and we use the normal min-fill heuristic thereafter
    G = nx.Graph()
    G.add_nodes_from(v for interaction in poly for v in interaction)
    G.add_edges_from((u, v)
                     for interaction in poly
                     for (u, v) in itertools.combinations(interaction, 2)
                     if u != v)

    if elimination_order is None:
        tw, elimination_order = dnx.min_fill_heuristic(G)
    elif not isinstance(elimination_order, Sequence):
        raise TypeError("expected elimination_order to be a list")
    else:
        # this also checks the order against the polynomial
        tw = dnx.elimination_order_width(G, elimination_order)

    label_to_idx = {v: idx for idx, v in enumerate(elimination_order)}

    if vartype is dimod.SPIN:
        ordered_poly = _ordered_spin_polynomial(poly, label_to_idx)
    else:
        ordered_poly = _ordered_binary_polynomial(poly, label_to_idx)

    samples, energies = solve(ordered_poly, len(elimination_order), tw + 2.,
                              vartype,
                              num_reads)

    return dimod.SampleSet.from_samples((samples, elimination_order), vartype,
                                        energy=energies)


def _ordered_binary_polynomial(poly, label_to_idx):
    """In Binary space v^n == v so just make them unique"""
    return {tuple(sorted({label_to_idx[v] for v in interaction})): bias
            for interaction, bias in poly.items()}


def _ordered_spin_polynomial(poly, label_to_idx):
    ordered_poly = {}

    for interaction, bias in poly.items():
        ordered = tuple(sorted({label_to_idx[v] for v in interaction}))

        if len(ordered) != len(interaction):
            if not hasattr(interaction, "count"):
                interaction = list(interaction)

            orderedset = set()
            for v in interaction:
                if v in orderedset:
                    continue

                count = interaction.count(v)

                if count % 2:
                    orderedset.add(label_to_idx[v])
                # else: even number and we can skip since -1^2 = 1

            ordered = tuple(sorted(orderedset))

        if ordered in ordered_poly:
            ordered_poly[ordered] += bias
        else:
            ordered_poly[ordered] = bias

    return ordered_poly
