# distutils: language = c++
# cython: language_level = 3

from cython.operator cimport preincrement as inc, dereference as deref

from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set

import dimod
cimport dimod


__all__ = ['elimination_order_width', 'min_fill_heuristic']


ctypedef unordered_map[Py_ssize_t, unordered_set[Py_ssize_t]] adj_t

cdef adj_t _cybqm_to_adj(dimod.cyBQM_float64 cybqm):
    cdef adj_t adj

    # add the nodes, even if they don't have any edges
    cdef Py_ssize_t vi
    for vi in range(cybqm.num_variables()):
        adj[vi]

    it = cybqm.data().cbegin_quadratic()
    while it != cybqm.data().cend_quadratic():
        adj[deref(it).u].insert(deref(it).v)
        adj[deref(it).v].insert(deref(it).u)
        inc(it)

    return adj


cdef void _elim_adj(adj_t& adj, Py_ssize_t vi) except +:
    """Remove vi from adj and make its neighborhood a clique."""

    # make the neighborhood of vi a clique
    uit = adj[vi].begin()
    while uit != adj[vi].end():
        vit = uit
        inc(vit)  # no self-loops
        while vit != adj[vi].end():
            adj[deref(vit)].insert(deref(uit))
            adj[deref(uit)].insert(deref(vit))
            inc(vit)

        # remove vi from its neighbors
        adj[deref(uit)].erase(vi)

        inc(uit)

    # finally remove vi
    adj.erase(vi)


def elimination_order_width(bqm, order):
    """Calculates the width of the tree decomposition induced by a
    variable elimination order.

    order must contain exactly the variables of the bqm
    """

    if len(bqm) != len(order):
        raise ValueError("bqm and order must have the name variables")

    cdef dimod.cyBQM_float64 cybqm = dimod.as_bqm(bqm, dtype=float).data

    cdef adj_t adj = _cybqm_to_adj(cybqm)

    # if there is at least one node then the treewidth is at least 1
    cdef Py_ssize_t treewidth = cybqm.num_variables() > 0

    cdef Py_ssize_t vi
    for v in order:
        vi = cybqm.variables.index(v)

        if adj[vi].size() > treewidth:
            treewidth = adj[vi].size()

        _elim_adj(adj, vi)

    return treewidth


# dev note: adj is actually a const reference, but Cython does not like
# that (fixed in Cython 3.0)
cdef Py_ssize_t _min_num_edges(adj_t& adj):
    """Get the node that would need to add the fewest edges when eliminated.

    Only defined for len(adj) > 0.
    """
    # The goal is to go through each node in adj, and determine how many
    # edges we would need to add in order to eliminate it.

    # C++ lambdas don't work so well in Cython so we do 'min' the hard way...

    cdef Py_ssize_t min_num_edges = adj.size() * adj.size()
    cdef Py_ssize_t min_node  # our return value


    cdef Py_ssize_t num_edges
    cdef Py_ssize_t vi

    it = adj.begin()
    while it != adj.end():
        vi = deref(it).first

        # for all pairs of nodes in the neighborhood of vi, count the missing edges
        num_edges = 0
        uit = deref(it).second.begin()
        while uit != deref(it).second.end():
            vit = uit
            while vit != deref(it).second.end():
                if not adj[deref(vit)].count(deref(uit)):
                    num_edges += 1
                inc(vit)
            inc(uit)

        if num_edges < min_num_edges:
            min_num_edges = num_edges
            min_node = vi

        inc(it)

    return min_node


def min_fill_heuristic(bqm):
    """Compute an upper bound on the treewidth of the given bqm based on
    the min-fill heuristic for the elimination ordering.

    Args:
        bqm: a binary quadratic model

    Returns:
        A 2-tuple containing the bound on the treewidth and the elimination 
        order.

    """
    cdef dimod.cyBQM_float64 cybqm = dimod.as_bqm(bqm, dtype=float).data

    cdef adj_t adj = _cybqm_to_adj(cybqm)

    cdef vector[Py_ssize_t] order
    order.reserve(cybqm.num_variables())

    # if there is at least one node then the treewidth is at least 1
    cdef Py_ssize_t upper_bound = cybqm.num_variables() > 0

    cdef Py_ssize_t vi
    while adj.size():
        vi = _min_num_edges(adj)

        if adj[vi].size() > upper_bound:
            upper_bound = adj[vi].size()

        # remove vi from adj
        _elim_adj(adj, vi)

        order.push_back(vi)

    cdef Py_ssize_t i
    variables = []
    for i in range(order.size()):
        variables.append(cybqm.variables.at(order[i]))

    return upper_bound, variables
