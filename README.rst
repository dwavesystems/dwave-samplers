.. image:: https://img.shields.io/pypi/v/dwave-samplers.svg
    :target: https://pypi.python.org/pypi/dwave-samplers

.. image:: https://img.shields.io/pypi/pyversions/dwave-samplers.svg
    :target: https://pypi.python.org/pypi/dwave-samplers

.. image:: https://codecov.io/gh/dwavesystems/dwave-samplers/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/dwavesystems/dwave-samplers

.. image:: https://circleci.com/gh/dwavesystems/dwave-samplers.svg?style=svg
    :target: https://circleci.com/gh/dwavesystems/dwave-samplers

dwave-samplers
==============

.. index-start-marker

Implements several classical algorithms for solving
`binary quadratic models <https://docs.ocean.dwavesys.com/en/stable/concepts/bqm.html>`_.

Algorithms
----------

The following algorithms are implemented in **dwave-samplers**:

* `Simulated Annealing <readme_simulated_annealing_>`_: a probabilistic heuristic for optimization and approximate Boltzmann sampling well suited to finding good solutions of large problems.
* `Steepest Descent <readme_steepest_descent_>`_: a discrete analogue of gradient descent, often used in machine learning, that quickly finds a local minimum.
* `Tabu <readme_tabu_>`_: a heuristic that employs local search with methods to escape local minima.
* `Tree Decomposition <readme_tree_decomposition_>`_: an exact solver for problems with low treewidth.

.. _readme_simulated_annealing:

Simulated Annealing
~~~~~~~~~~~~~~~~~~~

`Simulated annealing <https://en.wikipedia.org/wiki/Simulated_annealing>`_ can be used for approximate Boltzmann sampling or heuristic optimization. This implementation approaches the equilibrium distribution by performing updates at a sequence of increasing beta values, ``beta_schedule``, terminating at the target beta. Each spin is updated once in a fixed order per point in the ``beta_schedule`` according to a Metropolis- Hastings update. When beta is large the target distribution concentrates, at equilibrium, over ground states of the model. Samples are guaranteed to match the equilibrium for long 'smooth' beta schedules.

>>> from dwave.samplers import SimulatedAnnealingSampler
>>> sampler = SimulatedAnnealingSampler()

Create a random binary quadratic model (BQM).

>>> import dimod
>>> bqm = dimod.generators.gnp_random_bqm(100, .5, 'BINARY')

Sample using simulated annealing, both with the default beta range and a user-specified one.

>>> sampleset = sampler.sample(bqm)
>>> sampleset = sampler.sample(bqm, beta_range=[.1, 4.2])

.. _readme_steepest_descent:

Steepest Descent
~~~~~~~~~~~~~~~~

`Steepest descent <https://en.wikipedia.org/wiki/Gradient_descent>`_ is the discrete analogue of gradient descent, but the best move is computed using a local minimization rather rather than computing a gradient. At each step, we determine the dimension along which to descend based on the highest energy drop caused by a variable flip.

Steepest descent is fast and effective for unfrustrated problems, but it can get stuck in local minima. Consider a quadratic unconstrained binary optimization (QUBO) of the form `E(x, y) = x + y - 2.5 * x * y`. This problem has two local minima, `(0, 0)` which has an energy of `0`, and `(1, 1)` which has an energy of `-0.5`.

>>> from dwave.samplers import SteepestDescentSolver
>>> solver = SteepestDescentSolver()

Construct the QUBO:

>>> from dimod import Binaries
>>> x, y = Binaries(['x', 'y'])
>>> qubo = x + y - 2.5 * x * y

If the solver starts uphill from the global minimum, it takes the
steepest path, thereby finding the optimal solution.

>>> sampleset = solver.sample(qubo, initial_states={'x': 0, 'y': 1})
>>> print(sampleset)
   x  y energy num_oc. num_st.
0  1  1   -0.5       1       1
['BINARY', 1 rows, 1 samples, 2 variables]

If the solver starts in a local minimum, it gets stuck.

>>> sampleset = solver.sample(qubo, initial_states={'x': 0, 'y': 0})
>>> print(sampleset)
   x  y energy num_oc. num_st.
0  0  0    0.0       1       0
['BINARY', 1 rows, 1 samples, 2 variables]

.. _readme_tabu:

Tabu
~~~~

An implementation of the `MST2 multistart tabu search algorithm <https://link.springer.com/article/10.1023/B:ANOR.0000039522.58036.68>`_ for quadratic unconstrained binary optimization (QUBO) problems.

Tabu solvers maintain a "tabu list" of recently explored states. The algorithm will not revisit those states, thereby avoiding getting stuck in local minima. The length of the tabu list is called the "tenure".

Each read of the tabu algorithm consists of many starts. The solver will then take the best non-tabu step repeatedly until it does not improve its energy any more.

>>> from dwave.samplers import TabuSampler
>>> sampler = TabuSampler()

Construct a simple problem

>>> from dimod import Binaries
>>> a, b = Binaries(['a', 'b'])
>>> qubo = -.5 * a + b - a * b

We can sample using the ``TabuSampler``, either using the default tenure and number of restarts, or specifying them explicitly.

>>> sampleset0 = sampler.sample(qubo)
>>> sampleset1 = sampler.sample(qubo, tenure=1, num_restarts=1)

.. _readme_tree_decomposition:

Tree Decomposition
~~~~~~~~~~~~~~~~~~

`Tree decomposition <https://en.wikipedia.org/wiki/Tree_decomposition>`_-based solvers have a runtime that is exponential in the `treewidth <https://en.wikipedia.org/wiki/Treewidth>`_ of the problem graph. For problems with low treewidth, the solver can find ground states very quickly. However, the performance is very poor for even moderately dense problems.

>>> from dwave.samplers import TreeDecompositionSolver
>>> solver = TreeDecompositionSolver()

Construct a large, tree-shaped problem.

>>> import dimod
>>> import networkx as nx
>>> tree = nx.balanced_tree(2, 5)  # binary tree with a height of five
>>> bqm = dimod.BinaryQuadraticModel('SPIN')
>>> bqm.set_linear(0, .5)
>>> for u, v in tree.edges:
...     bqm.set_quadratic(u, v, 1)

Because the BQM is a binary tree, it has a tree width of 1 and can be solved exactly.

>>> sampleset = solver.sample(bqm)
>>> print(sampleset)
   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 ... 62 energy num_oc.
0 -1 +1 +1 -1 -1 -1 -1 +1 +1 +1 +1 +1 +1 +1 +1 -1 -1 -1 ... +1  -62.5       1
['SPIN', 1 rows, 1 samples, 63 variables]

.. index-end-marker

Installation
------------

To install the core package:

.. code-block:: bash

    pip install dwave-samplers

License
-------

Released under the Apache License 2.0

Contributing
------------

Ocean's `contributing guide <https://docs.ocean.dwavesys.com/en/stable/contributing.html>`_
has guidelines for contributing to Ocean packages.

Release Notes
~~~~~~~~~~~~~

**dwave-samplers** makes use of `reno <https://docs.openstack.org/reno/>`_ to manage its
release notes.

When making a contribution to **dwave-samplers** that will affect users, create a new
release note file by running

.. code-block:: bash

    reno new your-short-descriptor-here

You can then edit the file created under ``releasenotes/notes/``.
Remove any sections not relevant to your changes.
Commit the file along with your changes.
