.. image:: https://img.shields.io/pypi/v/dwave-samplers.svg
    :target: https://pypi.python.org/pypi/dwave-samplers

.. image:: https://img.shields.io/pypi/pyversions/dwave-samplers.svg
    :target: https://pypi.python.org/pypi/dwave-samplers

.. image:: https://codecov.io/gh/dwavesystems/dwave-samplers/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/dwavesystems/dwave-samplers

.. image:: https://circleci.com/gh/dwavesystems/dwave-samplers.svg?style=svg
    :target: https://circleci.com/gh/dwavesystems/dwave-samplers


==============
dwave-samplers
==============

.. start_samplers_about

Ocean software provides a variety of quantum, classical, and quantum-classical
hybrid `dimod <https://docs.dwavequantum.com/en/latest/ocean/api_ref_dimod>`_
`samplers <https://docs.dwavequantum.com/en/latest/concepts/samplers.html>`_
that run either remotely (for example, in the
`Leap <https://cloud.dwavesys.com/leap/>`_ service) or locally on your CPU.

Supported Samplers
==================

*dwave-samplers* implements the following
`classical <https://docs.dwavequantum.com/en/latest/quantum_research/classical_intro.html>`_
algorithms for solving
`binary quadratic models <https://docs.dwavequantum.com/en/latest/concepts/models.html>`_
(BQM):

*   `Planar`_: an exact solver for planar Ising problems with no linear biases.
*   `Random`_: a sampler that draws uniform random samples.
*   `Simulated Annealing`_: a probabilistic heuristic for optimization and
    approximate Boltzmann sampling well suited to finding good solutions of
    large problems.
*   `Simulated Quantum Annealing`_: a heuristic for optimization or approximate
    sampling.
*   `Steepest Descent`_: a discrete analogue of gradient descent, often used in
    machine learning, that quickly finds a local minimum.
*   `Tabu`_: a heuristic that employs local search with methods to escape local
    minima.
*   `Tree Decomposition`_: an exact solver for problems with low treewidth.

Planar
------

There are polynomial-time algorithms for finding the ground state of a planar
Ising model [#]_.

.. [#] Nicol Schraudolph, Dmitry Kamenetsky. *Efficient Exact Inference in
    Planar Ising Models*.
    Advances in Neural Information Processing Systems 21 (NIPS 2008).

>>> from dwave.samplers import PlanarGraphSolver
>>> solver = PlanarGraphSolver()

Get the ground state of a planar Ising model

>>> h = {}
>>> J = {(0, 1): -1, (1, 2): -1, (0, 2): 1}

>>> sampleset = solver.sample_ising(h, J)

Random
------

Random samplers provide a useful baseline performance comparison. The variable
assignments in each sample are chosen by a coin flip.

>>> from dwave.samplers import RandomSampler
>>> sampler = RandomSampler()

Create a random binary quadratic model.

>>> import dimod
>>> bqm = dimod.generators.gnp_random_bqm(100, .5, 'BINARY')

Get the best 5 sample found in .1 seconds.

>>> sampleset = sampler.sample(bqm, time_limit=.1, max_num_samples=5)
>>> num_reads = sampleset.info['num_reads']  # the total number of samples generated

Simulated Annealing
-------------------

`Simulated annealing <https://en.wikipedia.org/wiki/Simulated_annealing>`_ can
be used for heuristic optimization or approximate Boltzmann sampling. The
*dwave-samplers* implementation approaches the equilibrium distribution by
performing updates at a sequence of decreasing temperatures, terminating at the
target `β`.\ [#]_ Each spin is updated once in a fixed (by default) or
randomized order per point per temperature according to a customizable (by default
Metropolis-Hastings) update. When the temperature
is low the target distribution concentrates, at equilibrium, over ground states
of the model. Samples are guaranteed to match the equilibrium for long, smooth
temperature schedules. Schedules are customizable so that the sampler can be used
for standard Metropolis, Gibbs, Block-Gibbs or reverse annealing dynamics.

.. [#] `β` represents the inverse temperature, `1/(k T)`, of a
    `Boltzmann distribution <https://en.wikipedia.org/wiki/Boltzmann_distribution>`_
    where `T` is the thermodynamic temperature in kelvin and `k` is Boltzmann's
    constant.

>>> from dwave.samplers import SimulatedAnnealingSampler
>>> sampler = SimulatedAnnealingSampler()

Create a random binary quadratic model.

>>> import dimod
>>> bqm = dimod.generators.gnp_random_bqm(100, .5, 'BINARY')

Sample using simulated annealing with both the default temperature schedule
and a custom one.

>>> sampleset = sampler.sample(bqm)
>>> sampleset = sampler.sample(bqm, beta_range=[.1, 4.2], beta_schedule_type='linear')

Simulated Quantum Annealing
---------------------------

Simulated quantum annealing can be used for heuristic optimization or approximate
sampling. The *dwave-samplers* implementation performs dynamics defined by a schedule
for the driver(transverse) and problem(diagonal, binary quadratic model) Hamiltonian terms.
Each spin is updated according to a model-appropriate update as the schedule is
stepped through in discretized time (by sweep).
Using these methods equilibrated (thermalized) quantum Boltzmann distributions can be
approached, or the QPU schedule can be provided to simulate by classical dynamics
the annealing process. Although out-of-equilibrium dynamics and dynamical timescales
cannot be simulated by these classical methods, some phenomena can be emulated beyond
what is possible with simulated annealing. For algorithm details see [1]

[1] https://doi.org/10.1038/s41467-021-20901-5

>>> from dwave.samplers import PathIntegralAnnealingSampler
>>> sampler = PathIntegralAnnealingSampler()  # or RotorModelAnnealingSampler()

Create a random binary quadratic model.

>>> import dimod
>>> bqm = dimod.generators.gnp_random_bqm(100, .5, 'BINARY')

Sample projected states from a quantum process with a linear schedule

>>> sampleset = sampler.sample(bqm, beta_schedule_type="custom", Hp_field=[0, 10],  Hd_field=[10, 0])


Steepest Descent
----------------

`Steepest descent <https://en.wikipedia.org/wiki/Gradient_descent>`_ is the
discrete analogue of gradient descent, but the best move is computed using a
local minimization rather rather than computing a gradient. The dimension along
which to descend is determined, at each step, by the variable flip that causes
the greatest reduction in energy.

Steepest descent is fast and effective for unfrustrated problems, but it can get
stuck in local minima.

The quadratic unconstrained binary optimization (QUBO)
`E(x, y) = x + y - 2.5 * x * y`, for example, has two local minima:
`(0, 0)` with an energy of `0` and `(1, 1)` with an energy of `-0.5`.

>>> from dwave.samplers import SteepestDescentSolver
>>> solver = SteepestDescentSolver()

Construct the QUBO:

>>> from dimod import Binaries
>>> x, y = Binaries(['x', 'y'])
>>> qubo = x + y - 2.5 * x * y

If the solver starts uphill from the global minimum, it takes the steepest path
and finds the optimal solution.

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

Tabu
----

`Tabu search <https://en.wikipedia.org/wiki/Tabu_search>`_ is a heuristic that
employs local search and can escape local minima by maintaining a "tabu list" of
recently explored states that it does not revisit. The length of this tabu list
is called the "tenure". *dwave-samplers* implements the
`MST2 multistart tabu search algorithm <https://link.springer.com/article/10.1023/B:ANOR.0000039522.58036.68>`_
for quadratic unconstrained binary optimization (QUBO) problems.

Each read of the tabu algorithm consists of many starts. The solver takes the
best non-tabu step repeatedly until it does not improve its energy any more.

>>> from dwave.samplers import TabuSampler
>>> sampler = TabuSampler()

Construct a simple problem.

>>> from dimod import Binaries
>>> a, b = Binaries(['a', 'b'])
>>> qubo = -.5 * a + b - a * b

Sample using both default and custom values of tenure and number of restarts.

>>> sampleset0 = sampler.sample(qubo)
>>> sampleset1 = sampler.sample(qubo, tenure=1, num_restarts=1)

Tree Decomposition
------------------

`Tree decomposition <https://en.wikipedia.org/wiki/Tree_decomposition>`_-based
solvers have a runtime that is exponential in the
`treewidth <https://en.wikipedia.org/wiki/Treewidth>`_ of the problem graph. For
problems with low treewidth, the solver can find ground states very quickly.
However, for even moderately dense problems, performance is very poor.

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

Because the BQM is a binary tree, it has a treewidth of 1 and can be solved
exactly.

>>> sampleset = solver.sample(bqm)
>>> print(sampleset)
   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 ... 62 energy num_oc.
0 -1 +1 +1 -1 -1 -1 -1 +1 +1 +1 +1 +1 +1 +1 +1 -1 -1 -1 ... +1  -62.5       1
['SPIN', 1 rows, 1 samples, 63 variables]

.. end_samplers_about

Installation
============

To install the core package:

.. code-block:: bash

    pip install dwave-samplers

License
=======

Released under the Apache License 2.0

Contributing
============

Ocean's `contributing guide <https://docs.dwavequantum.com/en/latest/ocean/contribute.html>`_
has guidelines for contributing to Ocean packages.

Release Notes
-------------

**dwave-samplers** makes use of `reno <https://docs.openstack.org/reno/>`_ to
manage its release notes.

When making a contribution to **dwave-samplers** that will affect users, create
a new release note file by running

.. code-block:: bash

    reno new your-short-descriptor-here

You can then edit the file created under ``releasenotes/notes/``.
Remove any sections not relevant to your changes.
Commit the file along with your changes.