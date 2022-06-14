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
`binary quadratic models <https://docs.ocean.dwavesys.com/en/stable/concepts/bqm.html>`_

Algorithms
----------

The following algorithms are implemented in **dwave-samplers**.

Simulated Annealing
~~~~~~~~~~~~~~~~~~~
A simulated annealing sampler can be used for approximate Boltzmann sampling or heuristic optimization. This implementation approaches the equilibrium distribution by performing updates at a sequence of increasing beta values, beta_schedule, terminating at the target beta. Each spin is updated once in a fixed order per point in the beta_schedule according to a Metropolis- Hastings update. When beta is large the target distribution concentrates, at equilibrium, over ground states of the model. Samples are guaranteed to match the equilibrium for long 'smooth' beta schedules.

.. todo: example(s) specific to simulated annealing

>>> from dwave.samplers import SimulatedAnnealingSampler
...
>>> sampler = SimulatedAnnealingSampler()
...
>>> h = {0: -1, 1: -1}
>>> J = {(0, 1): -1}
>>> sampleset = sampler.sample_ising(h, J)

Steepest Descent
~~~~~~~~~~~~~~~~

Steepest descent is the discrete analogue of gradient descent, but the best move is computed using a local minimization rather rather than computing a gradient. At each step, we determine the dimension along which to descend based on the highest energy drop caused by a variable flip.

.. todo: example(s) specific to steepest descent

>>> from dwave.samplers import SteepestDescentSolver
...
>>> solver = SteepestDescentSolver()
>>> sampleset = solver.sample_ising({0: 2, 1: 1}, {(0, 1): -1}, initial_state={0: 1, 1: 1})
...
>>> sampleset.first.energy
-4.0

Tabu
~~~~

An implementation of the `MST2 multistart tabu search algorithm
<https://link.springer.com/article/10.1023/B:ANOR.0000039522.58036.68>`_
for quadratic unconstrained binary optimization (QUBO) problems.

.. todo: example(s) specific to tabu


>>> from dwave.samplers import TabuSampler
>>> sampleset = TabuSampler().sample_ising({'a': -0.5, 'b': 1.0}, {('a', 'b'): -1})

Tree Decomposition
~~~~~~~~~~~~~~~~~~

An implementation of a `generic tree decomposition-based solver <https://en.wikipedia.org/wiki/Tree_decomposition>`_  for binary quadratic models, using bucket tree elimination.

.. todo: example(s) specific to tree decomposition

>>> from dwave.samplers import TreeDecompositionSolver
...
>>> sampler = TreeDecompositionSolver()
>>> sampleset = sampler.sample_ising({0: .1}, {(0, 1): -1})
>>> sampleset.first.sample
{0: -1, 1: -1}

.. index-end-marker

Installation
------------

.. installation-start-marker

To install the core package:

.. code-block:: bash

    pip install dwave-samplers

.. installation-end-marker

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
