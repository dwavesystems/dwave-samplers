.. image:: https://circleci.com/gh/dwavesystems/dwave-greedy.svg?style=svg
    :target: https://circleci.com/gh/dwavesystems/dwave-greedy
    :alt: Linux/Mac build status


============
dwave-greedy
============

.. index-start-marker

An implementation of a steepest descent solver for binary quadratic models.

Steepest descent is the discrete analogue of gradient descent, but the best
move is computed using a local minimization rather rather than computing a
gradient. At each step, we determine the dimension along which to descend based
on the highest energy drop caused by a variable flip.

.. code-block:: python

    >>> import greedy
    ...
    >>> solver = greedy.SteepestDescentSolver()
    >>> sampleset = solver.sample_ising({0: 2, 1: 2}, {(0, 1): -1})
    ...
    >>> print(sampleset)
        0  1 energy num_oc.
    0 -1 -1   -5.0       1
    ['SPIN', 1 rows, 1 samples, 2 variables]

.. index-end-marker


Installation
============

.. installation-start-marker

Install from a package on PyPI:

.. code-block:: bash

    pip install dwave-greedy

or install from source:

.. code-block:: bash

    USE_CYTHON=1 pip install -e git+https://github.com/dwavesystems/dwave-greedy.git#egg=dwave-greeedy

Note: ``USE_CYTHON=1`` forces Cythonization and proper build from source. When
building from *PyPI package* source (which includes Cythonized files), this is
not necessary.

To build from source:

.. code-block:: bash

    pip install -r requirements.txt
    python setup.py build_ext --inplace
    python setup.py install

.. installation-end-marker


Examples
========

.. example-start-marker

Simple frustrated Ising triangle:

.. code-block:: python

    import dimod
    import greedy

    # Construct a simple problem
    bqm = dimod.BQM.from_qubo({'ab': 1, 'bc': 1, 'ca': 1})

    # Instantiate the sampler
    sampler = greedy.SteepestDescentSampler()

    # Solve the problem
    result = sampler.sample(bqm)

Large RAN1 sparse problem (requires NetworkX_ package):

.. code-block:: python

    import dimod
    import greedy
    import networkx

    # Generate random Erdős-Rényi sparse graph with 10% density
    graph = networkx.fast_gnp_random_graph(n=1000, p=0.1)

    # Generate RAN1 problem on the sparse graph
    bqm = dimod.generators.random.ran_r(r=1, graph=graph)

    # Instantiate the sampler
    sampler = greedy.SteepestDescentSampler()

    # Run steepest descent for 100 times, each time from a random state
    sampleset = sampler.sample(bqm, num_reads=100)

    # Print the best energy
    print(min(sampleset.record.energy))

.. example-end-marker


License
=======

Released under the Apache License 2.0. See `<LICENSE>`_ file.


.. _NetworkX: https://networkx.github.io/