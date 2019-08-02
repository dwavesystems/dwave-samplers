======================
D-Wave Greedy Samplers
======================

.. index-start-marker

A collection of simple BQM samplers, including:

*  Steepest descent Ising solver
*  Random greedy descent Ising solver
*  Exact (combinatorial) Ising solver

.. index-end-marker


Installation or Building
========================

.. installation-start-marker

Install from a package on PyPI::

    pip install dwave-greedy

or from source::

    pip install -e git+https://github.com/dwavesystems/dwave-greedy.git#egg=dwave-greeedy

.. installation-end-marker


Example
=======

.. example-start-marker

Simple frustrated Ising triangle:

.. code-block:: python

    import dimod
    import greedy

    # Construct a simple problem
    bqm = dimod.BQM.from_ising({}, {'ab': 1, 'bc': 1, 'ca': 1})

    # Instantiate the sampler
    sampler = greedy.SteepestDescentSampler()

    # Solve the problem
    result = sampler.sample(bqm)

Large RAN1 sparse problem (requires `networkx` package):

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
