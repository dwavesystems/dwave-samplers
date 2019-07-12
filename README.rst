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

.. code-block:: python

    import dimod
    import greedy

    # Construct a problem
    bqm = dimod.BQM.from_ising({}, {'ab': 1, 'bc': -1, 'ca': 1})

    # Instantiate the sampler
    solver = greedy.SteepestDescent()

    # Solve the problem
    result = solver.sample()

.. example-end-marker


License
=======

Released under the Apache License 2.0. See `<LICENSE>`_ file.
