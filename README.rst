.. image:: https://img.shields.io/pypi/v/dwave-orang.svg
    :target: https://pypi.python.org/pypi/dwave-orang

.. image:: https://codecov.io/gh/dwavesystems/dwave-orang/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/dwavesystems/dwave-orang

.. image:: https://readthedocs.com/projects/d-wave-systems-orang/badge/?version=latest
  :target: https://docs.ocean.dwavesys.com/projects/dwave-orang/en/latest/?badge=latest

.. image:: https://circleci.com/gh/dwavesystems/dwave-orang.svg?style=svg
    :target: https://circleci.com/gh/dwavesystems/dwave-orang

===========
dwave-orang
===========

.. index-start-marker

An implementation of a `generic tree decomposition-based solver <https://en.wikipedia.org/wiki/Tree_decomposition>`_ 
for binary quadratic models, using bucket tree elimination.

This package contains both an exact solver (OrangSolver) and a Boltzmann sampler (OrangSampler).

.. index-end-marker

Example Usage
-------------

.. example-start-marker

>>> import orang
>>> sampler = orang.OrangSolver()
>>> sampleset = sampler.sample_ising({0: .1}, {(0, 1): -1})
>>> sampleset.first.sample
{0: -1, 1: -1}

.. example-end-marker

See the documentation for more examples.

Installation
------------

.. installation-start-marker

To build from source:

.. code-block:: bash

    pip install -r requirements.txt
    python setup.py build_ext --inplace

To install from source:

.. code-block:: bash

    pip install -r requirements.txt
    python setup.py build_ext --inplace
    pip install .


.. installation-end-marker

License
-------

Released under the Apache License 2.0. See `<LICENSE>`_ file.

Contributing
------------

Ocean's `contributing guide <https://docs.ocean.dwavesys.com/en/stable/contributing.html>`_
has guidelines for contributing to Ocean packages.
