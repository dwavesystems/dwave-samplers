.. image:: https://img.shields.io/pypi/v/dwave-orang.svg
    :target: https://pypi.python.org/pypi/dwave-orang

.. image:: https://ci.appveyor.com/api/projects/status/2oc8vrxxh15ecgo1/branch/master?svg=true
    :target: https://ci.appveyor.com/project/dwave-adtt/dwave-orang

.. image:: https://codecov.io/gh/dwavesystems/dwave-orang/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/dwavesystems/dwave-orang

.. image:: https://readthedocs.com/projects/d-wave-systems-orang/badge/?version=latest
  :target: https://docs.ocean.dwavesys.com/projects/dwave-orang/en/latest/?badge=latest

.. image:: https://circleci.com/gh/dwavesystems/dwave-orang.svg?style=svg
    :target: https://circleci.com/gh/dwavesystems/dwave-orang

orang
=====

.. index-start-marker1

Generic tree decomposition-based solver.

.. index-end-marker1

Example Usage
-------------

.. index-start-marker2


>>> import orang
>>> sampler = orang.OrangSolver()
>>> samples = sampler.sample_ising({}, {(0, 1): -1})

Using the polynomial interface

>>> import orang
>>> poly = {('a'): -1, ('a', 'b', 'c'): 1}
>>> samples = orang.solve_polynomial(poly, 'SPIN', num_reads=3)

.. index-end-marker2

See the documentation for more examples.

Installation
------------

.. installation-start-marker

Boost is required to build.  On Ubuntu, you can install Boost like this:

.. code-block:: bash
    
    sudo apt-get install libboost-dev

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

Released under the Apache License 2.0. See LICENSE file.

Contribution
------------

See CONTRIBUTING.rst file.
