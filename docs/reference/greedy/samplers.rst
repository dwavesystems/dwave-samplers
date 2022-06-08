.. _included_samplers:

========
Samplers
========

The `dwave-greedy` package currently includes just one sampler,
:class:`~greedy.sampler.SteepestDescentSampler`, which is an alias for
:class:`~greedy.sampler.SteepestDescentSolver`.

A :term:`sampler` accepts a :term:`binary quadratic model` (BQM) and returns
variable assignments. Samplers generally try to find minimizing values but can
also sample from distributions defined by the BQM.

.. currentmodule:: greedy.sampler


SteepestDescentSolver
=====================

Class
-----

.. autoclass:: SteepestDescentSolver

Attributes
----------

.. autosummary::
   :toctree: generated/

   SteepestDescentSolver.properties
   SteepestDescentSolver.parameters

Methods
-------

.. autosummary::
   :toctree: generated/

   SteepestDescentSolver.sample
   SteepestDescentSolver.sample_ising
   SteepestDescentSolver.sample_qubo


SteepestDescentSampler
======================

Class
-----

.. autoclass:: SteepestDescentSampler
