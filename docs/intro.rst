============
Introduction
============

*Samplers* are processes that sample from low energy states of a problem’s objective function.
A binary quadratic model (BQM) sampler samples from low energy states in models such as those
defined by an Ising equation or a Quadratic Unconstrained Binary Optimization (QUBO) problem
and returns an iterable of samples, in order of increasing energy. A dimod :term:`sampler` provides
‘sample_qubo’ and ‘sample_ising’ methods as well as the generic BQM sampler method.

The :class:`~orag.OrangSampler` and :class:`~orag.OrangSolver` implement a low-treewidth exact
solving algorithm.

Example
=======

.. include:: ../README.rst
  :start-after: example-start-marker
  :end-before: example-end-marker
