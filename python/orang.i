%module orang

%{
#define SWIG_FILE_WITH_INIT
#include <stdexcept>
#include "python-api.h"
#include <orang/exception.h>
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (std::bad_alloc&) {
    SWIG_exception(SWIG_MemoryError, "Out of memory");
  } catch (std::invalid_argument& e) {
    SWIG_exception(SWIG_ValueError, e.what());
  } catch (std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (orang::Exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what().c_str());
  } catch (...) {
    SWIG_exception(SWIG_UnknownError, "Unknown error");
  }
}

%include "numpy.i"

%init %{
import_array();
%}

%define SAMPLE_ISING_DOCSTRING
"Draw samples and compute marginals of an Ising problem.

log_pf, samples, single_mrg, pair_mrg, pairs = sample_ising(
    h, j, var_order, max_complexity, num_samples, marginals, beta, rng_seed)

    Args:
        h: linear Ising coefficients as a numpy array or sequence
        j: quadratic Ising coefficients as a numpy array or nested sequences
        var_order: variable elimination order as a numpy array or sequence
        max_complexity: sanity check.  Upper bound on exponent of space
            requirements at each step of the algorithm.
        num_samples: number of samples to draw
        marginals: whether or not to compute marginals (boolean).
            Marginals are always returned but will be empty if this
            parameter is False.
        beta: Boltzmann distribution inverse temperature parameter
        rng_seed: random number generator seed.  Negative values will cause
            a time-based seed to be used.
    Returns:
        (log_pf, samples, single_mrg, pair_mrg, pairs)
        log_pf: natural log of the partition function
        samples: numpy array of samples (each row is a sample)
        single_mrg: 1D numpy array of single-variable marginals.  Each
            value is the probability that the corresponding variable is 1.
        pair_mrg: numpy array of pairwise marginals.  Corresponding
            indices are in pairs return value.  Each row contains four
            values: [m_i-_j-, m_i+_j-, m_i-_j+, m_i+_j+] where m_ix_jy
            is the probability that variables i and j have respective
            signs x and y, where [i, j] is the corresponding row of the
            pairs return value.
        pairs: numpy array of pairwise marginal indices.  These are
            precisely the pairs with nonzero J coefficient."
%enddef

%define SAMPLE_QUBO_DOCSTRING
"Draw samples and compute marginals of a QUBO.

log_pf, samples, single_mrg, pair_mrg, pairs = sample_ising(
    q, var_order, max_complexity, num_samples, marginals, beta, rng_seed)

    Args:
        q: QUBO matrix as a numpy array or sequence
        var_order: variable elimination order as a numpy array or sequence
        max_complexity: sanity check.  Upper bound on exponent of space
            requirements at each step of the algorithm.
        num_samples: number of samples to draw
        marginals: whether or not to compute marginals (boolean).
            Marginals are always returned but will be empty if this
            parameter is False.
        beta: Boltzmann distribution inverse temperature parameter
        rng_seed: random number generator seed.  Negative values will cause
            a time-based seed to be used.
    Returns:
        (log_pf, samples, single_mrg, pair_mrg, pairs)
        log_pf: natural log of the partition function
        samples: numpy array of samples (each row is a sample)
        single_mrg: 1D numpy array of single-variable marginals.  Each
            value is the probability that the corresponding variable is 1.
        pair_mrg: numpy array of pairwise marginals.  Corresponding
            indices are in pairs return value.  Each row contains four
            values: [m_i0_j0, m_i1_j0, m_i0_j1, m_i1_j1] where m_ix_jy
            is the probability that variables i and j have respective
            values x and y, where [i, j] is the corresponding row of the
            pairs return value.
        pairs: numpy array of pairwise marginal indices.  These are
            precisely the pairs with nonzero Q coefficient."
%enddef

%define SOLVE_ISING_DOCSTRING
"Optimize an Ising problem.

log_pf, samples, single_mrg, pair_mrg, pairs = sample_ising(
    h, j, var_order, max_complexity, max_solutions)

    Args:
        h: linear Ising coefficients as a numpy array or sequence
        j: quadratic Ising coefficients as a numpy array or nested sequences
        var_order: variable elimination order as a numpy array or sequence
        max_complexity: sanity check.  Upper bound on exponent of space
            requirements at each step of the algorithm.
        max_solutions: maximum number of solutions to return
    Returns: (energies, solutions)
        energies: numpy array of energy values for corresponding
            solutions.  If max_solutions < 1, the optimum energy will
            still be returned.
        solutions: numpy array of solutions (each row is a solution)"
%enddef

%define SOLVE_QUBO_DOCSTRING
"Optimize a QUBO.

log_pf, samples, single_mrg, pair_mrg, pairs = sample_ising(
    q, var_order, max_complexity, max_solutions)

    Args:
        q: QUBO as a numpy array or sequence
        var_order: variable elimination order as a numpy array or sequence
        max_complexity: sanity check.  Upper bound on exponent of space
            requirements at each step of the algorithm.
        max_solutions: maximum number of solutions to return
    Returns: (energies, solutions)
        energies: numpy array of energy values for corresponding
            solutions.  If max_solutions < 1, the optimum energy will
            still be returned.
        solutions: numpy array of solutions (each row is a solution)"
%enddef

%feature("autodoc", SAMPLE_ISING_DOCSTRING) sample_ising;
%feature("autodoc", SAMPLE_QUBO_DOCSTRING) sample_qubo;
%feature("autodoc", SOLVE_ISING_DOCSTRING) solve_ising;
%feature("autodoc", SOLVE_QUBO_DOCSTRING) solve_qubo;


%numpy_typemaps(double, NPY_DOUBLE, int)
%numpy_typemaps(int, NPY_INT, int)

%apply double* OUTPUT {double* logPf}
%apply (double* IN_ARRAY1, int DIM1) {(double* hData, int hLen)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* jData, int jRows, int jCols)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* qData, int qRows, int qCols)}
%apply (int* IN_ARRAY1, int DIM1) {(int* voData, int voLen)}
%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** energiesData, int* energiesLen)}
%apply (int** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(int** solsData, int* solsRows, int* solsCols)};
%apply (int** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(int** samplesData, int* samplesRows, int* samplesCols)};
%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** singleMrgData, int* singleMrgLen)}
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** pairMrgData, int* pairMrgRows, int* pairMrgCols)};
%apply (int** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(int** pairData, int* pairRows, int* pairCols)};
%include "python-api.h"
