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
