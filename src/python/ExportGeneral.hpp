#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <mpi.h>
// see here: https://github.com/mpi4py/mpi4py/issues/19#issuecomment-768143143
#ifdef MSMPI_VER
#define PyMPI_HAVE_MPI_Message 1
#endif
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ExportHelpers.hpp"

namespace py = pybind11;

namespace vlx_general {  // vlx_general namespace

/**
 Gets pointer to MPI communicator from mpi4py communicator object.

 @param py_comm the mpi4py communictor object.
 @return the pointer to MPI communicator.
 */
MPI_Comm* get_mpi_comm(py::object py_comm);

/**
 Gets numpy array from pointer and shape.

 @param ptr pointer to data.
 @param dimension the shape of numpy array.
 @return numpy array.
 */
auto pointer_to_numpy(const double* ptr, const std::vector<int64_t>& dimension) -> py::array_t<double>;

/**
 Exports classes/functions in src/general to python.
 */
auto export_general(py::module& m) -> void;

}  // namespace vlx_general

#endif /* ExportGeneral_hpp */
