//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <mpi.h>

namespace py = pybind11;

namespace vlx_general {  // vlx_general namespace

/**
 Gets MPI_Comm pointer from a mpi4py communicator object.

 @param py_comm mpi4py communicator object.
 @return the MPI_Comm pointer.
 */
MPI_Comm* get_mpi_comm(py::object py_comm);

/**
 Gets numpy array from double pointer and int32_t dimension.

 @param ptr the double pointer.
 @param dimension the shape of numpy array.
 @return numpy array.
 */
py::array_t<double> pointer_to_numpy(const double* ptr, const std::vector<int32_t>& dimension);

/**
 Gets 1d numpy array from double pointer and int32_t dimension.

 @param ptr the double pointer.
 @param nElements number of elements.
 @return numpy array.
 */
py::array_t<double> pointer_to_numpy(const double* ptr, const int32_t nElements);

/**
 Gets 2d numpy array from double pointer and int32_t dimension.

 @param ptr the double pointer.
 @param nRows number of rows.
 @param nColumns number of columns.
 @return numpy array.
 */
py::array_t<double> pointer_to_numpy(const double* ptr, const int32_t nRows, const int32_t nColumns);

/**
 Gets numpy array from int32_t pointer and dimension.

 @param ptr the int32_t pointer.
 @param dimension the shape of numpy array.
 @return numpy array.
 */
py::array_t<int32_t> pointer_to_numpy(const int32_t* ptr, const std::vector<int32_t>& dimension);

/**
 Gets 1d numpy array from int32_t pointer and dimension.

 @param ptr the int32_t pointer.
 @param nElements number of elements.
 @return numpy array.
 */
py::array_t<int32_t> pointer_to_numpy(const int32_t* ptr, const int32_t nElements);

/**
 Gets 2d numpy array from int32_t pointer and dimension.

 @param ptr the int32_t pointer.
 @param nRows number of rows.
 @param nColumns number of columns.
 @return numpy array.
 */
py::array_t<int32_t> pointer_to_numpy(const int32_t* ptr, const int32_t nRows, const int32_t nColumns);

/**
 Exports classes/functions in src/general to python.
 */
void export_general(py::module& m);

}  // namespace vlx_general

#endif /* ExportGeneral_hpp */
