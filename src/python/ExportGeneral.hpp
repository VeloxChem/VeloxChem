//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <mpi.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_general {  // vlx_general namespace

/**
 Gets MPI_Comm pointer from a mpi4py communicator object.

 @param py_comm mpi4py communicator object.
 @return the MPI_Comm pointer.
 */
MPI_Comm get_mpi_comm(py::object py_comm);

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
 Bind overloaded functions in a less verbose fashion

 Use as:

    vlx_general::overload_cast_<argument-list>()(name-of-function)

 NOTE This is a workaround for C++11 usage of pybind11::overload_cast.
 See: https://pybind11.readthedocs.io/en/stable/classes.html#overloaded-methods
 */
template <typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

/**
 Exports classes/functions in src/general to python.
 */
void export_general(py::module& m);

}  // namespace vlx_general

#endif /* ExportGeneral_hpp */
