//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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
auto get_mpi_comm(py::object py_comm) -> MPI_Comm*;

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
