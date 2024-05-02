//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>

#include "DenseMatrix.hpp"
#include "ExportHelpers.hpp"

namespace vlx_math {  // vlx_math namespace

/**
 Converts numpy array to CDenseMatrix

 @param arr the numpy array.
 @return a CDenseMatrix object.
 */
auto CDenseMatrix_from_numpy(const py::array_t<double>& arr) -> std::shared_ptr<CDenseMatrix>;

/**
 Exports classes/functions in src/math to python.
 */
auto export_math(py::module& m) -> void;

}  // namespace vlx_math

#endif /* ExportMath_hpp */
