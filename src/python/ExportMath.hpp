//
//                           VELOXCHEM 1.0-RC2
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

#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>

#include "DenseMatrix.hpp"
#include "Dense4DTensor.hpp"

namespace py = pybind11;

namespace vlx_math {  // vlx_math namespace

/**
 Converts numpy array to CDenseMatrix

 @param arr the numpy array.
 @return a CDenseMatrix object.
 */
std::shared_ptr<CDenseMatrix> CDenseMatrix_from_numpy(const py::array_t<double>& arr);

/**
 Converts numpy array to CDense4DTensor

 @param arr the numpy array.
 @return a CDense4DTensor object.
 */
std::shared_ptr<CDense4DTensor> CDense4DTensor_from_numpy(const py::array_t<double>& arr);

/**
 Exports classes/functions in src/math to python.
 */
void export_math(py::module& m);

}  // namespace vlx_math

#endif /* ExportMath_hpp */
