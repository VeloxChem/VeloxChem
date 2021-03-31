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

#ifndef ExportOneInts_hpp
#define ExportOneInts_hpp

#include <pybind11/pybind11.h>

#include "CartesianComponents.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"

namespace py = pybind11;

namespace vlx_oneints {  // vlx_oneints namespace
/** Convert a VeloxChem matrix object, i.e. COverlapMatrix, to NumPy array.
 *
 * @tparam T type of the VeloxChem matrix object.
 * @param obj VeloxChem matrix object.
 * @return The NumPy array.
 */
template <typename T>
inline py::array_t<double>
matrix_to_numpy(const T& obj)
{
    return vlx_general::pointer_to_numpy(obj.values(), obj.getNumberOfRows(), obj.getNumberOfColumns());
}

/** Convert a VeloxChem matrix object with Cartesian components, i.e. CElectricDipoleMomentMatrix, to NumPy array.
 *
 * @tparam T type of the VeloxChem matrix object.
 * @tparam cart requested Cartesian component.
 * @param obj VeloxChem matrix object.
 * @return The NumPy array.
 */
template <typename T, cartesians cart>
inline py::array_t<double>
matrix_to_numpy(const T& obj)
{
    return vlx_general::pointer_to_numpy(obj.values(cart), obj.getNumberOfRows(), obj.getNumberOfColumns());
}

/** Convert a VeloxChem matrix object with Cartesian components, i.e. CElectricDipoleMomentMatrix, to NumPy array.
 *
 * @tparam T type of the VeloxChem matrix object.
 * @param obj VeloxChem matrix object.
 * @param cart requested Cartesian component.
 * @return The NumPy array.
 */
template <typename T>
inline py::array_t<double>
matrix_to_numpy(const T& obj, cartesians cart)
{
    return vlx_general::pointer_to_numpy(obj.values(cart), obj.getNumberOfRows(), obj.getNumberOfColumns());
}

/** Convert a NumPy array to a VeloxChem matrix object, i.e. COverlapMatrix.
 *
 * @tparam T type of the VeloxChem matrix object.
 * @param np the NumPy array.
 * @return The VeloxChem matrix object.
 */
template <typename T>
inline std::shared_ptr<T>
matrix_from_numpy(const py::array_t<double>& np)
{
    auto mp = vlx_math::CDenseMatrix_from_numpy(np);
    return std::make_shared<T>(*mp);
}

/**
 Exports classes/functions in src/oneints to python.
 */
void export_oneints(py::module& m);
}  // namespace vlx_oneints

#endif /* ExportOneInts_hpp */
