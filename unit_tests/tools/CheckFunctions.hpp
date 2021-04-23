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

#ifndef CheckFunctions_hpp
#define CheckFunctions_hpp

#include <cstdint>
#include <vector>

namespace vlxtest {

void compare(const std::vector<double>& aVector, const double* bVector);

void compare(const std::vector<double>& aVector, const double* bVector, const double threshod);

void compare(const std::vector<int32_t>& aVector, const int32_t* bVector);

void compare(const int32_t* aVector, const int32_t* bVector, const int32_t nElements);

void compare(const double* aVector, const double* bVector, const int32_t nElements);

void compare(const std::vector<double>& aVector, const std::vector<double>& bVector);

void compare(const std::vector<int32_t>& aVector, const std::vector<int32_t>& bVector);

void checkNorm(const double* aVector, const int32_t nElements);
}  // namespace vlxtest

#endif /* CheckFunctions_hpp */
