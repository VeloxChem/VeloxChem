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

#ifndef StringFormat_hpp
#define StringFormat_hpp

#include <cstdint>
#include <string>
#include <vector>

#include "FmtType.hpp"

namespace fstr {  // fstr namespace

/**
 Creates uppercased string from string.

 @param source the string.
 @return the uppercased string.
 */
auto upcase(const std::string& source) -> std::string;

/**
 Creates formatted string with requested width from string.

 @param source the string.
 @param width the width of formatted string.
 @param aligment the alignment of formatted string.
 @return the formatted string.
 */
auto format(const std::string& source, const size_t width, const fmt_t aligment) -> std::string;

/**
 Creates formatted string with requested width from real number.

 @param source the real number.
 @param presicion the conversion precision in decimal places.
 @param width the width of formatted string.
 @param aligment the alignment of formatted string.
 @return the formatted string.
 */
auto to_string(const double source, const size_t presicion, const size_t width, const fmt_t aligment) -> std::string;

/**
 Creates string with requested precision from real number.

 @param source the real number.
 @param presicion the conversion precision in decimal places.
 @return the formatted string.
 */
auto to_string(const double source, const size_t presicion) -> std::string;

/**
 Creates formatted string with requested width from integer number.

 @param source the integer number.
 @param width the width of formatted string.
 @param aligment the alignment of formatted string.
 @return the formatted string.
 */
auto to_string(const int64_t source, const size_t width, const fmt_t aligment) -> std::string;

/**
 Creates formatted string from boolean.

 @param source the boolean.
 @return the formatted string.
 */
auto to_string(const bool source) -> std::string;

/**
 Converts angular momentum label to angular momentum quantum number.
 Supported angular momentum:  S - I.

 @param label the angular momentum label.
 @return the angular momentum quantum number.
 */
auto to_AngularMomentum(const std::string& label) -> int64_t;

/**
 Converts angular momentum quantum number to angular momentum label.
 Supported angular momentum: S - I.

 @param angmom the angular momentum quantum number.
 @return the angular momentum label.
 */
auto to_AngularMomentum(const int64_t angmom) -> std::string;

/**
 Converts tensor of given order to vector of it's component labels.

 @param torder the order of tensor.
 @return the vector of tensor component labels.
 */
auto to_TensorComponents(const int64_t torder) -> std::vector<std::string>;

/**
 Converts tensor component label to it's standard index.

 @param tlabel the tensor component label.
 @return the index of tensor component.
 */
auto to_TensorComponent(const std::string& tlabel) -> int64_t;

}  // namespace fstr

#endif /* StringFormat_hpp */
