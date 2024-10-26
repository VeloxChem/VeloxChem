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

#ifndef ErrorHandler_hpp
#define ErrorHandler_hpp

#include <string>

namespace errors {  // errors namespace

/** Print message with header to `std::cerr`.
 *
 * @param[in] message message to print.
 * @param[in] header title of message.
 */
auto msg(const std::string& message, const std::string& header) -> void;

/** Prints message and aborts in case of a critical error.
 *
 * @param[in] condition what to check.
 * @param[in] message what to print.
 */
auto assertMsgCritical(const bool condition, const std::string& message) -> void;

}  // namespace errors

#endif /* ErrorHandler_hpp */
