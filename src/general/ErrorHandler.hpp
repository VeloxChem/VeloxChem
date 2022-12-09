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

#ifndef ErrorHandler_hpp
#define ErrorHandler_hpp

#include <string>

namespace errors {  // errors namespace

namespace detail {
/** Print message with header to `std::cerr`.
 *
 * @param[in] message message to print.
 * @param[in] header title of message.
 *
 * A call like `msg("hello", "Log");` results in the following printout:
 *
 *    **** Log (process 1) ****
 *         hello
 */
auto msg(const std::string& message, const std::string& header) -> void;
}  // namespace detail

/** Prints message to `std::cerr` and aborts.
 *
 * @param[in] message what to print.
 */
auto msgCritical(const std::string& message) -> void;

/** Prints message and aborts in case of a critical error.
 *
 * @param[in] condition what to check.
 * @param[in] message what to print.
 */
auto assertMsgCritical(const bool condition, const std::string& message) -> void;

/** Prints warning message to `std::cerr`.
 *
 * @param[in] message what to print.
 */
auto msgWarning(const std::string& message) -> void;

/** Prints warning message if `condition` is not fulfilled.
 *
 * @param[in] condition what to check.
 * @param[in] message what to print.
 */
auto assertMsgWarning(const bool condition, const std::string& message) -> void;
}  // namespace errors

#endif /* ErrorHandler_hpp */
