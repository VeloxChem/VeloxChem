//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef TensorLabels_hpp
#define TensorLabels_hpp

#include <string>
#include <vector>

namespace tensor {  // tensor

/// @brief Generates all Cartesian labels for tensor.
/// @param order The order of tensor.
/// @return The vector of Cartesian tensor labels.
auto cartesian_labels(const int order) -> std::vector<std::string>;

/// @brief Generates all spherical labels for tensor.
/// @param order The order of tensor.
/// @return The vector of spherical tensor labels.
auto spherical_labels(const int order) -> std::vector<std::string>;

/// @brief Gets Cartesian index of canonical tensor component.
/// @param label The label of Cartesian tensor component.
/// @return The index of Cartesian tensor component.
auto cartesian_index(const std::string& label) -> int;

/// @brief Gets label of canonical tensor.
/// @param order The order of canonical tensor.
/// @return The label of canonical tensor.
auto label(const int order) -> char;

/// @brief Gets order of canonical tensor.
/// @param label The uppercased label of canonical tensor.
/// @return The order of canonical tensor.
auto order(const char label) -> int;

}  // namespace tensor

#endif /* TensorLabels_hpp */
