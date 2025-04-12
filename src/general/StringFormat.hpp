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
