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

#ifndef T4CDebug_hpp
#define T4CDebug_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "SimdArray.hpp"
#include "GtoPairBlock.hpp"

namespace t4cfunc {  // t4cfunc namespace

/// @brief Dumps buffer to output stream.
/// @param buffer The buffer to be dumped.
/// @param position The starting position in buffer to be dumped. 
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
/// @param bra_index The index of basis function pairs on bra side.
/// @param components The number of geometrical components.
auto
dump_buffer(const CSimdArray<double>&        buffer,
            const size_t                     position,
            const CGtoPairBlock&             bra_gto_pair_block,
            const CGtoPairBlock&             ket_gto_pair_block,
            const std::pair<size_t, size_t>& ket_indices,
            const size_t                     bra_index,
            const size_t                     components) -> void;

/// @brief Dumps buffer to output stream.
/// @param buffer The buffer to be dumped.
/// @param position The starting position in buffer to be dumped.
/// @param label The label of buffer.
/// @param elements The number of elements in row of buffer.
/// @param components The number of components to dump.
auto
dump_buffer(const CSimdArray<double>& buffer,
            const size_t              position,
            const std::string&        label,
            const size_t              elements,
            const size_t              components) -> void;

}  // namespace t4cfunc

#endif /* T4CDebug_hpp */
