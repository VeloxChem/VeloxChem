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

#ifndef GtoPairBlockFunc_hpp
#define GtoPairBlockFunc_hpp

#include <vector>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtofunc {  // gtofunc namespace

/// @brief Creates vector of basis function pairs blocks.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @return The vector of basis function pairs blocks.
auto make_gto_pair_blocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoPairBlock>;

/// @brief Creates vector of basis function pairs blocks for vector of basis functions blocks.
/// @param gto_blocks The vector of basis functions blocks on bra and ket sides.
/// @return The vector of basis function pairs blocks.
auto make_gto_pair_blocks(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<CGtoPairBlock>;

/// @brief Creates vector of basis function pairs blocks for pair of vectors of basis functions blocks.
/// @param bra_gto_blocks The vector of basis functions blocks on bra side.
/// @param ket_gto_blocks The vector of basis functions blocks on ket side.
/// @return The vector of basis function pairs blocks.
auto make_gto_pair_blocks(const std::vector<CGtoBlock>& bra_gto_blocks, const std::vector<CGtoBlock>& ket_gto_blocks) -> std::vector<CGtoPairBlock>; 

}  // namespace gtofunc

#endif /* GtoPairBlockFunc_hpp */
