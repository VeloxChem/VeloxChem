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

#ifndef newints_ShellPairSchwarzDriver_hpp
#define newints_ShellPairSchwarzDriver_hpp

#include "CauchySchwarzPairData.hpp"

class CMolecule;
class CMolecularBasis;

namespace newints {

/// @brief Builder for the orbital shell-pair Cauchy-Schwarz factors
/// Q_ij = max sqrt((mu nu|mu nu)) (see CauchySchwarzPairData).
///
/// STUB: the real builder requires four-center (mu nu|mu nu) ERI-diagonal kernels,
/// which are a separate effort not yet available in newints. The signature and the
/// intended algorithm are fixed here so the integration point is visible; compute()
/// throws until the kernels are wired in.
class ShellPairSchwarzDriver
{
   public:
    /// @brief The default constructor.
    ShellPairSchwarzDriver() = default;

    /// @brief Builds the significant orbital shell-pair Cauchy-Schwarz factors.
    /// @param molecule The molecule.
    /// @param basis The molecular (orbital) basis.
    /// @param threshold The screening threshold; pairs with Q_ij < threshold are discarded.
    /// @return The significant shell pairs and their factors.
    auto compute(const CMolecule &molecule, const CMolecularBasis &basis, const double threshold) const
        -> CauchySchwarzPairData;
};

}  // namespace newints

#endif /* newints_ShellPairSchwarzDriver_hpp */
