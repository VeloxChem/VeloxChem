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

#include "ShellPairSchwarzDriver.hpp"

#include <stdexcept>

namespace newints {

auto
ShellPairSchwarzDriver::compute([[maybe_unused]] const CMolecule       &molecule,
                                [[maybe_unused]] const CMolecularBasis  &basis,
                                [[maybe_unused]] const double            threshold) const -> CauchySchwarzPairData
{
    // TODO: this builder requires four-center (mu nu|mu nu) ERI-diagonal kernels,
    // which are a separate effort not yet present in newints. Intended algorithm:
    //
    //   gather the contracted shells in atom-major order (as the other drivers do);
    //   for each unique shell pair (i >= j):
    //       evaluate the (mu nu|mu nu) diagonal block (mu in shell i, nu in shell j)
    //       Q_ij = max over (mu, nu) of sqrt((mu nu|mu nu))
    //       if Q_ij >= threshold: append ShellPair{i, j, Q_ij}
    //   return CauchySchwarzPairData(number_of_shells, threshold, significant_pairs);
    //
    // The work is naturally parallel over shell pairs and stores only significant
    // pairs, matching CauchySchwarzPairData's sparse layout.

    throw std::logic_error(
        "newints::ShellPairSchwarzDriver: (mu nu|mu nu) ERI-diagonal kernels are not yet implemented");
}

}  // namespace newints
