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


#include "SimdKineticEnergyFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"

namespace simdkin {  // simdkin namespace

auto
one_center_kinetic_energy(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    // NOTE: the kernels of the kinetic energy are not written yet. The value is
    // not returned as zero, as a matrix of zeros is indistinguishable from a
    // computed one and would be taken for an answer.

    errors::assertMsgCritical(
        false, std::string("SimdKineticEnergyFunc.one_center_kinetic_energy: Kernels of the kinetic energy are not implemented"));

    return 0.0;
}

auto
compute_kinetic_energy(double                         *values,
                       const size_t                    nvalues,
                       const CBasisFunction           &bra,
                       const CBasisFunction           &ket,
                       const std::vector<CSimdMatrix> &harmonics,
                       const CSimdMatrix              &coordinates,
                       const double                    threshold) -> void
{
    // NOTE: the kernels of the kinetic energy are not written yet. The driver
    // around them is complete and is exercised by the sparsity pattern and the
    // dimensions of the blocks, so this is reached only when the integrals
    // themselves are asked for.

    errors::assertMsgCritical(
        false, std::string("SimdKineticEnergyFunc.compute_kinetic_energy: Kernels of the kinetic energy are not implemented"));
}

}  // namespace simdkin
