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



#ifndef SimdKineticEnergyFunc_hpp
#define SimdKineticEnergyFunc_hpp

#include <cmath>
#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "SimdMatrix.hpp"

namespace simdkin {  // simdkin namespace

/// @brief Computes the kinetic energy integral of two basis functions centered on
/// the same atom.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @return The kinetic energy integral of the basis functions.
/// @note The kinetic energy operator is spherically symmetric about an atom, so the
/// integral of two solid harmonic primitives on the same center is diagonal in the
/// angular components and does not depend on the position of the atom. It is
/// (2 l + 3) a b / (a + b) times their overlap, which reduces to the familiar
/// 3 a b / (a + b) for two S type primitives. Basis functions of different angular
/// momenta do not overlap and give zero.
inline auto
one_center_kinetic_energy(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if (lbra != lket) return 0.0;

    const auto fdfact = mathfunc::double_factorial(2 * lbra - 1);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    auto fsum = 0.0;

    for (size_t i = 0; i < a_exps.size(); i++)
    {
        for (size_t j = 0; j < b_exps.size(); j++)
        {
            const auto fab = 1.0 / (a_exps[i] + b_exps[j]);

            auto fovl = mathconst::pi_value() * fab;

            fovl = a_norms[i] * b_norms[j] * fovl * std::sqrt(fovl);

            // NOTE: repeated multiplication is used, as the angular momentum is
            // small.

            const auto fhalf = 0.5 * fab;

            for (int k = 0; k < lbra; k++)
            {
                fovl *= fhalf;
            }

            fsum += fdfact * fovl * static_cast<double>(2 * lbra + 3) * a_exps[i] * b_exps[j] * fab;
        }
    }

    return fsum;
}

/// @brief Computes the kinetic energy integrals of a combination of basis functions
/// by dispatching to the kernel of their angular momenta.
/// @param values The values of the combination of basis functions in the values
/// block of the sparsity pattern.
/// @param nvalues The number of values to compute, i.e. the number of atom pairs
/// surviving the screening of the combination of basis functions.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs, as ten rows ordered by
/// ascending interatomic distance, holding the vector between the atoms in rows six
/// to eight and its squared length in row nine.
/// @param buffer The buffer of the block, which spans every combination the block
/// carries and is sized by simdfunc from the highest angular momenta it holds.
/// @param threshold The screening threshold of the integrals.
/// @note The values of the combination of basis functions are stored as one row of
/// nvalues columns for each of the (2 l_bra + 1) (2 l_ket + 1) spherical components,
/// with the components of the bra side running slowest.
/// @note The kernels are generated up to angular momentum six, so a combination
/// above it stops with an error rather than returning values it did not compute. The
/// kinetic energy of two basis functions on the same atom is unaffected and is
/// computed by one_center_kinetic_energy above.
auto compute_kinetic_energy(double               *values,
                            const size_t          nvalues,
                            const CBasisFunction &bra,
                            const CBasisFunction &ket,
                            const CSimdMatrix    &coordinates,
                            CSimdMatrix          &buffer,
                            const double          threshold) -> void;

}  // namespace simdkin

#endif /* SimdKineticEnergyFunc_hpp */
