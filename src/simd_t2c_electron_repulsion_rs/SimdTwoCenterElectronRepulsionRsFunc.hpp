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




#ifndef SimdTwoCenterElectronRepulsionRsFunc_hpp
#define SimdTwoCenterElectronRepulsionRsFunc_hpp

#include <cstddef>
#include <utility>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

/// @brief Computes the Coulomb and the range separated Coulomb integrals of two
/// basis functions centered on the same atom.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param omega The range separation parameter.
/// @return The integral of 1 / r and the integral of erf(omega r) / r.
/// @note The plain integral is the closed formula the unattenuated dispatcher
/// carries, and the attenuated one is that formula with the pair of primitives
/// weighted by theta to the 2l + 1, where theta is omega over the root of omega
/// squared plus the reduced exponent of the pair. The weight sits inside the sum
/// over the primitives, as the reduced exponent differs between them, and it is
/// exact rather than fitted: checked against pyscf over the angular momenta zero to
/// four, three exponents and four values of omega, to 1.9e-11.
/// @note Basis functions of different angular momenta have no integral and give a
/// pair of zeros, as in the unattenuated case.
auto one_center_rs_electron_repulsion(const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const double          omega) -> std::pair<double, double>;

/// @brief Computes the Coulomb and the range separated Coulomb integrals of a
/// combination of basis functions over the atom pairs of a block.
/// @param values The values of the combination, as one row of nvalues columns for
/// each of the (2 l_bra + 1) (2 l_ket + 1) spherical components, the components of
/// the bra side running slowest. The block of the Coulomb integrals comes first and
/// the block of the range separated ones follows it, each of that length.
/// @param nvalues The number of atom pairs of the block.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs, as ten rows ordered by
/// ascending interatomic distance, holding the vector between the atoms in rows six
/// to eight and its squared length in row nine.
/// @param buffer The buffer of the combination of basis functions.
/// @param omega The range separation parameter.
/// @note The two operators are computed together and not one after the other. They
/// share the recurrence, the primitives and the transformation, and differ only in
/// the Boys function each of them carries, so a kernel which forms both costs far
/// less than two kernels which form one each.
/// @note There is no screening threshold, for the reason the unattenuated driver
/// records: the Coulomb operator reaches every atom pair of every molecule. The
/// attenuated operator decays faster and could be screened, but it is not computed
/// on its own here and the pair it would drop is a pair the other operator needs.
auto compute_rs_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates,
                                   CSimdMatrix          &buffer,
                                   const double          omega) -> void;

}  // namespace simdt2ceri

#endif /* SimdTwoCenterElectronRepulsionRsFunc_hpp */
