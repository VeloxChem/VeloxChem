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

#ifndef ElectronRepulsionPrimRecSDSH_hpp
#define ElectronRepulsionPrimRecSDSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SD|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sdsh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sssh The primitive integrals buffer.
/// @param idx_eri_1_sssh The primitive integrals buffer.
/// @param idx_eri_1_spsg The primitive integrals buffer.
/// @param idx_eri_0_spsh The primitive integrals buffer.
/// @param idx_eri_1_spsh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sdsh(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sdsh,
                                       size_t                idx_eri_0_sssh,
                                       size_t                idx_eri_1_sssh,
                                       size_t                idx_eri_1_spsg,
                                       size_t                idx_eri_0_spsh,
                                       size_t                idx_eri_1_spsh,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSDSH_hpp */
