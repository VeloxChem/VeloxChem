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


#ifndef SimdElectronRepulsionVrrRecPK_hpp
#define SimdElectronRepulsionVrrRecPK_hpp

#include <cstddef>
#include "SimdMatrix.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t si, const size_t sk,
                                     const size_t pi, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_pk_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t si, const size_t sk, const size_t ncols,
                                     const double p) -> void;

}  // namespace simdt2ceri

#endif /* SimdElectronRepulsionVrrRecPK_hpp */
