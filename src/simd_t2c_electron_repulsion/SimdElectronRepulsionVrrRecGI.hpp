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


#ifndef SimdElectronRepulsionVrrRecGI_hpp
#define SimdElectronRepulsionVrrRecGI_hpp

#include <cstddef>
#include "SimdMatrix.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t di0, const size_t di1,
                                     const size_t fh, const size_t fi, const size_t gg0,
                                     const size_t gg1, const size_t gh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_28(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_29(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_30(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_31(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_32(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_33(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_34(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_35(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_36(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_37(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_38(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_39(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_40(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_41(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_42(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_43(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_44(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_45(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_46(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_47(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

/// @brief Takes one step of the recurrence for one pair of primitives.
auto
compute_prim_gi_electron_repulsion_48(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t di0, const size_t di1,
                                      const size_t fh, const size_t fi, const size_t gg0,
                                      const size_t gg1, const size_t gh, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void;

}  // namespace simdt2ceri

#endif /* SimdElectronRepulsionVrrRecGI_hpp */
