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

#ifndef T4CLocalEriTensorDistributor_hpp
#define T4CLocalEriTensorDistributor_hpp

#include <array>
#include <cstddef>

#include "Dense4DTensor.hpp"
#include "Matrices.hpp"
#include "Matrix.hpp"
#include "SimdArray.hpp"

namespace t4cfunc {  // t4cfunc namespace

auto local_distribute_eri_tensor(CDense4DTensor*                  eri_tensor,
                                 const CSimdArray<double>&        buffer,
                                 const size_t                     offset,
                                 const std::vector<size_t>&       a_indices,
                                 const std::vector<size_t>&       b_indices,
                                 const std::vector<size_t>&       c_indices,
                                 const std::vector<size_t>&       d_indices,
                                 const std::vector<size_t>&       a_loc_indices,
                                 const std::vector<size_t>&       b_loc_indices,
                                 const std::vector<size_t>&       c_loc_indices,
                                 const std::vector<size_t>&       d_loc_indices,
                                 const int                        a_angmom,
                                 const int                        b_angmom,
                                 const int                        c_angmom,
                                 const int                        d_angmom,
                                 const size_t                     bra_igto,
                                 const std::pair<size_t, size_t>& ket_range) -> void;

}  // namespace t4cfunc

#endif /* T4CLocalEriTensorDistributor_hpp */
