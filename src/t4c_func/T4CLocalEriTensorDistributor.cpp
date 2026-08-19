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

#include "T4CLocalEriTensorDistributor.hpp"

#include <iomanip>
#include <iostream>

#include "TensorComponents.hpp"

namespace t4cfunc {  // t4cfunc namespace

auto
local_distribute_eri_tensor(CDense4DTensor*                  eri_tensor,
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
                            const std::pair<size_t, size_t>& ket_range) -> void
{
    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = static_cast<int>(i * adim + refp);

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = static_cast<int>(j * bdim + refq);

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = static_cast<int>(k * cdim + refr);

                        const auto s = static_cast<int>(l * ddim + refs);

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        // Note: updating 8-fold pqrs integrals without omp critical
                        // since only unique pqrs combinations are computed

                        eri_tensor->row(p, q, r)[s] = fval;
                        eri_tensor->row(p, q, s)[r] = fval;

                        eri_tensor->row(q, p, r)[s] = fval;
                        eri_tensor->row(q, p, s)[r] = fval;

                        eri_tensor->row(r, s, p)[q] = fval;
                        eri_tensor->row(r, s, q)[p] = fval;

                        eri_tensor->row(s, r, p)[q] = fval;
                        eri_tensor->row(s, r, q)[p] = fval;
                    }
                }
            }
        }
    }
}

}  // namespace t4cfunc
