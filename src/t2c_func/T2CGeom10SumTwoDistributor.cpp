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

#include "T2CGeom10SumTwoDistributor.hpp"

#include <iostream>

#include "TensorComponents.hpp"

CT2CGeom10SumTwoDistributor::CT2CGeom10SumTwoDistributor(double* values, const double* bra_gamma, const double* ket_gamma)
{
    _grad_values = values;
    
    _ptr_bra_gamma = bra_gamma;
    
    _ptr_ket_gamma = ket_gamma;
}

auto
CT2CGeom10SumTwoDistributor::distribute(const CSimdArray<double>&        buffer,
                                        const std::vector<size_t>&       bra_indices,
                                        const std::vector<size_t>&       ket_indices,
                                        const int                        bra_angmom,
                                        const int                        ket_angmom,
                                        const size_t                     bra_igto,
                                        const std::pair<size_t, size_t>& ket_range,
                                        const bool                       diagonal) -> void
{
    // reference indexes on bra side

    const auto refp = bra_indices[bra_igto + 1];
    
    // dimensions of bra and ket orbital indexes

    const auto adim = bra_indices[0];

    const auto bdim = ket_indices[0];
        
    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});
    
    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{ket_angmom});
    
    for (size_t n = 0; n < 3; n++)
    {
        for (int i = 0; i < acomps; i++)
        {
            const auto p = i * adim + refp;
            
            for (int j = 0; j < bcomps; j++)
            {
                auto curr_buffer = buffer.data(n * acomps * bcomps + i * bcomps + j);
                
                for (size_t k = ket_range.first; k < ket_range.second; k++)
                {
                    const auto q = bdim * j + ket_indices[k + 1];
                    
                    const auto fval = curr_buffer[k - ket_range.first];
                    
                    _grad_values[n] += _ptr_bra_gamma[p] * fval * _ptr_ket_gamma[q] + _ptr_ket_gamma[p] * fval * _ptr_bra_gamma[q];
                }
            }
        }
    }
}
