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

#include "T3CLocalReducedDistributor.hpp"

#include "TensorComponents.hpp"
#include "MathFunc.hpp"

CT3CLocalReducedDistributor::CT3CLocalReducedDistributor(CT3FlatBuffer<double>* values, const std::vector<size_t>& origins)
{
    _t3_values = values;
    
    _origins = origins;
}

auto
CT3CLocalReducedDistributor::distribute(const CSimdArray<double>&        buffer,
                                        const size_t                     offset,
                                        const std::vector<size_t>&       a_indices,
                                        const std::vector<size_t>&       c_indices,
                                        const std::vector<size_t>&       d_indices,
                                        const int                        a_angmom,
                                        const int                        c_angmom,
                                        const int                        d_angmom,
                                        const size_t                     ibra_gto,
                                        const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up masked indices
    
    const auto mask_indices = _t3_values->mask_indices();
    
    // set up size of buffer
    
    const auto nrows = _t3_values->width();
    
    // reference indexes on bra side

    const auto refp = a_indices[ibra_gto + 1];
    
    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];
    
    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});
    
    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});
    
    for (int i = 0; i < acomps; i++)
    {
        auto ptr_values = _t3_values->data(mask_indices.at(i * adim + refp));
        
        for (auto j = ket_range.first; j < ket_range.second; j++)
        {
            // reference indexes on ket side

            const auto refr = c_indices[j + 1];

            const auto refs = d_indices[j + 1];
            
            // loop over ket components
            
            size_t index = 0;
            
            for (int k = 0; k < ccomps; k++)
            {
                const auto lstart = (refr == refs) ? k : 0;
                
                for (int l = lstart; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * ccomps * dcomps + k * dcomps + l);
                    
                    ptr_values[_origins[j] + index] = curr_buffer[j - ket_range.first];
                    
                    index++;
                }
            }
        }
    }
}
