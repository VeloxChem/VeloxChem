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

#include "T3CGeom0X0Distributor.hpp"

#include "TensorComponents.hpp"

CT3CGeom0X0Distributor::CT3CGeom0X0Distributor(CT3RectFlatBuffer<double>* values)
{
    _t3_values = values;
}

auto
CT3CGeom0X0Distributor::distribute(const CSimdArray<double>&        buffer,
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
    
    // mask indices
    
    const auto mask_indices = _t3_values->mask_indices();
    
    // number of rows along z axis
    
    const auto nrows = _t3_values->width();
    
    // set up geom. deriv. blocks data
    
    const auto gblocks = _t3_values->aux_blocks();
    
    const auto grows = _t3_values->aux_width();
    
    for (size_t n = 0; n < gblocks; n++)
    {
        for (int i = 0; i < acomps; i++)
        {
            auto ptr_values = _t3_values->data(n * grows + i * adim + refp);
            
            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + n * acomps * ccomps * dcomps + i * ccomps * dcomps + k * dcomps + l);
                    
                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side
                    
                        const auto refr = c_indices[m + 1];
                    
                        const auto refs = d_indices[m + 1];
                    
                        const auto r = k * cdim + refr;
                    
                        const auto s = l * ddim + refs;
                    
                        // assign integrals
                    
                        ptr_values[mask_indices.at(r) * nrows + s] = curr_buffer[m - ket_range.first];
                       
                    }
                }
            }
        }
    }
}

