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

#include "T4CDebug.hpp"

#include <iostream>
#include <iomanip>

namespace t4cfunc {  // t4cfunc namespace

auto
dump_buffer(const CSimdArray<double>&        buffer,
            const size_t                     position,
            const CGtoPairBlock&             bra_gto_pair_block,
            const CGtoPairBlock&             ket_gto_pair_block,
            const std::pair<size_t, size_t>& ket_indices,
            const size_t                     bra_index,
            const size_t                     components) -> void
{
    const auto bramom = bra_gto_pair_block.angular_momentums();
    
    const auto ketmom = ket_gto_pair_block.angular_momentums();
    
    const auto acomps = 2 * bramom.first + 1;
    
    const auto bcomps = 2 * bramom.second + 1;
    
    const auto ccomps = 2 * ketmom.first + 1;
    
    const auto dcomps = 2 * ketmom.second + 1;
    
    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();
    
    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();
    
    //std::cout << " * Angular Momentum: " << bramom.first << "," << bramom.second << "," << ketmom.first << "," << ketmom.second << std::endl;
    
    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    const auto bkoff = i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l;
                    
                    for (size_t m = 0; m < components; m++)
                    {
                        auto tints = buffer.data(position + m * acomps * bcomps * ccomps * dcomps + bkoff);
                     
                        for (auto n = ket_indices.first; n < ket_indices.second; n++)
                        {
                            std::cout << m << " ";
                            
                            std::cout << a_indices[0] * i + a_indices[bra_index + 1] << " ";
                            
                            std::cout << b_indices[0] * j + b_indices[bra_index + 1] << " ";
                            
                            std::cout << c_indices[0] * k + c_indices[n + 1] << " ";
                            
                            std::cout << d_indices[0] * l + d_indices[n + 1] << " ";
                            
                            std::cout << std::setprecision(15);
                            
                            std::cout << tints[n - ket_indices.first] << std::endl;
                        }
                    }
                }
            }
        }
    }
}

auto
dump_buffer(const CSimdArray<double>& buffer,
            const size_t              position,
            const std::string&        label,
            const size_t              elements,
            const size_t              components) -> void
{
    std::cout << " *** BUFFER : " << label << std::endl;
    
    for (size_t i = 0; i < components; i++)
    {
        std::cout << " * " << i << " : ";
        
        const auto tints = buffer.data(position + i);
        
        for (size_t j = 0; j < elements; j++)
        {
            std::cout << std::setprecision(15);
            
            std::cout << tints[j] << " ";
        }
        
        std::cout << std::endl;
    }
}

}  // namespace t4cfunc
