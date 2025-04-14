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

#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_sxxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_sxxx,
                                            const size_t idx_sxxx,
                                            const size_t idx_dxxx,
                                            const int b_angmom,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});
    
    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});
    
    for (int i = 0; i < bcomps; i++)
    {
        for (int j = 0; j < ccomps; j++)
        {
            for (int k = 0; k < dcomps; k++)
            {
                /// Set up components of auxilary buffer : SSSS

                const auto sx_off = idx_sxxx + i * ccomps * dcomps  + j * dcomps + k;

                auto g_0_0 = cbuffer.data(sx_off);
                
                /// Set up components of auxilary buffer : DSSS

                const auto dx_off = idx_dxxx + i * ccomps * dcomps  + j * dcomps + k;

                auto g_xx_0 = cbuffer.data(dx_off);
                
                auto g_xy_0 = cbuffer.data(dx_off + bcomps * ccomps * dcomps);
                
                auto g_xz_0 = cbuffer.data(dx_off + 2 * bcomps * ccomps * dcomps);
                
                auto g_yy_0 = cbuffer.data(dx_off + 3 * bcomps * ccomps * dcomps);
                
                auto g_yz_0 = cbuffer.data(dx_off + 4 * bcomps * ccomps * dcomps);
                
                auto g_zz_0 = cbuffer.data(dx_off + 5 * bcomps * ccomps * dcomps);
                
                /// set up bra offset for contr_buffer_ssxx

                const auto ss_geom_20_off = idx_geom_20_sxxx + i * ccomps * dcomps  + j * dcomps + k;
                
                auto g_xx_0_0_0 = cbuffer.data(ss_geom_20_off);
                
                auto g_xy_0_0_0 = cbuffer.data(ss_geom_20_off + bcomps * ccomps * dcomps);
                
                auto g_xz_0_0_0 = cbuffer.data(ss_geom_20_off + 2 * bcomps * ccomps * dcomps);
                
                auto g_yy_0_0_0 = cbuffer.data(ss_geom_20_off + 3 * bcomps * ccomps * dcomps);
                
                auto g_yz_0_0_0 = cbuffer.data(ss_geom_20_off + 4 * bcomps * ccomps * dcomps);
                
                auto g_zz_0_0_0 = cbuffer.data(ss_geom_20_off + 5 * bcomps * ccomps * dcomps);
                
                #pragma omp simd aligned(g_0_0, g_xx_0_0_0, g_xy_0_0_0, g_xz_0_0_0, g_yy_0_0_0, g_yz_0_0_0, g_zz_0_0_0, g_xx_0, g_xy_0, g_xz_0, g_yy_0, g_yz_0, g_zz_0 : 64)
                for (size_t l = 0; l < nelems; l++)
                {
                    const auto fact = g_0_0[l];
                    
                    g_xx_0_0_0[l] = g_xx_0[l] - fact;
                    
                    g_xy_0_0_0[l] = g_xy_0[l];
                    
                    g_xz_0_0_0[l] = g_xz_0[l];
                    
                    g_yy_0_0_0[l] = g_yy_0[l] - fact;
                    
                    g_yz_0_0_0[l] = g_yz_0[l];
                    
                    g_zz_0_0_0[l] = g_zz_0[l] - fact;
                }
            }
        }
    }
}

} // erirec namespace
