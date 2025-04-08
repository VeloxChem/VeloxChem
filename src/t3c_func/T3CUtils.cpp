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

#include "T3CUtils.hpp"

#include <ranges>

#include "MathConst.hpp"
#include "CustomViews.hpp"
#include "TensorComponents.hpp"

#include <iostream>

namespace t3cfunc {  // t3cfunc namespace


auto
unique_indices(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<size_t>
{
    std::vector<size_t> indices;
    
    std::ranges::for_each(gto_blocks, [&] (const auto& gblock) {
        const size_t ngtos = gblock.number_of_basis_functions();
        const size_t comps = tensor::number_of_spherical_components(std::array<int, 1>({gblock.angular_momentum(), }));
        const auto orb_ids = gblock.orbital_indices();
        std::ranges::for_each(views::rectangular(ngtos, comps), [&](const auto& index) {
            indices.push_back(orb_ids[index.first + 1] + orb_ids[0] * index.second);
        });
    });
    
    return indices;
}

auto
mask_indices(const std::vector<CGtoBlock>& gto_blocks) -> std::map<size_t, size_t>
{
    std::map<size_t, size_t> mask;
    
    size_t loc_idx = 0;
    
    std::ranges::for_each(gto_blocks, [&] (const auto& gblock) {
        const size_t ngtos = gblock.number_of_basis_functions();
        const size_t comps = tensor::number_of_spherical_components(std::array<int, 1>({gblock.angular_momentum(), }));
        const auto orb_ids = gblock.orbital_indices();
        std::ranges::for_each(views::rectangular(comps, ngtos), [&](const auto& index) {
            mask.insert({orb_ids[index.second + 1] + orb_ids[0] * index.first, loc_idx});
            loc_idx++;
        });
    });
    
    return mask;
}

auto
comp_coordinates_w(CSimdArray<double>&   buffer,
                   const size_t          index_w,
                   const size_t          index_q,
                   const TPoint<double>& r_a,
                   const double          a_exp) -> void
{
    // Set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up Cartesian W coordinates

    auto w_x = buffer.data(index_w);

    auto w_y = buffer.data(index_w + 1);

    auto w_z = buffer.data(index_w + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute Cartesian W center coordinates

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(w_x, w_y, w_z, q_x, q_y, q_z, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double cd_exp = c_exps[i] + d_exps[i];

        const double fact = 1.0 / (a_exp + cd_exp);

        w_x[i] = fact * (a_x * a_exp + q_x[i] * cd_exp);

        w_y[i] = fact * (a_y * a_exp + q_y[i] * cd_exp);

        w_z[i] = fact * (a_z * a_exp + q_z[i] * cd_exp);
    }
}

auto
comp_distances_aq(CSimdArray<double>&   buffer,
                  const size_t          index_aq,
                  const size_t          index_q,
                  const TPoint<double>& r_a) -> void
{
    // set up R(AQ) distances

    auto aq_x = buffer.data(index_aq);

    auto aq_y = buffer.data(index_aq + 1);

    auto aq_z = buffer.data(index_aq + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian P coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(PQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(aq_x, aq_y, aq_z, q_x, q_y, q_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        aq_x[i] = a_x - q_x[i];

        aq_y[i] = a_y - q_y[i];

        aq_z[i] = a_z - q_z[i];
    }
}

auto
comp_distances_wa(CSimdArray<double>&   buffer,
                  const size_t          index_wa,
                  const size_t          index_w,
                  const TPoint<double>& r_a) -> void
{
    // set up R(WA) distances

    auto wa_x = buffer.data(index_wa);

    auto wa_y = buffer.data(index_wa + 1);

    auto wa_z = buffer.data(index_wa + 2);

    // set up Cartesian W coordinates

    auto w_x = buffer.data(index_w);

    auto w_y = buffer.data(index_w + 1);

    auto w_z = buffer.data(index_w + 2);

    // set up Cartesian P coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(WQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(wa_x, wa_y, wa_z, w_x, w_y, w_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        wa_x[i] = w_x[i] - a_x;

        wa_y[i] = w_y[i] - a_y;

        wa_z[i] = w_z[i] - a_z;
    }
}

auto
comp_boys_args(CSimdArray<double>&       bf_data,
               const size_t              index_args,
               const CSimdArray<double>& buffer,
               const size_t              index_aq,
               const double              a_exp) -> void
{
    // Set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up R(PQ) distances

    auto aq_x = buffer.data(index_aq);

    auto aq_y = buffer.data(index_aq + 1);

    auto aq_z = buffer.data(index_aq + 2);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, aq_x, aq_y, aq_z, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double cd_exp = c_exps[i] + d_exps[i];

        bargs[i] = a_exp * cd_exp * (aq_x[i] * aq_x[i] + aq_y[i] * aq_y[i] + aq_z[i] * aq_z[i]) / (a_exp + cd_exp);
    }
}

auto
comp_ovl_factors(CSimdArray<double>& buffer,
                 const size_t        index_ovl,
                 const size_t        index_ket_ovl,
                 const size_t        index_ket_norm,
                 const double        a_norm,
                 const double        a_exp) -> void
{
    // set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up combined overlap

    auto fss = buffer.data(index_ovl);

    // set up ket data

    auto ket_ovls = buffer.data(index_ket_ovl);

    auto ket_norms = buffer.data(index_ket_norm);

    // set up pi constant

    const auto fpi = mathconst::pi_value();

    // compute combined overlap factors

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(fss, ket_ovls, ket_norms, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double cd_exp = c_exps[i] + d_exps[i];

        fss[i] = 2.0 * fpi * a_norm * ket_norms[i] * ket_ovls[i] * std::sqrt(cd_exp / (a_exp + cd_exp)) / a_exp;
    }
}


}  // namespace t3cfunc
