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

#ifndef T2CTransform_hpp
#define T2CTransform_hpp

#include <algorithm>
#include <array>
#include <ranges>

#include "CustomViews.hpp"
#include "SphericalMomentum.hpp"
#include "TensorComponents.hpp"

namespace t2cfunc {  // t2cfunc namespace

/// @brief Transforms Cartesian integrals buffer to spherical integrals buffer.
/// @tparam N The order of angular momentum tensor on bra side.
/// @tparam M The order of angular momentum tensor on ket side.
/// @param sbuffer The spherical integrals buffer.
/// @param cbuffer The Cartesian integrals array.
template <int N, int M>
inline auto
transform(CSimdArray<double>& sbuffer, const CSimdArray<double>& cbuffer) -> void
{
    const auto ndims = sbuffer.number_of_active_elements();

    const auto bra_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{N});

    const auto ket_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{M});

    const auto bra_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{N});

    const auto ket_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{M});

    const auto nblocks = sbuffer.number_of_rows() / (bra_spher_comps * ket_spher_comps);

    std::ranges::for_each(std::views::iota(size_t{0}, nblocks), [&](const auto n) {
        const auto cart_off  = n * bra_cart_comps * ket_cart_comps;
        const auto spher_off = n * bra_spher_comps * ket_spher_comps;
        std::ranges::for_each(views::rectangular(bra_spher_comps, ket_spher_comps), [&](const auto& index) {
            const auto [i, j] = index;
            auto dst_ptr      = sbuffer.data(spher_off + i * ket_spher_comps + j);
            for (const auto& [bra_idx, bra_fact] : spher_mom::transformation_factors<N>(i))
            {
                for (const auto& [ket_idx, ket_fact] : spher_mom::transformation_factors<M>(j))
                {
                    auto         src_ptr = cbuffer.data(cart_off + bra_idx * ket_cart_comps + ket_idx);
                    const double fact    = bra_fact * ket_fact;
#pragma omp simd aligned(dst_ptr, src_ptr : 64)
                    for (size_t k = 0; k < ndims; k++)
                    {
                        dst_ptr[k] += fact * src_ptr[k];
                    }
                }
            }
        });
    });
}

}  // namespace t2cfunc

#endif /* T2CTransform_hpp */
