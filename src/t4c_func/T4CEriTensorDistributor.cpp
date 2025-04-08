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

#include "T4CEriTensorDistributor.hpp"

#include "StringFormat.hpp"
#include "T4CLocalDistributor.hpp"
#include "T4CLocalEriTensorDistributor.hpp"
#include "T4CUtils.hpp"
#include "TensorComponents.hpp"

CT4CEriTensorDistributor::CT4CEriTensorDistributor(CDense4DTensor* eri_tensor)

    : _eri_tensor(eri_tensor)

    , _a_loc_indices(std::vector<size_t>())

    , _b_loc_indices(std::vector<size_t>())

    , _c_loc_indices(std::vector<size_t>())

    , _d_loc_indices(std::vector<size_t>())

    , _a_glob_indices(std::vector<size_t>())

    , _b_glob_indices(std::vector<size_t>())

    , _c_glob_indices(std::vector<size_t>())

    , _d_glob_indices(std::vector<size_t>())
{
}

auto
CT4CEriTensorDistributor::set_indices(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
{
    // set up local indices

    _a_loc_indices = t4cfunc::masked_indices(bra_gto_pair_block.bra_orbital_indices());

    _b_loc_indices = t4cfunc::masked_indices(bra_gto_pair_block.ket_orbital_indices());

    _c_loc_indices = t4cfunc::masked_indices(ket_gto_pair_block.bra_orbital_indices());

    _d_loc_indices = t4cfunc::masked_indices(ket_gto_pair_block.ket_orbital_indices());

    // set up global indices

    _a_glob_indices = t4cfunc::compresed_indices(bra_gto_pair_block.bra_orbital_indices());

    _b_glob_indices = t4cfunc::compresed_indices(bra_gto_pair_block.ket_orbital_indices());

    _c_glob_indices = t4cfunc::compresed_indices(ket_gto_pair_block.bra_orbital_indices());

    _d_glob_indices = t4cfunc::compresed_indices(ket_gto_pair_block.ket_orbital_indices());

    // set up local matrices

    const auto bra_ang_moms = bra_gto_pair_block.angular_momentums();

    const auto ket_ang_moms = ket_gto_pair_block.angular_momentums();

    const auto a_dims = _a_loc_indices[0] * tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.first});

    const auto b_dims = _b_loc_indices[0] * tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.second});

    const auto c_dims = _c_loc_indices[0] * tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.first});

    const auto d_dims = _d_loc_indices[0] * tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.second});
}

auto
CT4CEriTensorDistributor::distribute(const CSimdArray<double>&        buffer,
                                     const size_t                     offset,
                                     const std::vector<size_t>&       a_indices,
                                     const std::vector<size_t>&       b_indices,
                                     const std::vector<size_t>&       c_indices,
                                     const std::vector<size_t>&       d_indices,
                                     const int                        a_angmom,
                                     const int                        b_angmom,
                                     const int                        c_angmom,
                                     const int                        d_angmom,
                                     const size_t                     ibra_gto,
                                     const std::pair<size_t, size_t>& ket_range,
                                     const bool                       diagonal) -> void
{
    t4cfunc::local_distribute_eri_tensor(_eri_tensor,
                                         buffer,
                                         offset,
                                         a_indices,
                                         b_indices,
                                         c_indices,
                                         d_indices,
                                         _a_loc_indices,
                                         _b_loc_indices,
                                         _c_loc_indices,
                                         _d_loc_indices,
                                         a_angmom,
                                         b_angmom,
                                         c_angmom,
                                         d_angmom,
                                         ibra_gto,
                                         ket_range);
}

auto
CT4CEriTensorDistributor::need_omega() const -> bool
{
    return false;
}

auto
CT4CEriTensorDistributor::get_omega() const -> double
{
    return 0.0;
}
