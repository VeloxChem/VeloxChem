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

#include "T4CGeomXYMatricesDistributor.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>

#include "StringFormat.hpp"
#include "T4CLocalDistributor.hpp"
#include "T4CUtils.hpp"
#include "TensorComponents.hpp"

CT4CGeomXYMatricesDistributor::CT4CGeomXYMatricesDistributor(CMatrices*         focks,
                                                             const CMatrix*     density,
                                                             const std::string& label,
                                                             const double       exchange_factor,
                                                             const double       omega)

    : _focks(focks)

    , _density(density)

    , _density_2(nullptr)

    , _label(format::lower_case(label))

    , _exchange_factor(exchange_factor)

    , _omega(omega)

    , _matrices(CMatrices())

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

CT4CGeomXYMatricesDistributor::CT4CGeomXYMatricesDistributor(const CMatrix*     density,
                                                             const CMatrix*     density_2,
                                                             const std::string& label,
                                                             const double       exchange_factor,
                                                             const double       omega)

    : _focks(nullptr)

    , _density(density)

    , _density_2(density_2)

    , _label(format::lower_case(label))

    , _exchange_factor(exchange_factor)

    , _omega(omega)

    , _matrices(CMatrices())

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
CT4CGeomXYMatricesDistributor::get_omega() const -> double
{
    return _omega;
}

auto
CT4CGeomXYMatricesDistributor::need_omega() const -> bool
{
    return (_label == "j_rs") || (_label == "k_rs") || (_label == "kx_rs");
}

auto
CT4CGeomXYMatricesDistributor::get_values() const -> std::vector<double>
{
    return _values;
}

auto
CT4CGeomXYMatricesDistributor::set_num_values(int num_values) -> void
{
    _values = std::vector<double>(num_values, 0.0);
}

auto
CT4CGeomXYMatricesDistributor::set_indices(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
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

    // adds submatrices storage

    if (_focks != nullptr)
    {
        auto keys = _focks->keys();

        if (const auto nkeys = keys.size(); nkeys > 0)
        {
            std::ranges::for_each(std::views::iota(size_t{0}, nkeys), [&](const auto i) {
                t4cfunc::add_local_matrices(_matrices, _label, mat_t::general, std::to_string(i), a_dims, b_dims, c_dims, d_dims);
            });
        }

        _matrices.zero();
    }
}

auto
CT4CGeomXYMatricesDistributor::distribute(const CSimdArray<double>&        buffer,
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
                                          const std::pair<size_t, size_t>& ket_range) -> void
{
    if (_focks != nullptr)
    {
        auto keys = _focks->keys();

        const auto nkeys = keys.size();

        const auto ncomps = tensor::number_of_spherical_components(std::array<int, 4>{a_angmom, b_angmom, c_angmom, d_angmom});

        std::ranges::for_each(std::views::iota(size_t{0}, nkeys), [&](const auto i) {
            t4cfunc::local_distribute_geom_ket_no_symm(_matrices,
                                                       std::to_string(i),
                                                       _density,
                                                       _label,
                                                       _exchange_factor,
                                                       buffer,
                                                       offset + i * ncomps,
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
        });
    }
    else
    {
        const auto ncomps = tensor::number_of_spherical_components(std::array<int, 4>{a_angmom, b_angmom, c_angmom, d_angmom});

        for (int i = 0; i < static_cast<int>(_values.size()); i++)
        {
            t4cfunc::local_distribute_geom_no_symm(_values,
                                           i,
                                           _density,
                                           _density_2,
                                           _label,
                                           _exchange_factor,
                                           buffer,
                                           offset + i * ncomps,
                                           a_indices,
                                           b_indices,
                                           c_indices,
                                           d_indices,
                                           a_angmom,
                                           b_angmom,
                                           c_angmom,
                                           d_angmom,
                                           ibra_gto,
                                           ket_range);
        };
    }
}

auto
CT4CGeomXYMatricesDistributor::accumulate(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
{
#pragma omp critical
    {
        const auto bra_ang_moms = bra_gto_pair_block.angular_momentums();

        const auto ket_ang_moms = ket_gto_pair_block.angular_momentums();

        const auto a_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.first});

        const auto b_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.second});

        const auto c_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.first});

        const auto d_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.second});

        auto keys = _focks->keys();

        if (const auto nkeys = keys.size(); nkeys > 0)
        {
            for (size_t i = 0; i < nkeys; i++)
            {
                const auto suffix = std::to_string(i);

                auto fock = _focks->matrix(keys[i]);

                if ((_label == "2jk") || (_label == "2jkx") || (_label == "j") || (_label == "j_rs"))
                {
                    const auto qp_pair = std::pair<int, int>({bra_ang_moms.second, bra_ang_moms.first});

                    const auto sr_pair = std::pair<int, int>({ket_ang_moms.second, ket_ang_moms.first});

                    // acummulate contributions to Fock matrix

                    t4cfunc::accumulate(fock->sub_matrix(bra_ang_moms),
                                        _matrices.matrix("PQ_" + suffix)->sub_matrix({0, 0}),
                                        _a_loc_indices,
                                        _b_loc_indices,
                                        _a_glob_indices,
                                        _b_glob_indices,
                                        a_comps,
                                        b_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(ket_ang_moms),
                                        _matrices.matrix("RS_" + suffix)->sub_matrix({0, 0}),
                                        _c_loc_indices,
                                        _d_loc_indices,
                                        _c_glob_indices,
                                        _d_glob_indices,
                                        c_comps,
                                        d_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(qp_pair),
                                        _matrices.matrix("QP_" + suffix)->sub_matrix({0, 0}),
                                        _b_loc_indices,
                                        _a_loc_indices,
                                        _b_glob_indices,
                                        _a_glob_indices,
                                        b_comps,
                                        a_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(sr_pair),
                                        _matrices.matrix("SR_" + suffix)->sub_matrix({0, 0}),
                                        _d_loc_indices,
                                        _c_loc_indices,
                                        _d_glob_indices,
                                        _c_glob_indices,
                                        d_comps,
                                        c_comps,
                                        true);
                }

                if ((_label == "2jk") || (_label == "2jkx") || (_label == "k") || (_label == "kx") || (_label == "k_rs") || (_label == "kx_rs"))
                {
                    // set up angular pairs

                    const auto pr_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.first});

                    const auto ps_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.second});

                    const auto qr_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.first});

                    const auto qs_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.second});

                    const auto rp_pair = std::pair<int, int>({ket_ang_moms.first, bra_ang_moms.first});

                    const auto sp_pair = std::pair<int, int>({ket_ang_moms.second, bra_ang_moms.first});

                    const auto rq_pair = std::pair<int, int>({ket_ang_moms.first, bra_ang_moms.second});

                    const auto sq_pair = std::pair<int, int>({ket_ang_moms.second, bra_ang_moms.second});

                    // acummulate contributions to Fock matrix

                    t4cfunc::accumulate(fock->sub_matrix(pr_pair),
                                        _matrices.matrix("PR_" + suffix)->sub_matrix({0, 0}),
                                        _a_loc_indices,
                                        _c_loc_indices,
                                        _a_glob_indices,
                                        _c_glob_indices,
                                        a_comps,
                                        c_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(ps_pair),
                                        _matrices.matrix("PS_" + suffix)->sub_matrix({0, 0}),
                                        _a_loc_indices,
                                        _d_loc_indices,
                                        _a_glob_indices,
                                        _d_glob_indices,
                                        a_comps,
                                        d_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(qr_pair),
                                        _matrices.matrix("QR_" + suffix)->sub_matrix({0, 0}),
                                        _b_loc_indices,
                                        _c_loc_indices,
                                        _b_glob_indices,
                                        _c_glob_indices,
                                        b_comps,
                                        c_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(qs_pair),
                                        _matrices.matrix("QS_" + suffix)->sub_matrix({0, 0}),
                                        _b_loc_indices,
                                        _d_loc_indices,
                                        _b_glob_indices,
                                        _d_glob_indices,
                                        b_comps,
                                        d_comps,
                                        true);

                    // acummulate contributions to Fock matrix

                    t4cfunc::accumulate(fock->sub_matrix(rp_pair),
                                        _matrices.matrix("RP_" + suffix)->sub_matrix({0, 0}),
                                        _c_loc_indices,
                                        _a_loc_indices,
                                        _c_glob_indices,
                                        _a_glob_indices,
                                        c_comps,
                                        a_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(sp_pair),
                                        _matrices.matrix("SP_" + suffix)->sub_matrix({0, 0}),
                                        _d_loc_indices,
                                        _a_loc_indices,
                                        _d_glob_indices,
                                        _a_glob_indices,
                                        d_comps,
                                        a_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(rq_pair),
                                        _matrices.matrix("RQ_" + suffix)->sub_matrix({0, 0}),
                                        _c_loc_indices,
                                        _b_loc_indices,
                                        _c_glob_indices,
                                        _b_glob_indices,
                                        c_comps,
                                        b_comps,
                                        true);

                    t4cfunc::accumulate(fock->sub_matrix(sq_pair),
                                        _matrices.matrix("SQ_" + suffix)->sub_matrix({0, 0}),
                                        _d_loc_indices,
                                        _b_loc_indices,
                                        _d_glob_indices,
                                        _b_glob_indices,
                                        d_comps,
                                        b_comps,
                                        true);
                }
            }
        }
    }
}
