#include "T4CGeomX0MatricesDistributor.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>

#include "StringFormat.hpp"
#include "T4CLocalDistributor.hpp"
#include "T4CUtils.hpp"
#include "TensorComponents.hpp"

CT4CGeomX0MatricesDistributor::CT4CGeomX0MatricesDistributor(CMatrices*         focks,
                                                             const CMatrix*     density,
                                                             const std::string& label,
                                                             const double       exchange_factor,
                                                             const double       omega)

    : _focks(focks)

    , _density(density)

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
CT4CGeomX0MatricesDistributor::get_omega() const -> double
{
    return _omega;
}

auto
CT4CGeomX0MatricesDistributor::need_omega() const -> bool
{
    return (_label == "j_rs") || (_label == "k_rs") || (_label == "kx_rs");
}

auto
CT4CGeomX0MatricesDistributor::set_indices(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
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

    auto keys = _focks->keys();

    if (const auto nkeys = keys.size(); nkeys > 0)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, nkeys), [&](const auto i) {
            t4cfunc::add_local_matrices(_matrices, _label, mat_t::general, std::to_string(i), a_dims, b_dims, c_dims, d_dims);
        });
    }

    _matrices.zero();
}

auto
CT4CGeomX0MatricesDistributor::distribute(const CSimdArray<double>&        buffer,
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
    auto keys = _focks->keys();

    if (const auto nkeys = keys.size(); nkeys > 0)
    {
        const auto ncomps = tensor::number_of_spherical_components(std::array<int, 4>{a_angmom, b_angmom, c_angmom, d_angmom});

        std::ranges::for_each(std::views::iota(size_t{0}, nkeys), [&](const auto i) {
            std::cout << " *** Offset : " << keys[i] << " : " << offset + i * ncomps;

            std::cout << " ang_mom : " << a_angmom << " " << b_angmom << " " << c_angmom << " " << d_angmom << std::endl;

            t4cfunc::local_distribute_geom_ket_symm(_matrices,
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
}

auto
CT4CGeomX0MatricesDistributor::accumulate(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
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
