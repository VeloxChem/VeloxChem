#include "T4CMatricesDistributor.hpp"

#include <algorithm>
#include <ranges>

#include "StringFormat.hpp"
#include "T4CLocalDistributor.hpp"
#include "T4CUtils.hpp"
#include "TensorComponents.hpp"

CT4CMatricesDistributor::CT4CMatricesDistributor(CMatrices*                      focks,
                                                 const CMatrices*                densities,
                                                 const std::vector<std::string>& labels,
                                                 const double                    exchange_factor,
                                                 const double                    omega)

    : _focks(focks)

    , _densities(densities)

    , _labels(labels)

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
    std::ranges::for_each(_labels, [](auto& label) { label = format::lower_case(label); });
}

auto
CT4CMatricesDistributor::get_omega() const -> double
{
    return _omega;
}

auto
CT4CMatricesDistributor::need_omega() const -> bool
{
    for (const auto& label : _labels)
    {
        if ((label != "j_rs") || (label != "k_rs") || (label != "kx_rs"))
        {
            return false;
        }
    }

    return true;
}

auto
CT4CMatricesDistributor::set_indices(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
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

    auto keys = _densities->keys();

    if (const auto nkeys = keys.size(); (nkeys > 0) && (nkeys == _labels.size()))
    {
        std::ranges::for_each(std::views::iota(size_t{0}, nkeys), [&](const auto i) {
            t4cfunc::add_local_matrices(
                _matrices, _labels[i], _densities->matrix(keys[i])->get_type(), std::to_string(i), a_dims, b_dims, c_dims, d_dims);
        });
    }

    _matrices.zero();
}

auto
CT4CMatricesDistributor::distribute(const CSimdArray<double>&        buffer,
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
    auto keys = _densities->keys();

    if (const auto nkeys = keys.size(); nkeys > 0)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, nkeys), [&](const auto i) {
            t4cfunc::local_distribute(_matrices,
                                      std::to_string(i),
                                      _densities->matrix(keys[i]),
                                      _labels[i],
                                      _exchange_factor,
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
                                      ket_range,
                                      diagonal);
        });
    }
}

auto
CT4CMatricesDistributor::accumulate(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
{
#pragma omp critical
    {
        const auto bra_ang_moms = bra_gto_pair_block.angular_momentums();

        const auto ket_ang_moms = ket_gto_pair_block.angular_momentums();

        const auto a_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.first});

        const auto b_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.second});

        const auto c_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.first});

        const auto d_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.second});

        auto keys = _densities->keys();

        if (const auto nkeys = keys.size(); nkeys > 0)
        {
            for (size_t i = 0; i < nkeys; i++)
            {
                const auto suffix = std::to_string(i);

                const auto label = _labels[i];

                auto fock = _focks->matrix(keys[i]);

                auto density = _densities->matrix(keys[i]);

                if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
                {
                    // acummulate contributions to Fock matrix

                    t4cfunc::accumulate(fock->sub_matrix(bra_ang_moms),
                                        _matrices.matrix("PQ_" + suffix)->sub_matrix({0, 0}),
                                        _a_loc_indices,
                                        _b_loc_indices,
                                        _a_glob_indices,
                                        _b_glob_indices,
                                        a_comps,
                                        b_comps,
                                        fock->is_angular_order(bra_ang_moms));

                    t4cfunc::accumulate(fock->sub_matrix(ket_ang_moms),
                                        _matrices.matrix("RS_" + suffix)->sub_matrix({0, 0}),
                                        _c_loc_indices,
                                        _d_loc_indices,
                                        _c_glob_indices,
                                        _d_glob_indices,
                                        c_comps,
                                        d_comps,
                                        fock->is_angular_order(ket_ang_moms));

                    if (density->get_type() == mat_t::general)
                    {
                        // set up angular pairs

                        const auto qp_pair = std::pair<int, int>({bra_ang_moms.second, bra_ang_moms.first});

                        const auto sr_pair = std::pair<int, int>({ket_ang_moms.second, ket_ang_moms.first});

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
                }

                if ((label == "2jk") || (label == "2jkx") || (label == "k") || (label == "kx") || (label == "k_rs") || (label == "kx_rs"))
                {
                    // set up angular pairs

                    const auto pr_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.first});

                    const auto ps_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.second});

                    const auto qr_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.first});

                    const auto qs_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.second});

                    // acummulate contributions to Fock matrix

                    t4cfunc::accumulate(fock->sub_matrix(pr_pair),
                                        _matrices.matrix("PR_" + suffix)->sub_matrix({0, 0}),
                                        _a_loc_indices,
                                        _c_loc_indices,
                                        _a_glob_indices,
                                        _c_glob_indices,
                                        a_comps,
                                        c_comps,
                                        fock->is_angular_order(pr_pair));

                    t4cfunc::accumulate(fock->sub_matrix(ps_pair),
                                        _matrices.matrix("PS_" + suffix)->sub_matrix({0, 0}),
                                        _a_loc_indices,
                                        _d_loc_indices,
                                        _a_glob_indices,
                                        _d_glob_indices,
                                        a_comps,
                                        d_comps,
                                        fock->is_angular_order(ps_pair));

                    t4cfunc::accumulate(fock->sub_matrix(qr_pair),
                                        _matrices.matrix("QR_" + suffix)->sub_matrix({0, 0}),
                                        _b_loc_indices,
                                        _c_loc_indices,
                                        _b_glob_indices,
                                        _c_glob_indices,
                                        b_comps,
                                        c_comps,
                                        fock->is_angular_order(qr_pair));

                    t4cfunc::accumulate(fock->sub_matrix(qs_pair),
                                        _matrices.matrix("QS_" + suffix)->sub_matrix({0, 0}),
                                        _b_loc_indices,
                                        _d_loc_indices,
                                        _b_glob_indices,
                                        _d_glob_indices,
                                        b_comps,
                                        d_comps,
                                        fock->is_angular_order(qs_pair));

                    if (density->get_type() == mat_t::general)
                    {
                        // set up angular pairs

                        const auto rp_pair = std::pair<int, int>({ket_ang_moms.first, bra_ang_moms.first});

                        const auto sp_pair = std::pair<int, int>({ket_ang_moms.second, bra_ang_moms.first});

                        const auto rq_pair = std::pair<int, int>({ket_ang_moms.first, bra_ang_moms.second});

                        const auto sq_pair = std::pair<int, int>({ket_ang_moms.second, bra_ang_moms.second});

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
}
