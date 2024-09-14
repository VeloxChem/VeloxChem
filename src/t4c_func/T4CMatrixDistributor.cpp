#include "T4CMatrixDistributor.hpp"

#include "StringFormat.hpp"
#include "T4CLocalDistributor.hpp"
#include "T4CUtils.hpp"
#include "TensorComponents.hpp"

CT4CMatrixDistributor::CT4CMatrixDistributor(CMatrix*           fock,
                                             const CMatrix*     density,
                                             const std::string& label,
                                             const double       exchange_factor,
                                             const double       omega)

    : _fock(fock)

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
CT4CMatrixDistributor::get_omega() const -> double
{
    return _omega;
}

auto
CT4CMatrixDistributor::set_indices(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
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

    t4cfunc::add_local_matrices(_matrices, _label, _density->get_type(), "0", a_dims, b_dims, c_dims, d_dims);

    _matrices.zero();
}

auto
CT4CMatrixDistributor::distribute(const CSimdArray<double>&        buffer,
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
    t4cfunc::local_distribute(_matrices,
                              "0",
                              _density,
                              _label,
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
}

auto
CT4CMatrixDistributor::accumulate(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void
{
#pragma omp critical
    {
        const auto bra_ang_moms = bra_gto_pair_block.angular_momentums();

        const auto ket_ang_moms = ket_gto_pair_block.angular_momentums();

        const auto a_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.first});

        const auto b_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_ang_moms.second});

        const auto c_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.first});

        const auto d_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_ang_moms.second});

        if ((_label == "2jk") || (_label == "2jkx") || (_label == "jk") || (_label == "jkx") || (_label == "j"))
        {
            // acummulate contributions to Fock matrix

            t4cfunc::accumulate(_fock->sub_matrix(bra_ang_moms),
                                _matrices.matrix("PQ_0")->sub_matrix({0, 0}),
                                _a_loc_indices,
                                _b_loc_indices,
                                _a_glob_indices,
                                _b_glob_indices,
                                a_comps,
                                b_comps,
                                _fock->is_angular_order(bra_ang_moms));

            t4cfunc::accumulate(_fock->sub_matrix(ket_ang_moms),
                                _matrices.matrix("RS_0")->sub_matrix({0, 0}),
                                _c_loc_indices,
                                _d_loc_indices,
                                _c_glob_indices,
                                _d_glob_indices,
                                c_comps,
                                d_comps,
                                _fock->is_angular_order(ket_ang_moms));

            if (_density->get_type() == mat_t::general)
            {
                // set up angular pairs

                const auto qp_pair = std::pair<int, int>({bra_ang_moms.second, bra_ang_moms.first});

                const auto sr_pair = std::pair<int, int>({ket_ang_moms.second, ket_ang_moms.first});

                t4cfunc::accumulate(_fock->sub_matrix(qp_pair),
                                    _matrices.matrix("QP_0")->sub_matrix({0, 0}),
                                    _b_loc_indices,
                                    _a_loc_indices,
                                    _b_glob_indices,
                                    _a_glob_indices,
                                    b_comps,
                                    a_comps,
                                    true);

                t4cfunc::accumulate(_fock->sub_matrix(sr_pair),
                                    _matrices.matrix("SR_0")->sub_matrix({0, 0}),
                                    _d_loc_indices,
                                    _c_loc_indices,
                                    _d_glob_indices,
                                    _c_glob_indices,
                                    d_comps,
                                    c_comps,
                                    true);
            }
        }

        if ((_label == "2jk") || (_label == "2jkx") || (_label == "jk") || (_label == "jkx") || (_label == "k") || (_label == "kx"))
        {
            // set up angular pairs

            const auto pr_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.first});

            const auto ps_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.second});

            const auto qr_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.first});

            const auto qs_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.second});

            // acummulate contributions to Fock matrix

            t4cfunc::accumulate(_fock->sub_matrix(pr_pair),
                                _matrices.matrix("PR_0")->sub_matrix({0, 0}),
                                _a_loc_indices,
                                _c_loc_indices,
                                _a_glob_indices,
                                _c_glob_indices,
                                a_comps,
                                c_comps,
                                _fock->is_angular_order(pr_pair));

            t4cfunc::accumulate(_fock->sub_matrix(ps_pair),
                                _matrices.matrix("PS_0")->sub_matrix({0, 0}),
                                _a_loc_indices,
                                _d_loc_indices,
                                _a_glob_indices,
                                _d_glob_indices,
                                a_comps,
                                d_comps,
                                _fock->is_angular_order(ps_pair));

            t4cfunc::accumulate(_fock->sub_matrix(qr_pair),
                                _matrices.matrix("QR_0")->sub_matrix({0, 0}),
                                _b_loc_indices,
                                _c_loc_indices,
                                _b_glob_indices,
                                _c_glob_indices,
                                b_comps,
                                c_comps,
                                _fock->is_angular_order(qr_pair));

            t4cfunc::accumulate(_fock->sub_matrix(qs_pair),
                                _matrices.matrix("QS_0")->sub_matrix({0, 0}),
                                _b_loc_indices,
                                _d_loc_indices,
                                _b_glob_indices,
                                _d_glob_indices,
                                b_comps,
                                d_comps,
                                _fock->is_angular_order(qs_pair));

            if (_density->get_type() == mat_t::general)
            {
                // set up angular pairs

                const auto rp_pair = std::pair<int, int>({ket_ang_moms.first, bra_ang_moms.first});

                const auto sp_pair = std::pair<int, int>({ket_ang_moms.second, bra_ang_moms.first});

                const auto rq_pair = std::pair<int, int>({ket_ang_moms.first, bra_ang_moms.second});

                const auto sq_pair = std::pair<int, int>({ket_ang_moms.second, bra_ang_moms.second});

                // acummulate contributions to Fock matrix

                t4cfunc::accumulate(_fock->sub_matrix(rp_pair),
                                    _matrices.matrix("RP_0")->sub_matrix({0, 0}),
                                    _c_loc_indices,
                                    _a_loc_indices,
                                    _c_glob_indices,
                                    _a_glob_indices,
                                    c_comps,
                                    a_comps,
                                    true);

                t4cfunc::accumulate(_fock->sub_matrix(sp_pair),
                                    _matrices.matrix("SP_0")->sub_matrix({0, 0}),
                                    _d_loc_indices,
                                    _a_loc_indices,
                                    _d_glob_indices,
                                    _a_glob_indices,
                                    d_comps,
                                    a_comps,
                                    true);

                t4cfunc::accumulate(_fock->sub_matrix(rq_pair),
                                    _matrices.matrix("RQ_0")->sub_matrix({0, 0}),
                                    _c_loc_indices,
                                    _b_loc_indices,
                                    _c_glob_indices,
                                    _b_glob_indices,
                                    c_comps,
                                    b_comps,
                                    true);

                t4cfunc::accumulate(_fock->sub_matrix(sq_pair),
                                    _matrices.matrix("SQ_0")->sub_matrix({0, 0}),
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
