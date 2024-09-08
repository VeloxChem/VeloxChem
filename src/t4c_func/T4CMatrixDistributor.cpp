#include "T4CMatrixDistributor.hpp"

#include <iostream>

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

    const auto ftype = _density->get_type();

    // Coulomb matrices

    if ((_label == "2jk") || (_label == "2jkx") || (_label == "j"))
    {
        _matrices.add(CMatrix(
                          {
                              {0, 0},
                          },
                          {
                              CSubMatrix({0, 0, a_dims, b_dims}),
                          },
                          ftype),
                      "PQ");

        _matrices.add(CMatrix(
                          {
                              {0, 0},
                          },
                          {
                              CSubMatrix({0, 0, c_dims, d_dims}),
                          },
                          ftype),
                      "RS");

        if (ftype == mat_t::general)
        {
            _matrices.add(CMatrix(
                              {
                                  {0, 0},
                              },
                              {
                                  CSubMatrix({0, 0, b_dims, a_dims}),
                              },
                              ftype),
                          "QP");

            _matrices.add(CMatrix(
                              {
                                  {0, 0},
                              },
                              {
                                  CSubMatrix({0, 0, d_dims, c_dims}),
                              },
                              ftype),
                          "SR");
        }
    }

    // Exchange matrices

    if ((_label == "2jk") || (_label == "2jkx") || (_label == "k") || (_label == "kx"))
    {
        _matrices.add(CMatrix(
                          {
                              {0, 0},
                          },
                          {
                              CSubMatrix({0, 0, a_dims, c_dims}),
                          },
                          ftype),
                      "PR");

        _matrices.add(CMatrix(
                          {
                              {0, 0},
                          },
                          {
                              CSubMatrix({0, 0, a_dims, d_dims}),
                          },
                          ftype),
                      "PS");

        _matrices.add(CMatrix(
                          {
                              {0, 0},
                          },
                          {
                              CSubMatrix({0, 0, b_dims, c_dims}),
                          },
                          ftype),
                      "QR");

        _matrices.add(CMatrix(
                          {
                              {0, 0},
                          },
                          {
                              CSubMatrix({0, 0, b_dims, d_dims}),
                          },
                          ftype),
                      "QS");

        if (ftype == mat_t::general)
        {
            _matrices.add(CMatrix(
                              {
                                  {0, 0},
                              },
                              {
                                  CSubMatrix({0, 0, c_dims, a_dims}),
                              },
                              ftype),
                          "RP");

            _matrices.add(CMatrix(
                              {
                                  {0, 0},
                              },
                              {
                                  CSubMatrix({0, 0, d_dims, a_dims}),
                              },
                              ftype),
                          "SP");

            _matrices.add(CMatrix(
                              {
                                  {0, 0},
                              },
                              {
                                  CSubMatrix({0, 0, c_dims, b_dims}),
                              },
                              ftype),
                          "RQ");

            _matrices.add(CMatrix(
                              {
                                  {0, 0},
                              },
                              {
                                  CSubMatrix({0, 0, d_dims, b_dims}),
                              },
                              ftype),
                          "SQ");
        }
    }

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
    if (_density->get_type() == mat_t::symmetric)
    {
        if (_label == "2jk")
        {
            t4cfunc::local_distribute_rest_jk(_matrices,
                                              _density,
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

        if (_label == "2jkx")
        {
            t4cfunc::local_distribute_rest_jkx(_matrices,
                                               _density,
                                               buffer,
                                               offset,
                                               _exchange_factor,
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

        if (_label == "j")
        {
            t4cfunc::local_distribute_rest_j(_matrices,
                                             _density,
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

        if (_label == "k")
        {
            t4cfunc::local_distribute_rest_k(_matrices,
                                             _density,
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

        if (_label == "kx")
        {
            t4cfunc::local_distribute_rest_kx(_matrices,
                                              _density,
                                              buffer,
                                              offset,
                                              _exchange_factor,
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
    }
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

        if (_density->get_type() == mat_t::symmetric)
        {
            if ((_label == "2jk") || (_label == "2jkx") || (_label == "j"))
            {
                // set up Fock submatrices

                auto submat_pq = _fock->sub_matrix(bra_ang_moms);

                auto submat_rs = _fock->sub_matrix(ket_ang_moms);

                // set up angular orders

                const auto angord_pq = _fock->is_angular_order(bra_ang_moms);

                const auto angord_rs = _fock->is_angular_order(ket_ang_moms);

                // acummulate contributions to Fock matrix

                t4cfunc::accumulate(submat_pq,
                                    _matrices.matrix("PQ")->sub_matrix({0, 0}),
                                    _a_loc_indices,
                                    _b_loc_indices,
                                    _a_glob_indices,
                                    _b_glob_indices,
                                    a_comps,
                                    b_comps,
                                    angord_pq);

                t4cfunc::accumulate(submat_rs,
                                    _matrices.matrix("RS")->sub_matrix({0, 0}),
                                    _c_loc_indices,
                                    _d_loc_indices,
                                    _c_glob_indices,
                                    _d_glob_indices,
                                    c_comps,
                                    d_comps,
                                    angord_rs);
            }

            if ((_label == "2jk") || (_label == "2jkx") || (_label == "k") || (_label == "kx"))
            {
                // set up angular pairs

                const auto pr_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.first});

                const auto ps_pair = std::pair<int, int>({bra_ang_moms.first, ket_ang_moms.second});

                const auto qr_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.first});

                const auto qs_pair = std::pair<int, int>({bra_ang_moms.second, ket_ang_moms.second});

                // set up Fock submatrices

                auto submat_pr = _fock->sub_matrix(pr_pair);

                auto submat_ps = _fock->sub_matrix(ps_pair);

                auto submat_qr = _fock->sub_matrix(qr_pair);

                auto submat_qs = _fock->sub_matrix(qs_pair);

                // set up angular orders

                const auto angord_pr = _fock->is_angular_order(pr_pair);

                const auto angord_ps = _fock->is_angular_order(ps_pair);

                const auto angord_qr = _fock->is_angular_order(qr_pair);

                const auto angord_qs = _fock->is_angular_order(qs_pair);

                // acummulate contributions to Fock matrix

                t4cfunc::accumulate(submat_pr,
                                    _matrices.matrix("PR")->sub_matrix({0, 0}),
                                    _a_loc_indices,
                                    _c_loc_indices,
                                    _a_glob_indices,
                                    _c_glob_indices,
                                    a_comps,
                                    c_comps,
                                    angord_pr);

                t4cfunc::accumulate(submat_ps,
                                    _matrices.matrix("PS")->sub_matrix({0, 0}),
                                    _a_loc_indices,
                                    _d_loc_indices,
                                    _a_glob_indices,
                                    _d_glob_indices,
                                    a_comps,
                                    d_comps,
                                    angord_ps);

                t4cfunc::accumulate(submat_qr,
                                    _matrices.matrix("QR")->sub_matrix({0, 0}),
                                    _b_loc_indices,
                                    _c_loc_indices,
                                    _b_glob_indices,
                                    _c_glob_indices,
                                    b_comps,
                                    c_comps,
                                    angord_qr);

                t4cfunc::accumulate(submat_qs,
                                    _matrices.matrix("QS")->sub_matrix({0, 0}),
                                    _b_loc_indices,
                                    _d_loc_indices,
                                    _b_glob_indices,
                                    _d_glob_indices,
                                    b_comps,
                                    d_comps,
                                    angord_qs);
            }
        }
    }
}
