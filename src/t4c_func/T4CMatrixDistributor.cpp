#include "T4CMatrixDistributor.hpp"

#include "StringFormat.hpp"
#include "T4CDistributor.hpp"
#include "T4CUtils.hpp"
#include "TensorComponents.hpp"
#include "T4CLocalDistributor.hpp"

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
    }
}
