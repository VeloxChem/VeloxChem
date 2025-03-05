#ifndef TwoCenterElectronRepulsionRecGD_hpp
#define TwoCenterElectronRepulsionRecGD_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecSD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (G|1/|r-r'||D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_gd(T& distributor,
                           const CGtoBlock& bra_gto_block,
                           const CGtoBlock& ket_gto_block,
                           const std::pair<size_t, size_t>& bra_indices,
                           const std::pair<size_t, size_t>& ket_indices,
                           const bool bra_eq_ket) -> void
{
    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.coordinates();

    const auto bra_gto_exps = bra_gto_block.exponents();

    const auto bra_gto_norms = bra_gto_block.normalization_factors();

    const auto bra_gto_indices = bra_gto_block.orbital_indices();

    const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();

    const auto bra_npgtos = bra_gto_block.number_of_primitives();

    // intialize GTOs data on ket side

    const auto ket_gto_coords = ket_gto_block.coordinates();

    const auto ket_gto_exps = ket_gto_block.exponents();

    const auto ket_gto_norms = ket_gto_block.normalization_factors();

    const auto ket_gto_indices = ket_gto_block.orbital_indices();

    const auto ket_npgtos = ket_gto_block.number_of_primitives();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> factors(14, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(448, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(90, 1);

    CSimdArray<double> sbuffer(45, 1);

    // setup Boys function data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

        factors.load(ket_gto_exps, ket_range, 0, ket_npgtos);

        factors.load(ket_gto_norms, ket_range, 1, ket_npgtos);

        factors.replicate_points(ket_gto_coords, ket_range, 2, ket_npgtos);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        sbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        pbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over contracted basis functions on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            cbuffer.zero();

            sbuffer.zero();

            const auto r_a = bra_gto_coords[j];

            t2cfunc::comp_distances_ab(factors, 5, 2, r_a);

            for (size_t k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];

                const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];

                t2cfunc::comp_distances_pa(factors, 8, 5, a_exp);

                t2cfunc::comp_distances_pb(factors, 11, 5, a_exp);

                t2cfunc::comp_boys_args_with_rho(bf_data, 7, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 7);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 0, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 7, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 10, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 13, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 16, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 19, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 22, 0, 1, 7, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 28, 1, 2, 10, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 34, 2, 3, 13, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 40, 3, 4, 16, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 46, 4, 5, 19, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 52, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 55, 2, 10, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 64, 3, 13, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 73, 4, 16, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 82, 10, 34, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 100, 13, 40, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 118, 16, 46, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 136, 2, 3, 52, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 142, 10, 13, 52, 73, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 160, 22, 28, 55, 82, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 196, 28, 34, 64, 100, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 232, 34, 40, 73, 118, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 268, 55, 64, 136, 142, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 298, 82, 100, 142, 232, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 358, 160, 196, 268, 298, factors, 8, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 358, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionRecGD_hpp */
