#ifndef TwoCenterElectronRepulsionGeom100RecIP_hpp
#define TwoCenterElectronRepulsionGeom100RecIP_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFS.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGS.hpp"
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecHS.hpp"
#include "TwoCenterElectronRepulsionPrimRecHP.hpp"
#include "TwoCenterElectronRepulsionPrimRecIS.hpp"
#include "TwoCenterElectronRepulsionPrimRecIP.hpp"
#include "TwoCenterElectronRepulsionPrimRecKP.hpp"
#include "GeometricalDerivatives1X0ForIY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)I|1/|r-r'||P)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_ip(T& distributor,
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

    CSimdArray<double> pbuffer(1269, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(252, 1);

    CSimdArray<double> sbuffer(117, 1);

    // setup Boys function data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 8, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 8);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 8, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 11, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 14, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 17, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 20, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 23, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 26, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 29, 2, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 32, 3, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 35, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 38, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 41, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 44, 0, 8, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 53, 1, 11, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 62, 2, 14, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 71, 3, 17, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 80, 4, 20, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 89, 5, 23, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 98, 6, 26, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 107, 0, 1, 29, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 113, 1, 2, 32, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 119, 2, 3, 35, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 125, 3, 4, 38, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 131, 4, 5, 41, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 137, 8, 11, 29, 62, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 155, 11, 14, 32, 71, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 173, 14, 17, 35, 80, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 191, 17, 20, 38, 89, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 209, 20, 23, 41, 98, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 227, 29, 32, 119, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 237, 32, 35, 125, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 247, 35, 38, 131, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 257, 44, 53, 107, 137, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 287, 53, 62, 113, 155, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 317, 62, 71, 119, 173, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 347, 71, 80, 125, 191, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 377, 80, 89, 131, 209, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 407, 107, 113, 227, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 422, 113, 119, 237, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 437, 119, 125, 247, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 452, 137, 155, 227, 317, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 497, 155, 173, 237, 347, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 542, 173, 191, 247, 377, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hs(pbuffer, 587, 227, 237, 437, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hp(pbuffer, 608, 257, 287, 407, 452, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hp(pbuffer, 671, 287, 317, 422, 497, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hp(pbuffer, 734, 317, 347, 437, 542, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_is(pbuffer, 797, 407, 422, 587, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ip(pbuffer, 825, 452, 497, 587, 734, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_kp(pbuffer, 909, 608, 671, 797, 825, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_ix(pbuffer, 1017, 608, 909, 1, 3, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 1017, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 1>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 1, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecIP_hpp */
