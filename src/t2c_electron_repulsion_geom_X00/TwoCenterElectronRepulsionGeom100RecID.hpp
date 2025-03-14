#ifndef TwoCenterElectronRepulsionGeom100RecID_hpp
#define TwoCenterElectronRepulsionGeom100RecID_hpp

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
#include "TwoCenterElectronRepulsionPrimRecFS.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGS.hpp"
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecHS.hpp"
#include "TwoCenterElectronRepulsionPrimRecHP.hpp"
#include "TwoCenterElectronRepulsionPrimRecHD.hpp"
#include "TwoCenterElectronRepulsionPrimRecIP.hpp"
#include "TwoCenterElectronRepulsionPrimRecID.hpp"
#include "TwoCenterElectronRepulsionPrimRecKD.hpp"
#include "GeometricalDerivatives1X0ForIY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)I|1/|r-r'||D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_id(T& distributor,
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

    CSimdArray<double> pbuffer(2823, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(504, 1);

    CSimdArray<double> sbuffer(195, 1);

    // setup Boys function data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 9, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 9);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 9, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 12, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 15, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 18, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 21, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 24, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 27, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 30, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 33, 0, 1, 12, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 39, 1, 2, 15, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 45, 2, 3, 18, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 51, 3, 4, 21, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 57, 4, 5, 24, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 63, 5, 6, 27, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 69, 6, 7, 30, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 75, 2, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 78, 3, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 81, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 84, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 87, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 90, 2, 15, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 99, 3, 18, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 108, 4, 21, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 117, 5, 24, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 126, 6, 27, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 135, 9, 33, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 153, 12, 39, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 171, 15, 45, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 189, 18, 51, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 207, 21, 57, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 225, 24, 63, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 243, 27, 69, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 261, 2, 3, 81, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 267, 3, 4, 84, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 273, 4, 5, 87, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 279, 9, 12, 75, 90, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 297, 12, 15, 78, 99, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 315, 15, 18, 81, 108, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 333, 18, 21, 84, 117, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 351, 21, 24, 87, 126, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 369, 33, 39, 90, 171, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 405, 39, 45, 99, 189, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 441, 45, 51, 108, 207, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 477, 51, 57, 117, 225, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 513, 57, 63, 126, 243, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 549, 75, 78, 261, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 559, 78, 81, 267, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 569, 81, 84, 273, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 579, 90, 99, 261, 315, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 609, 99, 108, 267, 333, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 639, 108, 117, 273, 351, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 669, 135, 153, 279, 369, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 729, 153, 171, 297, 405, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 789, 171, 189, 315, 441, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 849, 189, 207, 333, 477, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 909, 207, 225, 351, 513, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 969, 261, 267, 569, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 984, 279, 297, 549, 579, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 1029, 297, 315, 559, 609, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 1074, 315, 333, 569, 639, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 1119, 369, 405, 579, 789, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 1209, 405, 441, 609, 849, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 1299, 441, 477, 639, 909, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hs(pbuffer, 1389, 549, 559, 969, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hp(pbuffer, 1410, 579, 609, 969, 1074, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 1473, 669, 729, 984, 1119, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 1599, 729, 789, 1029, 1209, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 1725, 789, 849, 1074, 1299, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ip(pbuffer, 1851, 984, 1029, 1389, 1410, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_id(pbuffer, 1935, 1119, 1209, 1410, 1725, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_kd(pbuffer, 2103, 1473, 1599, 1851, 1935, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_ix(pbuffer, 2319, 1473, 2103, 1, 6, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 2319, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 2, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecID_hpp */
