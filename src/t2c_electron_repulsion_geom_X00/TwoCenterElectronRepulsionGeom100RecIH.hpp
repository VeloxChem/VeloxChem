#ifndef TwoCenterElectronRepulsionGeom100RecIH_hpp
#define TwoCenterElectronRepulsionGeom100RecIH_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecSD.hpp"
#include "TwoCenterElectronRepulsionPrimRecSF.hpp"
#include "TwoCenterElectronRepulsionPrimRecSG.hpp"
#include "TwoCenterElectronRepulsionPrimRecSH.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGH.hpp"
#include "TwoCenterElectronRepulsionPrimRecHF.hpp"
#include "TwoCenterElectronRepulsionPrimRecHG.hpp"
#include "TwoCenterElectronRepulsionPrimRecHH.hpp"
#include "TwoCenterElectronRepulsionPrimRecIG.hpp"
#include "TwoCenterElectronRepulsionPrimRecIH.hpp"
#include "TwoCenterElectronRepulsionPrimRecKH.hpp"
#include "GeometricalDerivatives1X0ForIY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)I|1/|r-r'||H)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_ih(T& distributor,
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

    CSimdArray<double> pbuffer(11880, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(1764, 1);

    CSimdArray<double> sbuffer(429, 1);

    // setup Boys function data

    const CBoysFunc<12> bf_table;

    CSimdArray<double> bf_data(14, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 12, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 12);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 9, bf_data, 10, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 10, bf_data, 11, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 11, bf_data, 12, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 12, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 15, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 18, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 21, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 24, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 27, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 30, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 33, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 36, 9, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 39, 10, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 42, 11, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 45, 0, 1, 15, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 51, 1, 2, 18, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 57, 2, 3, 21, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 63, 3, 4, 24, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 69, 4, 5, 27, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 75, 5, 6, 30, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 81, 6, 7, 33, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 87, 7, 8, 36, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 93, 8, 9, 39, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 99, 9, 10, 42, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 105, 12, 15, 51, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 115, 15, 18, 57, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 125, 18, 21, 63, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 135, 21, 24, 69, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 145, 24, 27, 75, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 155, 27, 30, 81, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 165, 30, 33, 87, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 175, 33, 36, 93, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 185, 36, 39, 99, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 195, 45, 51, 115, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 210, 51, 57, 125, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 225, 57, 63, 135, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 240, 63, 69, 145, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 255, 69, 75, 155, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 270, 75, 81, 165, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 285, 81, 87, 175, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 300, 87, 93, 185, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 315, 105, 115, 210, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 336, 115, 125, 225, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 357, 125, 135, 240, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 378, 135, 145, 255, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 399, 145, 155, 270, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 420, 155, 165, 285, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 441, 165, 175, 300, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 462, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 465, 4, 24, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 474, 5, 27, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 483, 6, 30, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 492, 24, 69, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 510, 27, 75, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 528, 30, 81, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 546, 57, 125, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 576, 63, 135, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 606, 69, 145, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 636, 75, 155, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 666, 81, 165, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 696, 125, 225, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 741, 135, 240, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 786, 145, 255, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 831, 155, 270, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 876, 165, 285, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 921, 195, 315, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 984, 210, 336, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1047, 225, 357, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1110, 240, 378, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1173, 255, 399, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1236, 270, 420, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1299, 285, 441, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 1362, 4, 5, 462, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 1368, 24, 27, 462, 483, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 1386, 57, 63, 465, 492, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 1422, 63, 69, 474, 510, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 1458, 69, 75, 483, 528, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1494, 125, 135, 492, 606, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1554, 135, 145, 510, 636, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1614, 145, 155, 528, 666, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1674, 195, 210, 546, 696, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1764, 210, 225, 576, 741, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1854, 225, 240, 606, 786, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1944, 240, 255, 636, 831, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 2034, 255, 270, 666, 876, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2124, 315, 336, 696, 1047, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2250, 336, 357, 741, 1110, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2376, 357, 378, 786, 1173, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2502, 378, 399, 831, 1236, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2628, 399, 420, 876, 1299, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 2754, 465, 474, 1362, 1368, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 2784, 492, 510, 1368, 1458, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 2844, 546, 576, 1386, 1494, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 2944, 576, 606, 1422, 1554, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 3044, 606, 636, 1458, 1614, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 3144, 696, 741, 1494, 1854, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 3294, 741, 786, 1554, 1944, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 3444, 786, 831, 1614, 2034, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 3594, 921, 984, 1674, 2124, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 3804, 984, 1047, 1764, 2250, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 4014, 1047, 1110, 1854, 2376, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 4224, 1110, 1173, 1944, 2502, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 4434, 1173, 1236, 2034, 2628, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 4644, 1386, 1422, 2754, 2784, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 4734, 1494, 1554, 2784, 3044, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 4884, 1674, 1764, 2844, 3144, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 5109, 1764, 1854, 2944, 3294, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 5334, 1854, 1944, 3044, 3444, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 5559, 2124, 2250, 3144, 4014, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 5874, 2250, 2376, 3294, 4224, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 6189, 2376, 2502, 3444, 4434, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hf(pbuffer, 6504, 2844, 2944, 4644, 4734, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 6714, 3144, 3294, 4734, 5334, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hh(pbuffer, 7029, 3594, 3804, 4884, 5559, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hh(pbuffer, 7470, 3804, 4014, 5109, 5874, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hh(pbuffer, 7911, 4014, 4224, 5334, 6189, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ig(pbuffer, 8352, 4884, 5109, 6504, 6714, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ih(pbuffer, 8772, 5559, 5874, 6714, 7911, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_kh(pbuffer, 9360, 7029, 7470, 8352, 8772, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_ix(pbuffer, 10116, 7029, 9360, 1, 21, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 10116, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 5>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 5, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecIH_hpp */
