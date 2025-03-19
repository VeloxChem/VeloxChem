#ifndef TwoCenterElectronRepulsionGeom100RecIG_hpp
#define TwoCenterElectronRepulsionGeom100RecIG_hpp

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
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFS.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "TwoCenterElectronRepulsionPrimRecHD.hpp"
#include "TwoCenterElectronRepulsionPrimRecHF.hpp"
#include "TwoCenterElectronRepulsionPrimRecHG.hpp"
#include "TwoCenterElectronRepulsionPrimRecIF.hpp"
#include "TwoCenterElectronRepulsionPrimRecIG.hpp"
#include "TwoCenterElectronRepulsionPrimRecKG.hpp"
#include "GeometricalDerivatives1X0ForIY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)I|1/|r-r'||G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_ig(T& distributor,
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

    CSimdArray<double> pbuffer(8100, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(1260, 1);

    CSimdArray<double> sbuffer(351, 1);

    // setup Boys function data

    const CBoysFunc<11> bf_table;

    CSimdArray<double> bf_data(13, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 11, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 11);

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

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 11, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 14, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 17, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 20, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 23, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 26, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 29, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 32, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 35, 9, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 38, 10, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 41, 0, 1, 14, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 47, 1, 2, 17, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 53, 2, 3, 20, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 59, 3, 4, 23, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 65, 4, 5, 26, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 71, 5, 6, 29, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 77, 6, 7, 32, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 83, 7, 8, 35, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 89, 8, 9, 38, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 95, 11, 14, 47, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 105, 14, 17, 53, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 115, 17, 20, 59, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 125, 20, 23, 65, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 135, 23, 26, 71, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 145, 26, 29, 77, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 155, 29, 32, 83, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 165, 32, 35, 89, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 175, 41, 47, 105, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 190, 47, 53, 115, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 205, 53, 59, 125, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 220, 59, 65, 135, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 235, 65, 71, 145, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 250, 71, 77, 155, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 265, 77, 83, 165, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 280, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 283, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 286, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 289, 4, 23, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 298, 5, 26, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 307, 6, 29, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 316, 17, 53, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 334, 20, 59, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 352, 23, 65, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 370, 26, 71, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 388, 29, 77, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 406, 53, 115, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 436, 59, 125, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 466, 65, 135, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 496, 71, 145, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 526, 77, 155, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 556, 95, 175, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 601, 105, 190, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 646, 115, 205, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 691, 125, 220, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 736, 135, 235, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 781, 145, 250, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 826, 155, 265, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 871, 4, 5, 286, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 877, 17, 20, 280, 289, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 895, 20, 23, 283, 298, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 913, 23, 26, 286, 307, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 931, 53, 59, 289, 352, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 967, 59, 65, 298, 370, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 1003, 65, 71, 307, 388, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1039, 95, 105, 316, 406, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1099, 105, 115, 334, 436, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1159, 115, 125, 352, 466, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1219, 125, 135, 370, 496, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1279, 135, 145, 388, 526, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1339, 175, 190, 406, 646, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1429, 190, 205, 436, 691, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1519, 205, 220, 466, 736, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1609, 220, 235, 496, 781, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1699, 235, 250, 526, 826, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 1789, 280, 283, 871, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 1799, 289, 298, 871, 913, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1829, 316, 334, 877, 931, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1889, 334, 352, 895, 967, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1949, 352, 370, 913, 1003, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 2009, 406, 436, 931, 1159, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 2109, 436, 466, 967, 1219, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 2209, 466, 496, 1003, 1279, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2309, 556, 601, 1039, 1339, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2459, 601, 646, 1099, 1429, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2609, 646, 691, 1159, 1519, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2759, 691, 736, 1219, 1609, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2909, 736, 781, 1279, 1699, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 3059, 877, 895, 1789, 1799, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 3104, 931, 967, 1799, 1949, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 3194, 1039, 1099, 1829, 2009, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 3344, 1099, 1159, 1889, 2109, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 3494, 1159, 1219, 1949, 2209, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 3644, 1339, 1429, 2009, 2609, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 3869, 1429, 1519, 2109, 2759, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 4094, 1519, 1609, 2209, 2909, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 4319, 1829, 1889, 3059, 3104, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hf(pbuffer, 4445, 2009, 2109, 3104, 3494, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 4655, 2309, 2459, 3194, 3644, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 4970, 2459, 2609, 3344, 3869, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 5285, 2609, 2759, 3494, 4094, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_if(pbuffer, 5600, 3194, 3344, 4319, 4445, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ig(pbuffer, 5880, 3644, 3869, 4445, 5285, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_kg(pbuffer, 6300, 4655, 4970, 5600, 5880, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_ix(pbuffer, 6840, 4655, 6300, 1, 15, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 6840, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecIG_hpp */
