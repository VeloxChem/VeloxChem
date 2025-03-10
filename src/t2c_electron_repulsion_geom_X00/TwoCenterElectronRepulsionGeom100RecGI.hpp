#ifndef TwoCenterElectronRepulsionGeom100RecGI_hpp
#define TwoCenterElectronRepulsionGeom100RecGI_hpp

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
#include "TwoCenterElectronRepulsionPrimRecSI.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecPI.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDI.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFI.hpp"
#include "TwoCenterElectronRepulsionPrimRecGH.hpp"
#include "TwoCenterElectronRepulsionPrimRecGI.hpp"
#include "TwoCenterElectronRepulsionPrimRecHI.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)G|1/|r-r'||I)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_gi(T& distributor,
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

    CSimdArray<double> pbuffer(6153, ket_npgtos);

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

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 280, 95, 105, 190, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 301, 105, 115, 205, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 322, 115, 125, 220, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 343, 125, 135, 235, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 364, 135, 145, 250, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 385, 145, 155, 265, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 406, 175, 190, 301, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 434, 190, 205, 322, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 462, 205, 220, 343, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 490, 220, 235, 364, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 518, 235, 250, 385, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 546, 23, 65, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 564, 65, 135, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 594, 115, 205, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 639, 125, 220, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 684, 135, 235, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 729, 205, 322, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 792, 220, 343, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 855, 235, 364, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 918, 280, 406, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1002, 301, 434, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1086, 322, 462, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1170, 343, 490, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1254, 364, 518, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1338, 115, 125, 546, 564, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1398, 205, 220, 564, 684, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1488, 280, 301, 594, 729, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1614, 301, 322, 639, 792, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1740, 322, 343, 684, 855, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 1866, 406, 434, 729, 1086, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 2034, 434, 462, 792, 1170, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 2202, 462, 490, 855, 1254, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2370, 594, 639, 1338, 1398, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 2520, 729, 792, 1398, 1740, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 2730, 918, 1002, 1488, 1866, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 3010, 1002, 1086, 1614, 2034, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 3290, 1086, 1170, 1740, 2202, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 3570, 1488, 1614, 2370, 2520, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 3885, 1866, 2034, 2520, 3290, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hi(pbuffer, 4305, 2730, 3010, 3570, 3885, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 4893, 2730, 4305, 1, 28, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 4893, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 6>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 6, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecGI_hpp */
