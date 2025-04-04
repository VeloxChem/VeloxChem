#ifndef TwoCenterElectronRepulsionGeom100RecHI_hpp
#define TwoCenterElectronRepulsionGeom100RecHI_hpp

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
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecPI.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDI.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFI.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGH.hpp"
#include "TwoCenterElectronRepulsionPrimRecGI.hpp"
#include "TwoCenterElectronRepulsionPrimRecHH.hpp"
#include "TwoCenterElectronRepulsionPrimRecHI.hpp"
#include "TwoCenterElectronRepulsionPrimRecII.hpp"
#include "GeometricalDerivatives1X0ForHY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)H|1/|r-r'||I)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_hi(T& distributor,
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

    CSimdArray<double> pbuffer(10348, ket_npgtos);

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

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 0, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 9, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 10, bf_data, 10, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 11, bf_data, 11, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 12, bf_data, 12, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 13, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 16, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 19, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 22, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 25, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 28, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 31, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 34, 9, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 37, 10, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 40, 11, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 43, 12, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 46, 0, 1, 13, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 52, 1, 2, 16, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 58, 2, 3, 19, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 64, 3, 4, 22, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 70, 4, 5, 25, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 76, 5, 6, 28, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 82, 6, 7, 31, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 88, 7, 8, 34, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 94, 8, 9, 37, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 100, 9, 10, 40, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 106, 10, 11, 43, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 112, 13, 16, 58, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 122, 16, 19, 64, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 132, 19, 22, 70, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 142, 22, 25, 76, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 152, 25, 28, 82, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 162, 28, 31, 88, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 172, 31, 34, 94, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 182, 34, 37, 100, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 192, 37, 40, 106, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 202, 46, 52, 112, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 217, 52, 58, 122, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 232, 58, 64, 132, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 247, 64, 70, 142, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 262, 70, 76, 152, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 277, 76, 82, 162, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 292, 82, 88, 172, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 307, 88, 94, 182, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 322, 94, 100, 192, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 337, 112, 122, 232, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 358, 122, 132, 247, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 379, 132, 142, 262, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 400, 142, 152, 277, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 421, 152, 162, 292, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 442, 162, 172, 307, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 463, 172, 182, 322, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 484, 202, 217, 337, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 512, 217, 232, 358, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 540, 232, 247, 379, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 568, 247, 262, 400, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 596, 262, 277, 421, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 624, 277, 292, 442, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 652, 292, 307, 463, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 680, 6, 28, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 689, 28, 82, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 707, 70, 142, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 737, 76, 152, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 767, 82, 162, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 797, 142, 262, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 842, 152, 277, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 887, 162, 292, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 932, 232, 358, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 995, 247, 379, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1058, 262, 400, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1121, 277, 421, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1184, 292, 442, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1247, 358, 540, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1331, 379, 568, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1415, 400, 596, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1499, 421, 624, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1583, 442, 652, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 1667, 70, 76, 680, 689, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1703, 142, 152, 689, 767, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1763, 232, 247, 707, 797, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1853, 247, 262, 737, 842, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1943, 262, 277, 767, 887, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2033, 358, 379, 797, 1058, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2159, 379, 400, 842, 1121, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2285, 400, 421, 887, 1184, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 2411, 484, 512, 932, 1247, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 2579, 512, 540, 995, 1331, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 2747, 540, 568, 1058, 1415, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 2915, 568, 596, 1121, 1499, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 3083, 596, 624, 1184, 1583, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 3251, 707, 737, 1667, 1703, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 3351, 797, 842, 1703, 1943, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 3501, 932, 995, 1763, 2033, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 3711, 995, 1058, 1853, 2159, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 3921, 1058, 1121, 1943, 2285, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 4131, 1247, 1331, 2033, 2747, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 4411, 1331, 1415, 2159, 2915, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 4691, 1415, 1499, 2285, 3083, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 4971, 1763, 1853, 3251, 3351, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 5196, 2033, 2159, 3351, 3921, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 5511, 2411, 2579, 3501, 4131, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 5931, 2579, 2747, 3711, 4411, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 6351, 2747, 2915, 3921, 4691, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hh(pbuffer, 6771, 3501, 3711, 4971, 5196, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hi(pbuffer, 7212, 4131, 4411, 5196, 6351, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ii(pbuffer, 7800, 5511, 5931, 6771, 7212, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_hx(pbuffer, 8584, 5511, 7800, 1, 28, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 8584, ket_width, ket_npgtos);
            }

            t2cfunc::transform<5, 6>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 5, 6, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecHI_hpp */
