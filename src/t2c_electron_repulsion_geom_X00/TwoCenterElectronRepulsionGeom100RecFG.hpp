#ifndef TwoCenterElectronRepulsionGeom100RecFG_hpp
#define TwoCenterElectronRepulsionGeom100RecFG_hpp

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
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "GeometricalDerivatives1X0ForFY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"
#include "BoysFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)F|1/|r-r'||G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_fg(T& distributor,
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

    CSimdArray<double> pbuffer(1740, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(450, 1);

    CSimdArray<double> sbuffer(189, 1);

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

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 0, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 9, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 12, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 15, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 18, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 21, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 24, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 27, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 30, 0, 1, 9, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 36, 1, 2, 12, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 42, 2, 3, 15, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 48, 3, 4, 18, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 54, 4, 5, 21, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 60, 5, 6, 24, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 66, 6, 7, 27, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 72, 9, 12, 42, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 82, 12, 15, 48, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 92, 15, 18, 54, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 102, 18, 21, 60, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 112, 21, 24, 66, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 122, 30, 36, 72, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 137, 36, 42, 82, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 152, 42, 48, 92, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 167, 48, 54, 102, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 182, 54, 60, 112, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 197, 4, 18, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 206, 18, 54, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 224, 42, 82, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 254, 48, 92, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 284, 54, 102, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 314, 82, 152, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 359, 92, 167, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 404, 102, 182, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 449, 42, 48, 197, 206, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 485, 82, 92, 206, 284, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 545, 122, 137, 224, 314, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 635, 137, 152, 254, 359, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 725, 152, 167, 284, 404, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 815, 224, 254, 449, 485, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 915, 314, 359, 485, 725, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 1065, 545, 635, 815, 915, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_fx(pbuffer, 1290, 545, 1065, 1, 15, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 1290, ket_width, ket_npgtos);
            }

            t2cfunc::transform<3, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecFG_hpp */
