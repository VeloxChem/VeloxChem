#ifndef TwoCenterElectronRepulsionGeom100RecDH_hpp
#define TwoCenterElectronRepulsionGeom100RecDH_hpp

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
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "GeometricalDerivatives1X0ForDY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)D|1/|r-r'||H)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_dh(T& distributor,
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

    CSimdArray<double> pbuffer(1306, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(378, 1);

    CSimdArray<double> sbuffer(165, 1);

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

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 29, 0, 1, 11, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 35, 1, 2, 14, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 41, 2, 3, 17, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 47, 3, 4, 20, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 53, 4, 5, 23, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 59, 5, 6, 26, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 65, 8, 11, 35, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 75, 11, 14, 41, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 85, 14, 17, 47, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 95, 17, 20, 53, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 105, 20, 23, 59, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 115, 29, 35, 75, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 130, 35, 41, 85, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 145, 41, 47, 95, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 160, 47, 53, 105, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 175, 65, 75, 130, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 196, 75, 85, 145, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 217, 85, 95, 160, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 238, 41, 85, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 268, 85, 145, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 313, 115, 175, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 376, 130, 196, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 439, 145, 217, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 502, 115, 130, 238, 268, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 592, 175, 196, 268, 439, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 718, 313, 376, 502, 592, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_dx(pbuffer, 928, 313, 718, 1, 21, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 928, ket_width, ket_npgtos);
            }

            t2cfunc::transform<2, 5>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 5, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecDH_hpp */
