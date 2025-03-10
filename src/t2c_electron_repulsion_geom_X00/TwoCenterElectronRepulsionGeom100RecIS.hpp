#ifndef TwoCenterElectronRepulsionGeom100RecIS_hpp
#define TwoCenterElectronRepulsionGeom100RecIS_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecFS.hpp"
#include "TwoCenterElectronRepulsionPrimRecGS.hpp"
#include "TwoCenterElectronRepulsionPrimRecHS.hpp"
#include "TwoCenterElectronRepulsionPrimRecIS.hpp"
#include "TwoCenterElectronRepulsionPrimRecKS.hpp"
#include "GeometricalDerivatives1X0ForIY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)I|1/|r-r'||S)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_is(T& distributor,
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

    CSimdArray<double> factors(11, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(364, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(84, 1);

    CSimdArray<double> sbuffer(39, 1);

    // setup Boys function data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 7, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 7);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 7, 0, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 10, 1, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 13, 2, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 16, 3, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 19, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 22, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 25, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 28, 0, 1, 13, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 34, 1, 2, 16, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 40, 2, 3, 19, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 46, 3, 4, 22, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 52, 4, 5, 25, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 58, 7, 10, 28, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 68, 10, 13, 34, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 78, 13, 16, 40, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 88, 16, 19, 46, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 98, 19, 22, 52, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 108, 28, 34, 78, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 123, 34, 40, 88, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 138, 40, 46, 98, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hs(pbuffer, 153, 58, 68, 108, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hs(pbuffer, 174, 68, 78, 123, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hs(pbuffer, 195, 78, 88, 138, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_is(pbuffer, 216, 108, 123, 195, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ks(pbuffer, 244, 153, 174, 216, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_ix(pbuffer, 280, 153, 244, 1, 1, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 280, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 0>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 0, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecIS_hpp */
