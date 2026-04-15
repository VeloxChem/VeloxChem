#ifndef LocalCorePotentialGeom020GS_hpp
#define LocalCorePotentialGeom020GS_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "LocalCorePotentialPrimRecSS.hpp"
#include "LocalCorePotentialPrimRecSP.hpp"
#include "LocalCorePotentialPrimRecSD.hpp"
#include "LocalCorePotentialPrimRecPS.hpp"
#include "LocalCorePotentialPrimRecPP.hpp"
#include "LocalCorePotentialPrimRecPD.hpp"
#include "LocalCorePotentialPrimRecDS.hpp"
#include "LocalCorePotentialPrimRecDP.hpp"
#include "LocalCorePotentialPrimRecDD.hpp"
#include "LocalCorePotentialPrimRecFS.hpp"
#include "LocalCorePotentialPrimRecFP.hpp"
#include "LocalCorePotentialPrimRecFD.hpp"
#include "LocalCorePotentialPrimRecGS.hpp"
#include "LocalCorePotentialPrimRecGP.hpp"
#include "LocalCorePotentialPrimRecGD.hpp"
#include "LocalCorePotentialPrimRecHS.hpp"
#include "LocalCorePotentialPrimRecHP.hpp"
#include "LocalCorePotentialPrimRecIS.hpp"
#include "GeometricalDerivatives020ForGS.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes (G|U_L|S)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_local_core_potential_geom_020_gs(T& distributor,
                                      const CGtoBlock& bra_gto_block,
                                      const CGtoBlock& ket_gto_block,
                                      const CBaseCorePotential& ecp_potential,
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

    // intialize basic ECP data

    const auto ecp_nppt = ecp_potential.number_of_primitive_potentials();

    const auto ecp_exps = ecp_potential.get_exponents();

    const auto ecp_facts = ecp_potential.get_factors();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(15, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(552, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(90, 1);

    CSimdArray<double> sbuffer(54, 1);

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

        pfactors.load(ket_gto_exps, ket_range, 0, ket_npgtos);

        pfactors.load(ket_gto_norms, ket_range, 1, ket_npgtos);

        pfactors.replicate_points(ket_gto_coords, ket_range, 2, ket_npgtos);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        sbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        pbuffer.set_active_width(ket_width);

        // loop over contracted basis functions on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            cbuffer.zero();

            sbuffer.zero();

            const auto r_a = bra_gto_coords[j];

            for (size_t k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];

                const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];

                for (size_t l = 0; l < ecp_nppt; l++)
                {
                    const auto c_exp = ecp_exps[l];

                    const auto c_norm = ecp_facts[l];

                    t2cfunc::comp_coordinates_r(pfactors, 5, 2, r_a, a_exp, c_exp);

                    t2cfunc::comp_distances_ra(pfactors, 8 , 5, r_a);

                    t2cfunc::comp_distances_rb(pfactors, 11, 5, 2);

                    t2cfunc::comp_inverted_zeta(pfactors, 14, a_exp, c_exp);

                    t2lecp::comp_prim_local_core_potential_ss(pbuffer, 0, pfactors, 5, 14, r_a, a_exp, c_exp, a_norm, c_norm);

                    t2lecp::comp_prim_local_core_potential_sp(pbuffer, 1, 0, pfactors, 11, 14);

                    t2lecp::comp_prim_local_core_potential_sd(pbuffer, 4, 0, 1, pfactors, 11, 14);

                    t2lecp::comp_prim_local_core_potential_ps(pbuffer, 10, 0, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_pp(pbuffer, 13, 0, 1, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_pd(pbuffer, 22, 1, 4, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_ds(pbuffer, 40, 0, 10, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_dp(pbuffer, 46, 1, 10, 13, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_dd(pbuffer, 64, 4, 13, 22, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_fs(pbuffer, 100, 10, 40, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_fp(pbuffer, 110, 13, 40, 46, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_fd(pbuffer, 140, 22, 46, 64, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_gs(pbuffer, 200, 40, 100, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_gp(pbuffer, 215, 46, 100, 110, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_gd(pbuffer, 260, 64, 110, 140, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_hs(pbuffer, 350, 100, 200, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_hp(pbuffer, 371, 110, 200, 215, pfactors, 8, 14);

                    t2lecp::comp_prim_local_core_potential_is(pbuffer, 434, 200, 350, pfactors, 8, 14);

                    t2cgeom::comp_prim_op_geom_020_gs(pbuffer, 462, 40, 110, 200, 260, 371, 434, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 462, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 0>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 0, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2lecp namespace

#endif /* LocalCorePotentialGeom020GS_hpp */
