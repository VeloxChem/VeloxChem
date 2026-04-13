#ifndef ProjectedCorePotentialGeom010GGForG_hpp
#define ProjectedCorePotentialGeom010GGForG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "SimdArray.hpp"
#include "ProjectedCorePotentialPrimRecDDForD.hpp"
#include "ProjectedCorePotentialPrimRecDDForF.hpp"
#include "ProjectedCorePotentialPrimRecDDForP.hpp"
#include "ProjectedCorePotentialPrimRecDDForS.hpp"
#include "ProjectedCorePotentialPrimRecDFForD.hpp"
#include "ProjectedCorePotentialPrimRecDFForF.hpp"
#include "ProjectedCorePotentialPrimRecDFForG.hpp"
#include "ProjectedCorePotentialPrimRecDFForP.hpp"
#include "ProjectedCorePotentialPrimRecDFForS.hpp"
#include "ProjectedCorePotentialPrimRecDGForD.hpp"
#include "ProjectedCorePotentialPrimRecDGForF.hpp"
#include "ProjectedCorePotentialPrimRecDGForG.hpp"
#include "ProjectedCorePotentialPrimRecDGForP.hpp"
#include "ProjectedCorePotentialPrimRecDGForS.hpp"
#include "ProjectedCorePotentialPrimRecDHForD.hpp"
#include "ProjectedCorePotentialPrimRecDHForF.hpp"
#include "ProjectedCorePotentialPrimRecDHForG.hpp"
#include "ProjectedCorePotentialPrimRecDHForP.hpp"
#include "ProjectedCorePotentialPrimRecDHForS.hpp"
#include "ProjectedCorePotentialPrimRecDPForD.hpp"
#include "ProjectedCorePotentialPrimRecDPForP.hpp"
#include "ProjectedCorePotentialPrimRecDPForS.hpp"
#include "ProjectedCorePotentialPrimRecFDForD.hpp"
#include "ProjectedCorePotentialPrimRecFDForF.hpp"
#include "ProjectedCorePotentialPrimRecFDForP.hpp"
#include "ProjectedCorePotentialPrimRecFDForS.hpp"
#include "ProjectedCorePotentialPrimRecFFForD.hpp"
#include "ProjectedCorePotentialPrimRecFFForF.hpp"
#include "ProjectedCorePotentialPrimRecFFForG.hpp"
#include "ProjectedCorePotentialPrimRecFFForP.hpp"
#include "ProjectedCorePotentialPrimRecFFForS.hpp"
#include "ProjectedCorePotentialPrimRecFGForD.hpp"
#include "ProjectedCorePotentialPrimRecFGForF.hpp"
#include "ProjectedCorePotentialPrimRecFGForG.hpp"
#include "ProjectedCorePotentialPrimRecFGForP.hpp"
#include "ProjectedCorePotentialPrimRecFGForS.hpp"
#include "ProjectedCorePotentialPrimRecFHForD.hpp"
#include "ProjectedCorePotentialPrimRecFHForF.hpp"
#include "ProjectedCorePotentialPrimRecFHForG.hpp"
#include "ProjectedCorePotentialPrimRecFHForP.hpp"
#include "ProjectedCorePotentialPrimRecFHForS.hpp"
#include "ProjectedCorePotentialPrimRecGFForF.hpp"
#include "ProjectedCorePotentialPrimRecGFForG.hpp"
#include "ProjectedCorePotentialPrimRecGFForP.hpp"
#include "ProjectedCorePotentialPrimRecGGForD.hpp"
#include "ProjectedCorePotentialPrimRecGGForF.hpp"
#include "ProjectedCorePotentialPrimRecGGForG.hpp"
#include "ProjectedCorePotentialPrimRecGGForP.hpp"
#include "ProjectedCorePotentialPrimRecGGForS.hpp"
#include "ProjectedCorePotentialPrimRecGHForG.hpp"
#include "ProjectedCorePotentialPrimRecHGForG.hpp"
#include "ProjectedCorePotentialPrimRecPDForD.hpp"
#include "ProjectedCorePotentialPrimRecPDForF.hpp"
#include "ProjectedCorePotentialPrimRecPDForP.hpp"
#include "ProjectedCorePotentialPrimRecPDForS.hpp"
#include "ProjectedCorePotentialPrimRecPFForD.hpp"
#include "ProjectedCorePotentialPrimRecPFForF.hpp"
#include "ProjectedCorePotentialPrimRecPFForG.hpp"
#include "ProjectedCorePotentialPrimRecPFForP.hpp"
#include "ProjectedCorePotentialPrimRecPFForS.hpp"
#include "ProjectedCorePotentialPrimRecPGForD.hpp"
#include "ProjectedCorePotentialPrimRecPGForF.hpp"
#include "ProjectedCorePotentialPrimRecPGForG.hpp"
#include "ProjectedCorePotentialPrimRecPGForP.hpp"
#include "ProjectedCorePotentialPrimRecPGForS.hpp"
#include "ProjectedCorePotentialPrimRecPHForD.hpp"
#include "ProjectedCorePotentialPrimRecPHForF.hpp"
#include "ProjectedCorePotentialPrimRecPHForG.hpp"
#include "ProjectedCorePotentialPrimRecPHForP.hpp"
#include "ProjectedCorePotentialPrimRecPHForS.hpp"
#include "ProjectedCorePotentialPrimRecPPForD.hpp"
#include "ProjectedCorePotentialPrimRecPPForP.hpp"
#include "ProjectedCorePotentialPrimRecPPForS.hpp"
#include "ProjectedCorePotentialPrimRecPSForP.hpp"
#include "ProjectedCorePotentialPrimRecPSForS.hpp"
#include "ProjectedCorePotentialPrimRecSDForD.hpp"
#include "ProjectedCorePotentialPrimRecSDForF.hpp"
#include "ProjectedCorePotentialPrimRecSDForG.hpp"
#include "ProjectedCorePotentialPrimRecSDForP.hpp"
#include "ProjectedCorePotentialPrimRecSDForS.hpp"
#include "ProjectedCorePotentialPrimRecSFForD.hpp"
#include "ProjectedCorePotentialPrimRecSFForF.hpp"
#include "ProjectedCorePotentialPrimRecSFForG.hpp"
#include "ProjectedCorePotentialPrimRecSFForP.hpp"
#include "ProjectedCorePotentialPrimRecSFForS.hpp"
#include "ProjectedCorePotentialPrimRecSGForD.hpp"
#include "ProjectedCorePotentialPrimRecSGForF.hpp"
#include "ProjectedCorePotentialPrimRecSGForG.hpp"
#include "ProjectedCorePotentialPrimRecSGForP.hpp"
#include "ProjectedCorePotentialPrimRecSGForS.hpp"
#include "ProjectedCorePotentialPrimRecSHForD.hpp"
#include "ProjectedCorePotentialPrimRecSHForF.hpp"
#include "ProjectedCorePotentialPrimRecSHForG.hpp"
#include "ProjectedCorePotentialPrimRecSHForP.hpp"
#include "ProjectedCorePotentialPrimRecSHForS.hpp"
#include "ProjectedCorePotentialPrimRecSPForD.hpp"
#include "ProjectedCorePotentialPrimRecSPForF.hpp"
#include "ProjectedCorePotentialPrimRecSPForG.hpp"
#include "ProjectedCorePotentialPrimRecSPForP.hpp"
#include "ProjectedCorePotentialPrimRecSPForS.hpp"
#include "ProjectedCorePotentialPrimRecSS.hpp"
#include "GeometricalDerivatives010ForGG.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes (G|U_l|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_projected_core_potential_geom_010_gg_for_g(T& distributor,
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

    CSimdArray<double> pfactors(9, ket_npgtos);

    // allpcate I_n and L_n values

    CSimdArray<double> i_values(14, ket_npgtos);

    CSimdArray<double> l_values(5, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(26122, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(675, 1);

    CSimdArray<double> sbuffer(243, 1);

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

        i_values.set_active_width(ket_width);

        l_values.set_active_width(ket_width);

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

                    t2cfunc::comp_coordinates_norm(pfactors, 5, 2);

                    t2cfunc::comp_legendre_args(pfactors, 6, 2, 5, r_a);

                    t2cfunc::comp_gamma_factors(pfactors, 7, 5, r_a, a_exp, c_exp);

                    t2cfunc::comp_bessel_args(pfactors, 8, 5, r_a, a_exp, c_exp);

                    t2cfunc::comp_i_vals(i_values, 13, pfactors, 8);

                    t2cfunc::comp_l_vals(l_values, 4, pfactors, 8, 6);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 0, pbuffer, 0, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 1, pbuffer, 2610, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 1, pbuffer, 2611, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 2, pbuffer, 4064, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 2, pbuffer, 4065, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 2, pbuffer, 4066, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 3, pbuffer, 4980, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 3, pbuffer, 4981, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 3, pbuffer, 4982, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 3, pbuffer, 4983, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 4, pbuffer, 5473, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 4, pbuffer, 5474, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 4, pbuffer, 5475, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 4, pbuffer, 5476, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 4, pbuffer, 5477, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 5, pbuffer, 5703, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 5, pbuffer, 5704, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 5, pbuffer, 5705, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 5, pbuffer, 5706, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 6, pbuffer, 5738, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 6, pbuffer, 5739, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 0, pbuffer, 5743, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 1, pbuffer, 5778, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 1, pbuffer, 5779, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 1, pbuffer, 5780, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 2, pbuffer, 5868, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 2, pbuffer, 5869, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 2, pbuffer, 5870, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 2, pbuffer, 5871, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 3, pbuffer, 5968, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 3, pbuffer, 5969, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 3, pbuffer, 5970, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 3, pbuffer, 5971, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 3, pbuffer, 5972, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 4, pbuffer, 6072, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 4, pbuffer, 6073, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 4, pbuffer, 6074, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 4, pbuffer, 6075, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 4, pbuffer, 6076, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 5, pbuffer, 6142, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 5, pbuffer, 6143, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 5, pbuffer, 6144, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 6, pbuffer, 6157, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 0, pbuffer, 6158, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 1, pbuffer, 6178, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 1, pbuffer, 6179, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 1, pbuffer, 6180, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 2, pbuffer, 6228, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 2, pbuffer, 6229, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 2, pbuffer, 6230, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 2, pbuffer, 6231, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 2, pbuffer, 6232, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 3, pbuffer, 6302, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 3, pbuffer, 6303, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 3, pbuffer, 6304, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 3, pbuffer, 6305, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 3, pbuffer, 6306, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 4, pbuffer, 6357, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 4, pbuffer, 6358, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 4, pbuffer, 6359, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 4, pbuffer, 6360, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 5, pbuffer, 6392, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 5, pbuffer, 6393, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 3, 0, pbuffer, 6397, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 1, pbuffer, 6407, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 1, pbuffer, 6408, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 3, 1, pbuffer, 6409, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 2, pbuffer, 6431, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 2, pbuffer, 6432, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 2, pbuffer, 6433, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 2, pbuffer, 6434, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 3, 2, pbuffer, 6435, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 3, pbuffer, 6466, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 3, pbuffer, 6467, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 3, pbuffer, 6468, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 3, pbuffer, 6469, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 4, pbuffer, 6491, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 4, pbuffer, 6492, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 4, pbuffer, 6493, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 5, pbuffer, 6506, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 4, 0, pbuffer, 6507, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 1, pbuffer, 6511, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 1, pbuffer, 6512, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 4, 1, pbuffer, 6513, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 2, pbuffer, 6520, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 2, pbuffer, 6521, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 2, pbuffer, 6522, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 2, pbuffer, 6523, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 3, pbuffer, 6533, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 3, pbuffer, 6534, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 3, pbuffer, 6535, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 4, pbuffer, 6542, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 4, pbuffer, 6543, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 5, 0, pbuffer, 6547, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 5, 1, pbuffer, 6548, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 5, 1, pbuffer, 6549, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 5, 2, pbuffer, 6550, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 5, 2, pbuffer, 6551, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 5, 2, pbuffer, 6552, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 5, 3, pbuffer, 6553, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 5, 3, pbuffer, 6554, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 5, 4, pbuffer, 6555, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 0, pbuffer, 6556, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 1, pbuffer, 7711, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 1, pbuffer, 7712, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 1, pbuffer, 7713, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 2, pbuffer, 9675, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 2, pbuffer, 9676, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 2, pbuffer, 9677, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 2, pbuffer, 9678, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 3, pbuffer, 10911, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 3, pbuffer, 10912, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 3, pbuffer, 10913, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 3, pbuffer, 10914, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 3, pbuffer, 10915, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 4, pbuffer, 11541, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 4, pbuffer, 11542, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 4, pbuffer, 11543, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 4, pbuffer, 11544, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 4, pbuffer, 11545, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 5, pbuffer, 11611, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 5, pbuffer, 11612, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 5, pbuffer, 11613, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 6, pbuffer, 11626, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 0, pbuffer, 11627, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 1, pbuffer, 11662, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 1, pbuffer, 11663, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 1, pbuffer, 11664, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 1, pbuffer, 11665, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 2, pbuffer, 13171, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 2, pbuffer, 13172, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 2, pbuffer, 13173, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 2, pbuffer, 13174, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 2, pbuffer, 13175, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 3, pbuffer, 14132, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 3, pbuffer, 14133, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 3, pbuffer, 14134, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 3, pbuffer, 14135, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 3, pbuffer, 14136, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 4, pbuffer, 14236, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 4, pbuffer, 14237, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 4, pbuffer, 14238, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 4, pbuffer, 14239, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 5, pbuffer, 14271, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 5, pbuffer, 14272, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 0, pbuffer, 14276, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 1, pbuffer, 14296, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 1, pbuffer, 14297, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 1, pbuffer, 14298, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 1, pbuffer, 14299, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 1, pbuffer, 14300, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 2, pbuffer, 14416, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 2, pbuffer, 14417, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 2, pbuffer, 14418, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 2, pbuffer, 14419, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 2, pbuffer, 14420, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 3, pbuffer, 14505, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 3, pbuffer, 14506, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 3, pbuffer, 14507, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 3, pbuffer, 14508, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 3, pbuffer, 14509, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 4, pbuffer, 14560, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 4, pbuffer, 14561, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 4, pbuffer, 14562, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 5, pbuffer, 14575, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 3, 0, pbuffer, 14576, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 1, pbuffer, 14586, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 1, pbuffer, 14587, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 1, pbuffer, 14588, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 1, pbuffer, 14589, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 3, 1, pbuffer, 14590, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 2, pbuffer, 14650, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 2, pbuffer, 14651, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 2, pbuffer, 14652, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 2, pbuffer, 14653, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 3, 2, pbuffer, 14654, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 3, pbuffer, 14695, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 3, pbuffer, 14696, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 3, pbuffer, 14697, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 3, pbuffer, 14698, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 4, pbuffer, 14720, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 4, pbuffer, 14721, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 4, 0, pbuffer, 14725, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 1, pbuffer, 14729, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 1, pbuffer, 14730, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 1, pbuffer, 14731, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 4, 1, pbuffer, 14732, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 4, 1, pbuffer, 14733, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 2, pbuffer, 14758, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 2, pbuffer, 14759, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 2, pbuffer, 14760, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 4, 2, pbuffer, 14761, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 3, pbuffer, 14777, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 3, pbuffer, 14778, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 3, pbuffer, 14779, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 4, pbuffer, 14786, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 5, 0, pbuffer, 14787, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 1, pbuffer, 14788, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 1, pbuffer, 14789, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 5, 1, pbuffer, 14790, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 5, 1, pbuffer, 14791, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 2, pbuffer, 14798, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 2, pbuffer, 14799, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 5, 2, pbuffer, 14800, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 3, pbuffer, 14804, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 3, pbuffer, 14805, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 6, 1, pbuffer, 14806, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 6, 1, pbuffer, 14807, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 6, 2, pbuffer, 14808, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 0, pbuffer, 14809, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 1, pbuffer, 15429, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 1, pbuffer, 15430, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 1, pbuffer, 15431, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 2, pbuffer, 16446, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 2, pbuffer, 16447, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 2, pbuffer, 16448, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 2, pbuffer, 16449, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 2, pbuffer, 16450, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 3, pbuffer, 17632, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 3, pbuffer, 17633, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 3, pbuffer, 17634, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 3, pbuffer, 17635, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 3, pbuffer, 17636, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 4, pbuffer, 17850, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 4, pbuffer, 17851, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 4, pbuffer, 17852, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 4, pbuffer, 17853, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 5, pbuffer, 17885, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 5, pbuffer, 17886, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 0, pbuffer, 17890, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 1, pbuffer, 17925, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 1, pbuffer, 17926, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 1, pbuffer, 17927, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 1, pbuffer, 17928, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 1, pbuffer, 17929, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 2, pbuffer, 19944, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 2, pbuffer, 19945, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 2, pbuffer, 19946, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 2, pbuffer, 19947, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 2, pbuffer, 19948, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 3, pbuffer, 20451, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 3, pbuffer, 20452, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 3, pbuffer, 20453, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 3, pbuffer, 20454, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 3, pbuffer, 20455, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 4, pbuffer, 20521, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 4, pbuffer, 20522, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 4, pbuffer, 20523, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 5, pbuffer, 20536, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 2, 0, pbuffer, 20537, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 1, pbuffer, 20557, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 1, pbuffer, 20558, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 1, pbuffer, 20559, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 1, pbuffer, 20560, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 2, 1, pbuffer, 20561, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 2, pbuffer, 20677, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 2, pbuffer, 20678, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 2, pbuffer, 20679, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 2, pbuffer, 20680, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 2, 2, pbuffer, 20681, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 3, pbuffer, 20766, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 3, pbuffer, 20767, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 3, pbuffer, 20768, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 3, pbuffer, 20769, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 4, pbuffer, 20801, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 4, pbuffer, 20802, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 3, 0, pbuffer, 20806, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 1, pbuffer, 20816, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 1, pbuffer, 20817, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 1, pbuffer, 20818, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 1, pbuffer, 20819, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 3, 1, pbuffer, 20820, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 2, pbuffer, 20880, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 2, pbuffer, 20881, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 2, pbuffer, 20882, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 2, pbuffer, 20883, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 3, 2, pbuffer, 20884, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 3, pbuffer, 20925, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 3, pbuffer, 20926, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 3, pbuffer, 20927, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 4, pbuffer, 20940, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 4, 0, pbuffer, 20941, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 1, pbuffer, 20945, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 1, pbuffer, 20946, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 4, 1, pbuffer, 20947, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 4, 1, pbuffer, 20948, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 4, 1, pbuffer, 20949, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 2, pbuffer, 20974, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 2, pbuffer, 20975, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 4, 2, pbuffer, 20976, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 4, 2, pbuffer, 20977, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 3, pbuffer, 20993, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 3, pbuffer, 20994, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 5, 0, pbuffer, 20998, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 1, pbuffer, 20999, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 5, 1, pbuffer, 21000, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 5, 1, pbuffer, 21001, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 5, 1, pbuffer, 21002, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 2, pbuffer, 21009, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 5, 2, pbuffer, 21010, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 5, 2, pbuffer, 21011, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 3, pbuffer, 21015, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 6, 1, pbuffer, 21016, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 6, 1, pbuffer, 21017, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 6, 2, pbuffer, 21018, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 0, 0, pbuffer, 21019, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 1, pbuffer, 21303, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 1, pbuffer, 21304, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 0, 1, pbuffer, 21305, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 2, pbuffer, 21753, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 2, pbuffer, 21754, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 2, pbuffer, 21755, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 2, pbuffer, 21756, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 0, 2, pbuffer, 21757, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 3, pbuffer, 22233, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 3, pbuffer, 22234, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 3, pbuffer, 22235, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 3, pbuffer, 22236, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 4, pbuffer, 22302, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 4, pbuffer, 22303, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 4, pbuffer, 22304, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 5, pbuffer, 22317, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 1, 0, pbuffer, 22318, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 1, pbuffer, 22353, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 1, pbuffer, 22354, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 1, pbuffer, 22355, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 1, pbuffer, 22356, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 1, 1, pbuffer, 22357, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 2, pbuffer, 23425, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 2, pbuffer, 23426, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 2, pbuffer, 23427, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 2, pbuffer, 23428, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 1, 2, pbuffer, 23429, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 3, pbuffer, 23643, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 3, pbuffer, 23644, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 3, pbuffer, 23645, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 3, pbuffer, 23646, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 4, pbuffer, 23678, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 4, pbuffer, 23679, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 2, 0, pbuffer, 23683, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 1, pbuffer, 23703, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 1, pbuffer, 23704, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 1, pbuffer, 23705, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 2, 1, pbuffer, 23706, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 2, 1, pbuffer, 23707, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 2, pbuffer, 23823, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 2, pbuffer, 23824, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 2, pbuffer, 23825, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 2, 2, pbuffer, 23826, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 2, 2, pbuffer, 23827, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 3, pbuffer, 23893, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 3, pbuffer, 23894, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 3, pbuffer, 23895, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 4, pbuffer, 23908, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 3, 0, pbuffer, 23909, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 1, pbuffer, 23919, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 1, pbuffer, 23920, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 3, 1, pbuffer, 23921, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 3, 1, pbuffer, 23922, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 3, 1, pbuffer, 23923, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 2, pbuffer, 23983, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 2, pbuffer, 23984, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 3, 2, pbuffer, 23985, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 3, 2, pbuffer, 23986, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 3, pbuffer, 24018, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 3, pbuffer, 24019, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 4, 0, pbuffer, 24023, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 1, pbuffer, 24027, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 4, 1, pbuffer, 24028, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 4, 1, pbuffer, 24029, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 4, 1, pbuffer, 24030, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 4, 1, pbuffer, 24031, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 2, pbuffer, 24056, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 4, 2, pbuffer, 24057, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 4, 2, pbuffer, 24058, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 3, pbuffer, 24071, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 5, 0, pbuffer, 24072, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 5, 1, pbuffer, 24073, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 5, 1, pbuffer, 24074, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 5, 1, pbuffer, 24075, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 5, 1, pbuffer, 24076, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 5, 2, pbuffer, 24083, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 5, 2, pbuffer, 24084, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 6, 1, pbuffer, 24088, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 6, 1, pbuffer, 24089, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 6, 2, pbuffer, 24090, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 0, 0, pbuffer, 24091, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 0, 1, pbuffer, 24192, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 0, 1, pbuffer, 24193, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 0, 1, pbuffer, 24194, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 2, pbuffer, 24348, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 0, 2, pbuffer, 24349, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 0, 2, pbuffer, 24350, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 0, 2, pbuffer, 24351, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 3, pbuffer, 24514, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 0, 3, pbuffer, 24515, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 0, 3, pbuffer, 24516, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 4, pbuffer, 24548, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 0, 4, pbuffer, 24549, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 1, 0, pbuffer, 24553, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 1, pbuffer, 24588, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 1, 1, pbuffer, 24589, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 1, 1, pbuffer, 24590, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 1, 1, pbuffer, 24591, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 1, 1, pbuffer, 24592, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 2, pbuffer, 25068, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 1, 2, pbuffer, 25069, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 1, 2, pbuffer, 25070, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 1, 2, pbuffer, 25071, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 3, pbuffer, 25137, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 1, 3, pbuffer, 25138, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 1, 3, pbuffer, 25139, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 4, pbuffer, 25152, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 2, 0, pbuffer, 25153, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 1, pbuffer, 25173, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 2, 1, pbuffer, 25174, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 2, 1, pbuffer, 25175, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 2, 1, pbuffer, 25176, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 2, 1, pbuffer, 25177, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 2, pbuffer, 25277, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 2, 2, pbuffer, 25278, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 2, 2, pbuffer, 25279, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 2, 2, pbuffer, 25280, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 3, pbuffer, 25312, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 2, 3, pbuffer, 25313, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 3, 0, pbuffer, 25317, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 1, pbuffer, 25327, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 3, 1, pbuffer, 25328, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 3, 1, pbuffer, 25329, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 3, 1, pbuffer, 25330, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 3, 1, pbuffer, 25331, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 2, pbuffer, 25382, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 3, 2, pbuffer, 25383, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 3, 2, pbuffer, 25384, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 3, pbuffer, 25397, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 4, 0, pbuffer, 25398, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 4, 1, pbuffer, 25402, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 4, 1, pbuffer, 25403, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 4, 1, pbuffer, 25404, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 4, 1, pbuffer, 25405, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 4, 2, pbuffer, 25427, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 4, 2, pbuffer, 25428, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 5, 0, pbuffer, 25432, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 5, 1, pbuffer, 25433, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 5, 1, pbuffer, 25434, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 5, 1, pbuffer, 25435, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 5, 2, pbuffer, 25442, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 6, 1, pbuffer, 25443, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 6, 1, pbuffer, 25444, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 5, 0, 0, pbuffer, 25445, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 0, 1, pbuffer, 25480, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 5, 0, 1, pbuffer, 25481, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 0, 2, pbuffer, 25535, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 0, 2, pbuffer, 25536, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 0, 2, pbuffer, 25537, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 0, 3, pbuffer, 25600, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 0, 3, pbuffer, 25601, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 0, 4, pbuffer, 25614, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 5, 1, 0, pbuffer, 25615, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 1, pbuffer, 25635, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 1, 1, pbuffer, 25636, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 1, 1, pbuffer, 25637, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 5, 1, 1, pbuffer, 25638, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 2, pbuffer, 25801, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 1, 2, pbuffer, 25802, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 1, 2, pbuffer, 25803, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 3, pbuffer, 25835, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 1, 3, pbuffer, 25836, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 5, 2, 0, pbuffer, 25840, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 2, 1, pbuffer, 25850, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 2, 1, pbuffer, 25851, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 2, 1, pbuffer, 25852, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 5, 2, 1, pbuffer, 25853, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 2, 2, pbuffer, 25919, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 2, 2, pbuffer, 25920, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 2, 2, pbuffer, 25921, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 2, 3, pbuffer, 25934, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 5, 3, 0, pbuffer, 25935, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 3, 1, pbuffer, 25939, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 3, 1, pbuffer, 25940, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 3, 1, pbuffer, 25941, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 5, 3, 1, pbuffer, 25942, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 3, 2, pbuffer, 25974, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 3, 2, pbuffer, 25975, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 5, 4, 0, pbuffer, 25979, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 4, 1, pbuffer, 25980, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 4, 1, pbuffer, 25981, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 5, 4, 1, pbuffer, 25982, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 4, 2, pbuffer, 25995, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 5, 1, pbuffer, 25996, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 5, 1, pbuffer, 25997, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 6, 1, pbuffer, 26001, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 1, 1, pbuffer, 26002, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 6, 1, 1, pbuffer, 26003, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 1, 2, pbuffer, 26057, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 2, 1, pbuffer, 26067, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 6, 2, 1, pbuffer, 26068, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 2, 2, pbuffer, 26097, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 3, 1, pbuffer, 26101, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 6, 3, 1, pbuffer, 26102, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 3, 2, pbuffer, 26115, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 4, 1, pbuffer, 26116, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 6, 4, 1, pbuffer, 26117, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 5, 1, pbuffer, 26121, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5478, 5473, 6072, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5707, 5703, 6142, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5740, 5738, 6157, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5973, 5968, 6302, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6077, 6072, 6357, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6145, 6142, 6392, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6233, 6228, 6431, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6307, 6302, 6466, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6361, 6357, 6491, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6394, 6392, 6506, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6436, 6431, 6520, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6470, 6466, 6533, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6494, 6491, 6542, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6524, 6520, 6550, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6536, 6533, 6553, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6544, 6542, 6555, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 10916, 10911, 14132, 1, 5473, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11546, 11541, 14236, 1, 5703, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11614, 11611, 14271, 1, 5738, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13176, 13171, 14416, 1, 5968, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14137, 14132, 14505, 1, 6072, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14240, 14236, 14560, 1, 6142, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14273, 14271, 14575, 1, 6157, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14301, 14296, 14586, 1, 6228, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14421, 14416, 14650, 1, 6302, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14510, 14505, 14695, 1, 6357, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14563, 14560, 14720, 1, 6392, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14591, 14586, 14729, 1, 6431, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14655, 14650, 14758, 1, 6466, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14699, 14695, 14777, 1, 6491, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14722, 14720, 14786, 1, 6506, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14734, 14729, 14788, 1, 6520, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14762, 14758, 14798, 1, 6533, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14780, 14777, 14804, 1, 6542, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14792, 14788, 14806, 1, 6550, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 14801, 14798, 14808, 1, 6553, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 16451, 16446, 19944, 2, 10911, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 17637, 17632, 20451, 2, 11541, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 17854, 17850, 20521, 2, 11611, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 17887, 17885, 20536, 2, 11626, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 17930, 17925, 20557, 2, 13171, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 19949, 19944, 20677, 2, 14132, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20456, 20451, 20766, 2, 14236, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20524, 20521, 20801, 2, 14271, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20562, 20557, 20816, 2, 14416, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20682, 20677, 20880, 2, 14505, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20770, 20766, 20925, 2, 14560, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20803, 20801, 20940, 2, 14575, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20821, 20816, 20945, 2, 14650, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20885, 20880, 20974, 2, 14695, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20928, 20925, 20993, 2, 14720, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20950, 20945, 20999, 2, 14758, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20978, 20974, 21009, 2, 14777, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 20995, 20993, 21015, 2, 14786, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 21003, 20999, 21016, 2, 14798, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 21012, 21009, 21018, 2, 14804, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 21758, 21753, 23425, 3, 17632, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 22237, 22233, 23643, 3, 17850, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 22305, 22302, 23678, 3, 17885, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 22358, 22353, 23703, 3, 19944, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23430, 23425, 23823, 3, 20451, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23647, 23643, 23893, 3, 20521, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23680, 23678, 23908, 3, 20536, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23708, 23703, 23919, 3, 20677, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23828, 23823, 23983, 3, 20766, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23896, 23893, 24018, 3, 20801, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23924, 23919, 24027, 3, 20880, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 23987, 23983, 24056, 3, 20925, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24020, 24018, 24071, 3, 20940, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24032, 24027, 24073, 3, 20974, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24059, 24056, 24083, 3, 20993, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24077, 24073, 24088, 3, 21009, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24085, 24083, 24090, 3, 21015, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24352, 24348, 25068, 4, 22233, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24517, 24514, 25137, 4, 22302, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24550, 24548, 25152, 4, 22317, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 24593, 24588, 25173, 4, 23425, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25072, 25068, 25277, 4, 23643, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25140, 25137, 25312, 4, 23678, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25178, 25173, 25327, 4, 23823, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25281, 25277, 25382, 4, 23893, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25314, 25312, 25397, 4, 23908, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25332, 25327, 25402, 4, 23983, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25385, 25382, 25427, 4, 24018, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25406, 25402, 25433, 4, 24056, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25429, 25427, 25442, 4, 24071, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25436, 25433, 25443, 4, 24083, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25538, 25535, 25801, 5, 24514, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25602, 25600, 25835, 5, 24548, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25639, 25635, 25850, 5, 25068, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25804, 25801, 25919, 5, 25137, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25837, 25835, 25934, 5, 25152, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25854, 25850, 25939, 5, 25277, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25922, 25919, 25974, 5, 25312, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25943, 25939, 25980, 5, 25382, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25976, 25974, 25995, 5, 25397, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25983, 25980, 25996, 5, 25427, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 25998, 25996, 26001, 5, 25442, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 26004, 26002, 26067, 6, 25801, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 26058, 26057, 26097, 6, 25835, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 26069, 26067, 26101, 6, 25919, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 26098, 26097, 26115, 6, 25934, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 26103, 26101, 26116, 6, 25974, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 26118, 26116, 26121, 6, 25995, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4984, 4980, 5473, 5969, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5481, 5474, 5703, 6073, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5710, 5704, 5738, 6143, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5872, 5868, 5968, 6229, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5976, 5969, 6072, 6303, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6080, 6073, 6142, 6358, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6148, 6143, 6157, 6393, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6236, 6229, 6302, 6432, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6310, 6303, 6357, 6467, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6364, 6358, 6392, 6492, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6439, 6432, 6466, 6521, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6473, 6467, 6491, 6534, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6497, 6492, 6506, 6543, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6527, 6521, 6533, 6551, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6539, 6534, 6542, 6554, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9679, 9675, 10911, 13172, 1, 4980, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 10919, 10912, 11541, 14133, 1, 5474, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11549, 11542, 11611, 14237, 1, 5704, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11617, 11612, 11626, 14272, 1, 5739, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11666, 11662, 13171, 14297, 1, 5868, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13179, 13172, 14132, 14417, 1, 5969, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14140, 14133, 14236, 14506, 1, 6073, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14243, 14237, 14271, 14561, 1, 6143, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14304, 14297, 14416, 14587, 1, 6229, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14424, 14417, 14505, 14651, 1, 6303, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14513, 14506, 14560, 14696, 1, 6358, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14566, 14561, 14575, 14721, 1, 6393, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14594, 14587, 14650, 14730, 1, 6432, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14658, 14651, 14695, 14759, 1, 6467, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14702, 14696, 14720, 14778, 1, 6492, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14737, 14730, 14758, 14789, 1, 6521, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14765, 14759, 14777, 14799, 1, 6534, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14783, 14778, 14786, 14805, 1, 6543, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 14795, 14789, 14798, 14807, 1, 6551, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 16454, 16447, 17632, 19945, 2, 10912, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 17640, 17633, 17850, 20452, 2, 11542, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 17857, 17851, 17885, 20522, 2, 11612, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 17933, 17926, 19944, 20558, 2, 13172, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 19952, 19945, 20451, 20678, 2, 14133, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20459, 20452, 20521, 20767, 2, 14237, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20527, 20522, 20536, 20802, 2, 14272, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20565, 20558, 20677, 20817, 2, 14417, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20685, 20678, 20766, 20881, 2, 14506, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20773, 20767, 20801, 20926, 2, 14561, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20824, 20817, 20880, 20946, 2, 14651, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20888, 20881, 20925, 20975, 2, 14696, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20931, 20926, 20940, 20994, 2, 14721, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20953, 20946, 20974, 21000, 2, 14759, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 20981, 20975, 20993, 21010, 2, 14778, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 21006, 21000, 21009, 21017, 2, 14799, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 21761, 21754, 22233, 23426, 3, 17633, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 22240, 22234, 22302, 23644, 3, 17851, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 22308, 22303, 22317, 23679, 3, 17886, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 22361, 22354, 23425, 23704, 3, 19945, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 23433, 23426, 23643, 23824, 3, 20452, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 23650, 23644, 23678, 23894, 3, 20522, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 23711, 23704, 23823, 23920, 3, 20678, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 23831, 23824, 23893, 23984, 3, 20767, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 23899, 23894, 23908, 24019, 3, 20802, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 23927, 23920, 23983, 24028, 3, 20881, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 23990, 23984, 24018, 24057, 3, 20926, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 24035, 24028, 24056, 24074, 3, 20975, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 24062, 24057, 24071, 24084, 3, 20994, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 24080, 24074, 24083, 24089, 3, 21010, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 24355, 24349, 24514, 25069, 4, 22234, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 24520, 24515, 24548, 25138, 4, 22303, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 24596, 24589, 25068, 25174, 4, 23426, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25075, 25069, 25137, 25278, 4, 23644, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25143, 25138, 25152, 25313, 4, 23679, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25181, 25174, 25277, 25328, 4, 23824, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25284, 25278, 25312, 25383, 4, 23894, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25335, 25328, 25382, 25403, 4, 23984, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25388, 25383, 25397, 25428, 4, 24019, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25409, 25403, 25427, 25434, 4, 24057, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25439, 25434, 25442, 25444, 4, 24084, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25541, 25536, 25600, 25802, 5, 24515, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25605, 25601, 25614, 25836, 5, 24549, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25642, 25636, 25801, 25851, 5, 25069, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25807, 25802, 25835, 25920, 5, 25138, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25857, 25851, 25919, 25940, 5, 25278, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25925, 25920, 25934, 25975, 5, 25313, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25946, 25940, 25974, 25981, 5, 25383, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 25986, 25981, 25995, 25997, 5, 25428, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 26007, 26003, 26057, 26068, 6, 25802, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 26072, 26068, 26097, 26102, 6, 25920, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 26106, 26102, 26115, 26117, 6, 25975, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4067, 4064, 4980, 5869, 5968, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4987, 4981, 5474, 5970, 6072, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5484, 5475, 5704, 6074, 6142, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5713, 5705, 5739, 6144, 6157, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5781, 5778, 5868, 6178, 6228, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5875, 5869, 5969, 6230, 6302, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5979, 5970, 6073, 6304, 6357, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6083, 6074, 6143, 6359, 6392, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6181, 6178, 6229, 6407, 6431, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6239, 6230, 6303, 6433, 6466, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6313, 6304, 6358, 6468, 6491, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6367, 6359, 6393, 6493, 6506, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6410, 6407, 6432, 6511, 6520, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6442, 6433, 6467, 6522, 6533, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6476, 6468, 6492, 6535, 6542, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6514, 6511, 6521, 6548, 6550, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6530, 6522, 6534, 6552, 6553, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7714, 7711, 9675, 11663, 13171, 1, 4064, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9682, 9676, 10912, 13173, 14132, 1, 4981, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 10922, 10913, 11542, 14134, 14236, 1, 5475, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11552, 11543, 11612, 14238, 14271, 1, 5705, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11669, 11663, 13172, 14298, 14416, 1, 5869, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13182, 13173, 14133, 14418, 14505, 1, 5970, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14143, 14134, 14237, 14507, 14560, 1, 6074, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14246, 14238, 14272, 14562, 14575, 1, 6144, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14307, 14298, 14417, 14588, 14650, 1, 6230, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14427, 14418, 14506, 14652, 14695, 1, 6304, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14516, 14507, 14561, 14697, 14720, 1, 6359, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14597, 14588, 14651, 14731, 14758, 1, 6433, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14661, 14652, 14696, 14760, 14777, 1, 6468, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14705, 14697, 14721, 14779, 14786, 1, 6493, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14740, 14731, 14759, 14790, 14798, 1, 6522, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 14768, 14760, 14778, 14800, 14804, 1, 6535, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 15432, 15429, 16447, 17927, 19944, 2, 9676, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 16457, 16448, 17633, 19946, 20451, 2, 10913, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 17643, 17634, 17851, 20453, 20521, 2, 11543, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 17860, 17852, 17886, 20523, 20536, 2, 11613, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 17936, 17927, 19945, 20559, 20677, 2, 13173, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 19955, 19946, 20452, 20679, 20766, 2, 14134, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20462, 20453, 20522, 20768, 20801, 2, 14238, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20568, 20559, 20678, 20818, 20880, 2, 14418, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20688, 20679, 20767, 20882, 20925, 2, 14507, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20776, 20768, 20802, 20927, 20940, 2, 14562, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20827, 20818, 20881, 20947, 20974, 2, 14652, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20891, 20882, 20926, 20976, 20993, 2, 14697, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20956, 20947, 20975, 21001, 21009, 2, 14760, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 20984, 20976, 20994, 21011, 21015, 2, 14779, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 21306, 21303, 21754, 22355, 23425, 3, 16448, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 21764, 21755, 22234, 23427, 23643, 3, 17634, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 22243, 22235, 22303, 23645, 23678, 3, 17852, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 22364, 22355, 23426, 23705, 23823, 3, 19946, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 23436, 23427, 23644, 23825, 23893, 3, 20453, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 23653, 23645, 23679, 23895, 23908, 3, 20523, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 23714, 23705, 23824, 23921, 23983, 3, 20679, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 23834, 23825, 23894, 23985, 24018, 3, 20768, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 23930, 23921, 23984, 24029, 24056, 3, 20882, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 23993, 23985, 24019, 24058, 24071, 3, 20927, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 24038, 24029, 24057, 24075, 24083, 3, 20976, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 24195, 24192, 24349, 24590, 25068, 4, 21755, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 24358, 24350, 24515, 25070, 25137, 4, 22235, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 24523, 24516, 24549, 25139, 25152, 4, 22304, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 24599, 24590, 25069, 25175, 25277, 4, 23427, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25078, 25070, 25138, 25279, 25312, 4, 23645, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25184, 25175, 25278, 25329, 25382, 4, 23825, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25287, 25279, 25313, 25384, 25397, 4, 23895, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25338, 25329, 25383, 25404, 25427, 4, 23985, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25412, 25404, 25428, 25435, 25442, 4, 24058, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25482, 25480, 25536, 25637, 25801, 5, 24350, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25544, 25537, 25601, 25803, 25835, 5, 24516, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25645, 25637, 25802, 25852, 25919, 5, 25070, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25810, 25803, 25836, 25921, 25934, 5, 25139, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25860, 25852, 25920, 25941, 25974, 5, 25279, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 25949, 25941, 25975, 25982, 25995, 5, 25384, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2612, 2610, 4064, 5779, 5868, 13171, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4070, 4065, 4981, 5870, 5969, 14132, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4990, 4982, 5475, 5971, 6073, 14236, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5487, 5476, 5705, 6075, 6143, 14271, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5784, 5779, 5869, 6179, 6229, 14416, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5878, 5870, 5970, 6231, 6303, 14505, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5982, 5971, 6074, 6305, 6358, 14560, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6086, 6075, 6144, 6360, 6393, 14575, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6184, 6179, 6230, 6408, 6432, 14650, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6242, 6231, 6304, 6434, 6467, 14695, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6316, 6305, 6359, 6469, 6492, 14720, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6413, 6408, 6433, 6512, 6521, 14758, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6445, 6434, 6468, 6523, 6534, 14777, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6517, 6512, 6522, 6549, 6551, 14798, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7717, 7712, 9676, 11664, 13172, 19944, 1, 4065, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9685, 9677, 10913, 13174, 14133, 20451, 1, 4982, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 10925, 10914, 11543, 14135, 14237, 20521, 1, 5476, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11555, 11544, 11613, 14239, 14272, 20536, 1, 5706, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11672, 11664, 13173, 14299, 14417, 20677, 1, 5870, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13185, 13174, 14134, 14419, 14506, 20766, 1, 5971, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 14146, 14135, 14238, 14508, 14561, 20801, 1, 6075, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 14310, 14299, 14418, 14589, 14651, 20880, 1, 6231, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 14430, 14419, 14507, 14653, 14696, 20925, 1, 6305, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 14519, 14508, 14562, 14698, 14721, 20940, 1, 6360, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 14600, 14589, 14652, 14732, 14759, 20974, 1, 6434, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 14664, 14653, 14697, 14761, 14778, 20993, 1, 6469, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 14743, 14732, 14760, 14791, 14799, 21009, 1, 6523, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 15435, 15430, 16448, 17928, 19945, 23425, 2, 9677, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 16460, 16449, 17634, 19947, 20452, 23643, 2, 10914, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 17646, 17635, 17852, 20454, 20522, 23678, 2, 11544, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 17939, 17928, 19946, 20560, 20678, 23823, 2, 13174, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 19958, 19947, 20453, 20680, 20767, 23893, 2, 14135, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 20465, 20454, 20523, 20769, 20802, 23908, 2, 14239, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 20571, 20560, 20679, 20819, 20881, 23983, 2, 14419, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 20691, 20680, 20768, 20883, 20926, 24018, 2, 14508, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 20830, 20819, 20882, 20948, 20975, 24056, 2, 14653, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 20894, 20883, 20927, 20977, 20994, 24071, 2, 14698, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 20959, 20948, 20976, 21002, 21010, 24083, 2, 14761, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 21309, 21304, 21755, 22356, 23426, 25068, 3, 16449, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 21767, 21756, 22235, 23428, 23644, 25137, 3, 17635, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 22246, 22236, 22304, 23646, 23679, 25152, 3, 17853, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 22367, 22356, 23427, 23706, 23824, 25277, 3, 19947, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 23439, 23428, 23645, 23826, 23894, 25312, 3, 20454, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 23717, 23706, 23825, 23922, 23984, 25382, 3, 20680, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 23837, 23826, 23895, 23986, 24019, 25397, 3, 20769, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 23933, 23922, 23985, 24030, 24057, 25427, 3, 20883, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 24041, 24030, 24058, 24076, 24084, 25442, 3, 20977, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 24198, 24193, 24350, 24591, 25069, 25801, 4, 21756, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 24361, 24351, 24516, 25071, 25138, 25835, 4, 22236, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 24602, 24591, 25070, 25176, 25278, 25919, 4, 23428, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 25081, 25071, 25139, 25280, 25313, 25934, 4, 23646, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 25187, 25176, 25279, 25330, 25383, 25974, 4, 23826, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 25341, 25330, 25384, 25405, 25428, 25995, 4, 23986, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 25485, 25481, 25537, 25638, 25802, 26057, 5, 24351, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 25648, 25638, 25803, 25853, 25920, 26097, 5, 25071, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 25863, 25853, 25921, 25942, 25975, 26115, 5, 25280, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1, 0, 2610, 5743, 5778, 11662, 14296, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2615, 2611, 4065, 5780, 5869, 13172, 14416, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 4073, 4066, 4982, 5871, 5970, 14133, 14505, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 4993, 4983, 5476, 5972, 6074, 14237, 14560, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5490, 5477, 5706, 6076, 6144, 14272, 14575, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5744, 5743, 5779, 6158, 6178, 14297, 14586, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5787, 5780, 5870, 6180, 6230, 14417, 14650, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5881, 5871, 5971, 6232, 6304, 14506, 14695, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5985, 5972, 6075, 6306, 6359, 14561, 14720, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6159, 6158, 6179, 6397, 6407, 14587, 14729, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6187, 6180, 6231, 6409, 6433, 14651, 14758, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6245, 6232, 6305, 6435, 6468, 14696, 14777, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6398, 6397, 6408, 6507, 6511, 14730, 14788, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6416, 6409, 6434, 6513, 6522, 14759, 14798, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6508, 6507, 6512, 6547, 6548, 14789, 14806, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6557, 6556, 7712, 11627, 11663, 17926, 20557, 1, 2611, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7720, 7713, 9677, 11665, 13173, 19945, 20677, 1, 4066, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 9688, 9678, 10914, 13175, 14134, 20452, 20766, 1, 4983, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 10928, 10915, 11544, 14136, 14238, 20522, 20801, 1, 5477, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 11628, 11627, 11664, 14276, 14298, 20558, 20816, 1, 5780, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 11675, 11665, 13174, 14300, 14418, 20678, 20880, 1, 5871, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 13188, 13175, 14135, 14420, 14507, 20767, 20925, 1, 5972, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14149, 14136, 14239, 14509, 14562, 20802, 20940, 1, 6076, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14277, 14276, 14299, 14576, 14588, 20817, 20945, 1, 6180, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14313, 14300, 14419, 14590, 14652, 20881, 20974, 1, 6232, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14433, 14420, 14508, 14654, 14697, 20926, 20993, 1, 6306, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14577, 14576, 14589, 14725, 14731, 20946, 20999, 1, 6409, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14603, 14590, 14653, 14733, 14760, 20975, 21009, 1, 6435, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14726, 14725, 14732, 14787, 14790, 21000, 21016, 1, 6513, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 14810, 14809, 15430, 17890, 17927, 22354, 23703, 2, 7713, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 15438, 15431, 16449, 17929, 19946, 23426, 23823, 2, 9678, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 16463, 16450, 17635, 19948, 20453, 23644, 23893, 2, 10915, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 17649, 17636, 17853, 20455, 20523, 23679, 23908, 2, 11545, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 17891, 17890, 17928, 20537, 20559, 23704, 23919, 2, 11665, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 17942, 17929, 19947, 20561, 20679, 23824, 23983, 2, 13175, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 19961, 19948, 20454, 20681, 20768, 23894, 24018, 2, 14136, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 20538, 20537, 20560, 20806, 20818, 23920, 24027, 2, 14300, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 20574, 20561, 20680, 20820, 20882, 23984, 24056, 2, 14420, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 20694, 20681, 20769, 20884, 20927, 24019, 24071, 2, 14509, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 20807, 20806, 20819, 20941, 20947, 24028, 24073, 2, 14590, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 20833, 20820, 20883, 20949, 20976, 24057, 24083, 2, 14654, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 20942, 20941, 20948, 20998, 21001, 24074, 24088, 2, 14733, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 21020, 21019, 21304, 22318, 22355, 24589, 25173, 3, 15431, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 21312, 21305, 21756, 22357, 23427, 25069, 25277, 3, 16450, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 21770, 21757, 22236, 23429, 23645, 25138, 25312, 3, 17636, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 22319, 22318, 22356, 23683, 23705, 25174, 25327, 3, 17929, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 22370, 22357, 23428, 23707, 23825, 25278, 25382, 3, 19948, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 23442, 23429, 23646, 23827, 23895, 25313, 25397, 3, 20455, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 23684, 23683, 23706, 23909, 23921, 25328, 25402, 3, 20561, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 23720, 23707, 23826, 23923, 23985, 25383, 25427, 3, 20681, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 23910, 23909, 23922, 24023, 24029, 25403, 25433, 3, 20820, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 23936, 23923, 23986, 24031, 24058, 25428, 25442, 3, 20884, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 24024, 24023, 24030, 24072, 24075, 25434, 25443, 3, 20949, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 24092, 24091, 24193, 24553, 24590, 25636, 25850, 4, 21305, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 24201, 24194, 24351, 24592, 25070, 25802, 25919, 4, 21757, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 24554, 24553, 24591, 25153, 25175, 25851, 25939, 4, 22357, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 24605, 24592, 25071, 25177, 25279, 25920, 25974, 4, 23429, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25154, 25153, 25176, 25317, 25329, 25940, 25980, 4, 23707, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25190, 25177, 25280, 25331, 25384, 25975, 25995, 4, 23827, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25318, 25317, 25330, 25398, 25404, 25981, 25996, 4, 23923, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25399, 25398, 25405, 25432, 25435, 25997, 26001, 4, 24031, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25446, 25445, 25481, 25615, 25637, 26003, 26067, 5, 24194, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25616, 25615, 25638, 25840, 25852, 26068, 26101, 5, 24592, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25841, 25840, 25853, 25935, 25941, 26102, 26116, 5, 25177, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 25936, 25935, 25942, 25979, 25982, 26117, 26121, 5, 25331, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5493, 5473, 5478, 6072, 6077, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5716, 5703, 5707, 6142, 6145, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5988, 5968, 5973, 6302, 6307, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6089, 6072, 6077, 6357, 6361, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6151, 6142, 6145, 6392, 6394, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6248, 6228, 6233, 6431, 6436, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6319, 6302, 6307, 6466, 6470, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6370, 6357, 6361, 6491, 6494, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6448, 6431, 6436, 6520, 6524, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6479, 6466, 6470, 6533, 6536, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6500, 6491, 6494, 6542, 6544, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 10931, 10911, 10916, 14132, 14137, 1, 5473, 5478, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11558, 11541, 11546, 14236, 14240, 1, 5703, 5707, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11620, 11611, 11614, 14271, 14273, 1, 5738, 5740, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13191, 13171, 13176, 14416, 14421, 1, 5968, 5973, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14152, 14132, 14137, 14505, 14510, 1, 6072, 6077, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14249, 14236, 14240, 14560, 14563, 1, 6142, 6145, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14316, 14296, 14301, 14586, 14591, 1, 6228, 6233, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14436, 14416, 14421, 14650, 14655, 1, 6302, 6307, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14522, 14505, 14510, 14695, 14699, 1, 6357, 6361, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14569, 14560, 14563, 14720, 14722, 1, 6392, 6394, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14606, 14586, 14591, 14729, 14734, 1, 6431, 6436, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14667, 14650, 14655, 14758, 14762, 1, 6466, 6470, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14708, 14695, 14699, 14777, 14780, 1, 6491, 6494, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14746, 14729, 14734, 14788, 14792, 1, 6520, 6524, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 14771, 14758, 14762, 14798, 14801, 1, 6533, 6536, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 16466, 16446, 16451, 19944, 19949, 2, 10911, 10916, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 17652, 17632, 17637, 20451, 20456, 2, 11541, 11546, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 17863, 17850, 17854, 20521, 20524, 2, 11611, 11614, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 17945, 17925, 17930, 20557, 20562, 2, 13171, 13176, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 19964, 19944, 19949, 20677, 20682, 2, 14132, 14137, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20468, 20451, 20456, 20766, 20770, 2, 14236, 14240, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20530, 20521, 20524, 20801, 20803, 2, 14271, 14273, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20577, 20557, 20562, 20816, 20821, 2, 14416, 14421, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20697, 20677, 20682, 20880, 20885, 2, 14505, 14510, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20779, 20766, 20770, 20925, 20928, 2, 14560, 14563, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20836, 20816, 20821, 20945, 20950, 2, 14650, 14655, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20897, 20880, 20885, 20974, 20978, 2, 14695, 14699, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20934, 20925, 20928, 20993, 20995, 2, 14720, 14722, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20962, 20945, 20950, 20999, 21003, 2, 14758, 14762, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 20987, 20974, 20978, 21009, 21012, 2, 14777, 14780, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 21773, 21753, 21758, 23425, 23430, 3, 17632, 17637, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 22249, 22233, 22237, 23643, 23647, 3, 17850, 17854, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 22311, 22302, 22305, 23678, 23680, 3, 17885, 17887, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 22373, 22353, 22358, 23703, 23708, 3, 19944, 19949, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 23445, 23425, 23430, 23823, 23828, 3, 20451, 20456, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 23656, 23643, 23647, 23893, 23896, 3, 20521, 20524, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 23723, 23703, 23708, 23919, 23924, 3, 20677, 20682, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 23840, 23823, 23828, 23983, 23987, 3, 20766, 20770, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 23902, 23893, 23896, 24018, 24020, 3, 20801, 20803, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 23939, 23919, 23924, 24027, 24032, 3, 20880, 20885, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 23996, 23983, 23987, 24056, 24059, 3, 20925, 20928, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 24044, 24027, 24032, 24073, 24077, 3, 20974, 20978, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 24065, 24056, 24059, 24083, 24085, 3, 20993, 20995, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 24364, 24348, 24352, 25068, 25072, 4, 22233, 22237, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 24526, 24514, 24517, 25137, 25140, 4, 22302, 22305, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 24608, 24588, 24593, 25173, 25178, 4, 23425, 23430, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25084, 25068, 25072, 25277, 25281, 4, 23643, 23647, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25146, 25137, 25140, 25312, 25314, 4, 23678, 23680, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25193, 25173, 25178, 25327, 25332, 4, 23823, 23828, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25290, 25277, 25281, 25382, 25385, 4, 23893, 23896, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25344, 25327, 25332, 25402, 25406, 4, 23983, 23987, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25391, 25382, 25385, 25427, 25429, 4, 24018, 24020, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25415, 25402, 25406, 25433, 25436, 4, 24056, 24059, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25547, 25535, 25538, 25801, 25804, 5, 24514, 24517, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25608, 25600, 25602, 25835, 25837, 5, 24548, 24550, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25651, 25635, 25639, 25850, 25854, 5, 25068, 25072, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25813, 25801, 25804, 25919, 25922, 5, 25137, 25140, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25866, 25850, 25854, 25939, 25943, 5, 25277, 25281, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25928, 25919, 25922, 25974, 25976, 5, 25312, 25314, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25952, 25939, 25943, 25980, 25983, 5, 25382, 25385, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 25989, 25980, 25983, 25996, 25998, 5, 25427, 25429, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 26010, 26002, 26004, 26067, 26069, 6, 25801, 25804, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 26061, 26057, 26058, 26097, 26098, 6, 25835, 25837, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 26075, 26067, 26069, 26101, 26103, 6, 25919, 25922, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 26109, 26101, 26103, 26116, 26118, 6, 25974, 25976, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4996, 4980, 4984, 5478, 5969, 5976, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5499, 5474, 5481, 5707, 6073, 6080, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5722, 5704, 5710, 5740, 6143, 6148, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5884, 5868, 5872, 5973, 6229, 6236, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5994, 5969, 5976, 6077, 6303, 6310, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6095, 6073, 6080, 6145, 6358, 6364, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6254, 6229, 6236, 6307, 6432, 6439, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6325, 6303, 6310, 6361, 6467, 6473, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6376, 6358, 6364, 6394, 6492, 6497, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6454, 6432, 6439, 6470, 6521, 6527, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6485, 6467, 6473, 6494, 6534, 6539, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9691, 9675, 9679, 10916, 13172, 13179, 1, 4980, 4984, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 10937, 10912, 10919, 11546, 14133, 14140, 1, 5474, 5481, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11564, 11542, 11549, 11614, 14237, 14243, 1, 5704, 5710, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11678, 11662, 11666, 13176, 14297, 14304, 1, 5868, 5872, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13197, 13172, 13179, 14137, 14417, 14424, 1, 5969, 5976, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14158, 14133, 14140, 14240, 14506, 14513, 1, 6073, 6080, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14255, 14237, 14243, 14273, 14561, 14566, 1, 6143, 6148, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14322, 14297, 14304, 14421, 14587, 14594, 1, 6229, 6236, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14442, 14417, 14424, 14510, 14651, 14658, 1, 6303, 6310, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14528, 14506, 14513, 14563, 14696, 14702, 1, 6358, 6364, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14612, 14587, 14594, 14655, 14730, 14737, 1, 6432, 6439, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14673, 14651, 14658, 14699, 14759, 14765, 1, 6467, 6473, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14714, 14696, 14702, 14722, 14778, 14783, 1, 6492, 6497, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 14752, 14730, 14737, 14762, 14789, 14795, 1, 6521, 6527, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 16472, 16447, 16454, 17637, 19945, 19952, 2, 10912, 10919, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 17658, 17633, 17640, 17854, 20452, 20459, 2, 11542, 11549, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 17869, 17851, 17857, 17887, 20522, 20527, 2, 11612, 11617, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 17951, 17926, 17933, 19949, 20558, 20565, 2, 13172, 13179, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 19970, 19945, 19952, 20456, 20678, 20685, 2, 14133, 14140, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 20474, 20452, 20459, 20524, 20767, 20773, 2, 14237, 14243, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 20583, 20558, 20565, 20682, 20817, 20824, 2, 14417, 14424, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 20703, 20678, 20685, 20770, 20881, 20888, 2, 14506, 14513, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 20785, 20767, 20773, 20803, 20926, 20931, 2, 14561, 14566, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 20842, 20817, 20824, 20885, 20946, 20953, 2, 14651, 14658, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 20903, 20881, 20888, 20928, 20975, 20981, 2, 14696, 14702, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 20968, 20946, 20953, 20978, 21000, 21006, 2, 14759, 14765, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 21779, 21754, 21761, 22237, 23426, 23433, 3, 17633, 17640, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 22255, 22234, 22240, 22305, 23644, 23650, 3, 17851, 17857, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 22379, 22354, 22361, 23430, 23704, 23711, 3, 19945, 19952, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 23451, 23426, 23433, 23647, 23824, 23831, 3, 20452, 20459, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 23662, 23644, 23650, 23680, 23894, 23899, 3, 20522, 20527, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 23729, 23704, 23711, 23828, 23920, 23927, 3, 20678, 20685, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 23846, 23824, 23831, 23896, 23984, 23990, 3, 20767, 20773, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 23945, 23920, 23927, 23987, 24028, 24035, 3, 20881, 20888, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 24002, 23984, 23990, 24020, 24057, 24062, 3, 20926, 20931, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 24050, 24028, 24035, 24059, 24074, 24080, 3, 20975, 20981, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 24370, 24349, 24355, 24517, 25069, 25075, 4, 22234, 22240, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 24532, 24515, 24520, 24550, 25138, 25143, 4, 22303, 22308, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 24614, 24589, 24596, 25072, 25174, 25181, 4, 23426, 23433, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25090, 25069, 25075, 25140, 25278, 25284, 4, 23644, 23650, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25199, 25174, 25181, 25281, 25328, 25335, 4, 23824, 23831, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25296, 25278, 25284, 25314, 25383, 25388, 4, 23894, 23899, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25350, 25328, 25335, 25385, 25403, 25409, 4, 23984, 23990, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25421, 25403, 25409, 25429, 25434, 25439, 4, 24057, 24062, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25553, 25536, 25541, 25602, 25802, 25807, 5, 24515, 24520, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25657, 25636, 25642, 25804, 25851, 25857, 5, 25069, 25075, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25819, 25802, 25807, 25837, 25920, 25925, 5, 25138, 25143, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25872, 25851, 25857, 25922, 25940, 25946, 5, 25278, 25284, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 25958, 25940, 25946, 25976, 25981, 25986, 5, 25383, 25388, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 26016, 26003, 26007, 26058, 26068, 26072, 6, 25802, 25807, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 26081, 26068, 26072, 26098, 26102, 26106, 6, 25920, 25925, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4076, 4064, 4067, 4984, 5869, 5875, 5968, 5973, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5002, 4981, 4987, 5481, 5970, 5979, 6072, 6077, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5505, 5475, 5484, 5710, 6074, 6083, 6142, 6145, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5790, 5778, 5781, 5872, 6178, 6181, 6228, 6233, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5890, 5869, 5875, 5976, 6230, 6239, 6302, 6307, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6000, 5970, 5979, 6080, 6304, 6313, 6357, 6361, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6101, 6074, 6083, 6148, 6359, 6367, 6392, 6394, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6190, 6178, 6181, 6236, 6407, 6410, 6431, 6436, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6260, 6230, 6239, 6310, 6433, 6442, 6466, 6470, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6331, 6304, 6313, 6364, 6468, 6476, 6491, 6494, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6419, 6407, 6410, 6439, 6511, 6514, 6520, 6524, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6460, 6433, 6442, 6473, 6522, 6530, 6533, 6536, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7723, 7711, 7714, 9679, 11663, 11669, 13171, 13176, 1, 4064, 4067, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 9697, 9676, 9682, 10919, 13173, 13182, 14132, 14137, 1, 4981, 4987, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 10943, 10913, 10922, 11549, 14134, 14143, 14236, 14240, 1, 5475, 5484, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 11570, 11543, 11552, 11617, 14238, 14246, 14271, 14273, 1, 5705, 5713, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 11684, 11663, 11669, 13179, 14298, 14307, 14416, 14421, 1, 5869, 5875, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 13203, 13173, 13182, 14140, 14418, 14427, 14505, 14510, 1, 5970, 5979, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 14164, 14134, 14143, 14243, 14507, 14516, 14560, 14563, 1, 6074, 6083, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 14328, 14298, 14307, 14424, 14588, 14597, 14650, 14655, 1, 6230, 6239, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 14448, 14418, 14427, 14513, 14652, 14661, 14695, 14699, 1, 6304, 6313, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 14534, 14507, 14516, 14566, 14697, 14705, 14720, 14722, 1, 6359, 6367, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 14618, 14588, 14597, 14658, 14731, 14740, 14758, 14762, 1, 6433, 6442, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 14679, 14652, 14661, 14702, 14760, 14768, 14777, 14780, 1, 6468, 6476, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 15441, 15429, 15432, 16454, 17927, 17936, 19944, 19949, 2, 9676, 9682, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 16478, 16448, 16457, 17640, 19946, 19955, 20451, 20456, 2, 10913, 10922, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 17664, 17634, 17643, 17857, 20453, 20462, 20521, 20524, 2, 11543, 11552, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 17957, 17927, 17936, 19952, 20559, 20568, 20677, 20682, 2, 13173, 13182, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 19976, 19946, 19955, 20459, 20679, 20688, 20766, 20770, 2, 14134, 14143, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 20480, 20453, 20462, 20527, 20768, 20776, 20801, 20803, 2, 14238, 14246, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 20589, 20559, 20568, 20685, 20818, 20827, 20880, 20885, 2, 14418, 14427, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 20709, 20679, 20688, 20773, 20882, 20891, 20925, 20928, 2, 14507, 14516, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 20848, 20818, 20827, 20888, 20947, 20956, 20974, 20978, 2, 14652, 14661, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 20909, 20882, 20891, 20931, 20976, 20984, 20993, 20995, 2, 14697, 14705, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 21315, 21303, 21306, 21761, 22355, 22364, 23425, 23430, 3, 16448, 16457, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 21785, 21755, 21764, 22240, 23427, 23436, 23643, 23647, 3, 17634, 17643, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 22261, 22235, 22243, 22308, 23645, 23653, 23678, 23680, 3, 17852, 17860, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 22385, 22355, 22364, 23433, 23705, 23714, 23823, 23828, 3, 19946, 19955, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 23457, 23427, 23436, 23650, 23825, 23834, 23893, 23896, 3, 20453, 20462, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 23735, 23705, 23714, 23831, 23921, 23930, 23983, 23987, 3, 20679, 20688, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 23852, 23825, 23834, 23899, 23985, 23993, 24018, 24020, 3, 20768, 20776, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 23951, 23921, 23930, 23990, 24029, 24038, 24056, 24059, 3, 20882, 20891, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 24204, 24192, 24195, 24355, 24590, 24599, 25068, 25072, 4, 21755, 21764, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 24376, 24350, 24358, 24520, 25070, 25078, 25137, 25140, 4, 22235, 22243, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 24620, 24590, 24599, 25075, 25175, 25184, 25277, 25281, 4, 23427, 23436, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 25096, 25070, 25078, 25143, 25279, 25287, 25312, 25314, 4, 23645, 23653, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 25205, 25175, 25184, 25284, 25329, 25338, 25382, 25385, 4, 23825, 23834, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 25356, 25329, 25338, 25388, 25404, 25412, 25427, 25429, 4, 23985, 23993, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 25488, 25480, 25482, 25541, 25637, 25645, 25801, 25804, 5, 24350, 24358, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 25559, 25537, 25544, 25605, 25803, 25810, 25835, 25837, 5, 24516, 24523, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 25663, 25637, 25645, 25807, 25852, 25860, 25919, 25922, 5, 25070, 25078, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 25878, 25852, 25860, 25925, 25941, 25949, 25974, 25976, 5, 25279, 25287, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2618, 2610, 2612, 4067, 5779, 5784, 5868, 5872, 13176, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4082, 4065, 4070, 4987, 5870, 5878, 5969, 5976, 14137, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5008, 4982, 4990, 5484, 5971, 5982, 6073, 6080, 14240, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5511, 5476, 5487, 5713, 6075, 6086, 6143, 6148, 14273, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5796, 5779, 5784, 5875, 6179, 6184, 6229, 6236, 14421, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5896, 5870, 5878, 5979, 6231, 6242, 6303, 6310, 14510, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 6006, 5971, 5982, 6083, 6305, 6316, 6358, 6364, 14563, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 6196, 6179, 6184, 6239, 6408, 6413, 6432, 6439, 14655, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 6266, 6231, 6242, 6313, 6434, 6445, 6467, 6473, 14699, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 6425, 6408, 6413, 6442, 6512, 6517, 6521, 6527, 14762, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 7729, 7712, 7717, 9682, 11664, 11672, 13172, 13179, 19949, 1, 4065, 4070, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 9703, 9677, 9685, 10922, 13174, 13185, 14133, 14140, 20456, 1, 4982, 4990, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 10949, 10914, 10925, 11552, 14135, 14146, 14237, 14243, 20524, 1, 5476, 5487, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 11690, 11664, 11672, 13182, 14299, 14310, 14417, 14424, 20682, 1, 5870, 5878, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13209, 13174, 13185, 14143, 14419, 14430, 14506, 14513, 20770, 1, 5971, 5982, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 14170, 14135, 14146, 14246, 14508, 14519, 14561, 14566, 20803, 1, 6075, 6086, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 14334, 14299, 14310, 14427, 14589, 14600, 14651, 14658, 20885, 1, 6231, 6242, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 14454, 14419, 14430, 14516, 14653, 14664, 14696, 14702, 20928, 1, 6305, 6316, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 14624, 14589, 14600, 14661, 14732, 14743, 14759, 14765, 20978, 1, 6434, 6445, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 15447, 15430, 15435, 16457, 17928, 17939, 19945, 19952, 23430, 2, 9677, 9685, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 16484, 16449, 16460, 17643, 19947, 19958, 20452, 20459, 23647, 2, 10914, 10925, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 17670, 17635, 17646, 17860, 20454, 20465, 20522, 20527, 23680, 2, 11544, 11555, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 17963, 17928, 17939, 19955, 20560, 20571, 20678, 20685, 23828, 2, 13174, 13185, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 19982, 19947, 19958, 20462, 20680, 20691, 20767, 20773, 23896, 2, 14135, 14146, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 20595, 20560, 20571, 20688, 20819, 20830, 20881, 20888, 23987, 2, 14419, 14430, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 20715, 20680, 20691, 20776, 20883, 20894, 20926, 20931, 24020, 2, 14508, 14519, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 20854, 20819, 20830, 20891, 20948, 20959, 20975, 20981, 24059, 2, 14653, 14664, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 21321, 21304, 21309, 21764, 22356, 22367, 23426, 23433, 25072, 3, 16449, 16460, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 21791, 21756, 21767, 22243, 23428, 23439, 23644, 23650, 25140, 3, 17635, 17646, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 22391, 22356, 22367, 23436, 23706, 23717, 23824, 23831, 25281, 3, 19947, 19958, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 23463, 23428, 23439, 23653, 23826, 23837, 23894, 23899, 25314, 3, 20454, 20465, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 23741, 23706, 23717, 23834, 23922, 23933, 23984, 23990, 25385, 3, 20680, 20691, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 23957, 23922, 23933, 23993, 24030, 24041, 24057, 24062, 25429, 3, 20883, 20894, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 24210, 24193, 24198, 24358, 24591, 24602, 25069, 25075, 25804, 4, 21756, 21767, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 24382, 24351, 24361, 24523, 25071, 25081, 25138, 25143, 25837, 4, 22236, 22246, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 24626, 24591, 24602, 25078, 25176, 25187, 25278, 25284, 25922, 4, 23428, 23439, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 25211, 25176, 25187, 25287, 25330, 25341, 25383, 25388, 25976, 4, 23826, 23837, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 25494, 25481, 25485, 25544, 25638, 25648, 25802, 25807, 26058, 5, 24351, 24361, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 25669, 25638, 25648, 25810, 25853, 25863, 25920, 25925, 26098, 5, 25071, 25081, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 4, 0, 1, 2612, 5743, 5744, 5778, 5781, 11666, 14296, 14301, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 2624, 2611, 2615, 4070, 5780, 5787, 5869, 5875, 13179, 14416, 14421, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 4088, 4066, 4073, 4990, 5871, 5881, 5970, 5979, 14140, 14505, 14510, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5014, 4983, 4993, 5487, 5972, 5985, 6074, 6083, 14243, 14560, 14563, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5747, 5743, 5744, 5784, 6158, 6159, 6178, 6181, 14304, 14586, 14591, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5802, 5780, 5787, 5878, 6180, 6187, 6230, 6239, 14424, 14650, 14655, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5902, 5871, 5881, 5982, 6232, 6245, 6304, 6313, 14513, 14695, 14699, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 6162, 6158, 6159, 6184, 6397, 6398, 6407, 6410, 14594, 14729, 14734, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 6202, 6180, 6187, 6242, 6409, 6416, 6433, 6442, 14658, 14758, 14762, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 6401, 6397, 6398, 6413, 6507, 6508, 6511, 6514, 14737, 14788, 14792, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 6560, 6556, 6557, 7717, 11627, 11628, 11663, 11669, 17933, 20557, 20562, 1, 2611, 2615, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 7735, 7713, 7720, 9685, 11665, 11675, 13173, 13182, 19952, 20677, 20682, 1, 4066, 4073, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 9709, 9678, 9688, 10925, 13175, 13188, 14134, 14143, 20459, 20766, 20770, 1, 4983, 4993, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 10955, 10915, 10928, 11555, 14136, 14149, 14238, 14246, 20527, 20801, 20803, 1, 5477, 5490, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 11631, 11627, 11628, 11672, 14276, 14277, 14298, 14307, 20565, 20816, 20821, 1, 5780, 5787, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 11696, 11665, 11675, 13185, 14300, 14313, 14418, 14427, 20685, 20880, 20885, 1, 5871, 5881, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 13215, 13175, 13188, 14146, 14420, 14433, 14507, 14516, 20773, 20925, 20928, 1, 5972, 5985, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 14280, 14276, 14277, 14310, 14576, 14577, 14588, 14597, 20824, 20945, 20950, 1, 6180, 6187, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 14340, 14300, 14313, 14430, 14590, 14603, 14652, 14661, 20888, 20974, 20978, 1, 6232, 6245, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 14580, 14576, 14577, 14600, 14725, 14726, 14731, 14740, 20953, 20999, 21003, 1, 6409, 6416, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 14813, 14809, 14810, 15435, 17890, 17891, 17927, 17936, 22361, 23703, 23708, 2, 7713, 7720, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 15453, 15431, 15438, 16460, 17929, 17942, 19946, 19955, 23433, 23823, 23828, 2, 9678, 9688, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 16490, 16450, 16463, 17646, 19948, 19961, 20453, 20462, 23650, 23893, 23896, 2, 10915, 10928, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 17894, 17890, 17891, 17939, 20537, 20538, 20559, 20568, 23711, 23919, 23924, 2, 11665, 11675, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 17969, 17929, 17942, 19958, 20561, 20574, 20679, 20688, 23831, 23983, 23987, 2, 13175, 13188, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 19988, 19948, 19961, 20465, 20681, 20694, 20768, 20776, 23899, 24018, 24020, 2, 14136, 14149, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 20541, 20537, 20538, 20571, 20806, 20807, 20818, 20827, 23927, 24027, 24032, 2, 14300, 14313, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 20601, 20561, 20574, 20691, 20820, 20833, 20882, 20891, 23990, 24056, 24059, 2, 14420, 14433, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 20810, 20806, 20807, 20830, 20941, 20942, 20947, 20956, 24035, 24073, 24077, 2, 14590, 14603, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 21023, 21019, 21020, 21309, 22318, 22319, 22355, 22364, 24596, 25173, 25178, 3, 15431, 15438, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 21327, 21305, 21312, 21767, 22357, 22370, 23427, 23436, 25075, 25277, 25281, 3, 16450, 16463, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 21797, 21757, 21770, 22246, 23429, 23442, 23645, 23653, 25143, 25312, 25314, 3, 17636, 17649, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 22322, 22318, 22319, 22367, 23683, 23684, 23705, 23714, 25181, 25327, 25332, 3, 17929, 17942, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 22397, 22357, 22370, 23439, 23707, 23720, 23825, 23834, 25284, 25382, 25385, 3, 19948, 19961, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 23687, 23683, 23684, 23717, 23909, 23910, 23921, 23930, 25335, 25402, 25406, 3, 20561, 20574, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 23747, 23707, 23720, 23837, 23923, 23936, 23985, 23993, 25388, 25427, 25429, 3, 20681, 20694, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 23913, 23909, 23910, 23933, 24023, 24024, 24029, 24038, 25409, 25433, 25436, 3, 20820, 20833, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 24095, 24091, 24092, 24198, 24553, 24554, 24590, 24599, 25642, 25850, 25854, 4, 21305, 21312, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 24216, 24194, 24201, 24361, 24592, 24605, 25070, 25078, 25807, 25919, 25922, 4, 21757, 21770, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 24557, 24553, 24554, 24602, 25153, 25154, 25175, 25184, 25857, 25939, 25943, 4, 22357, 22370, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 24632, 24592, 24605, 25081, 25177, 25190, 25279, 25287, 25925, 25974, 25976, 4, 23429, 23442, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 25157, 25153, 25154, 25187, 25317, 25318, 25329, 25338, 25946, 25980, 25983, 4, 23707, 23720, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 25321, 25317, 25318, 25341, 25398, 25399, 25404, 25412, 25986, 25996, 25998, 4, 23923, 23936, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 25449, 25445, 25446, 25485, 25615, 25616, 25637, 25645, 26007, 26067, 26069, 5, 24194, 24201, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 25619, 25615, 25616, 25648, 25840, 25841, 25852, 25860, 26072, 26101, 26103, 5, 24592, 24605, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 25844, 25840, 25841, 25863, 25935, 25936, 25941, 25949, 26106, 26116, 26118, 5, 25177, 25190, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 5517, 5478, 5493, 6077, 6089, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 5728, 5707, 5716, 6145, 6151, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6012, 5973, 5988, 6307, 6319, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6107, 6077, 6089, 6361, 6370, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6272, 6233, 6248, 6436, 6448, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6337, 6307, 6319, 6470, 6479, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6382, 6361, 6370, 6494, 6500, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 10961, 10916, 10931, 14137, 14152, 1, 5478, 5493, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 11576, 11546, 11558, 14240, 14249, 1, 5707, 5716, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13221, 13176, 13191, 14421, 14436, 1, 5973, 5988, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 14176, 14137, 14152, 14510, 14522, 1, 6077, 6089, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 14261, 14240, 14249, 14563, 14569, 1, 6145, 6151, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 14346, 14301, 14316, 14591, 14606, 1, 6233, 6248, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 14460, 14421, 14436, 14655, 14667, 1, 6307, 6319, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 14540, 14510, 14522, 14699, 14708, 1, 6361, 6370, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 14630, 14591, 14606, 14734, 14746, 1, 6436, 6448, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 14685, 14655, 14667, 14762, 14771, 1, 6470, 6479, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 16496, 16451, 16466, 19949, 19964, 2, 10916, 10931, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 17676, 17637, 17652, 20456, 20468, 2, 11546, 11558, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 17875, 17854, 17863, 20524, 20530, 2, 11614, 11620, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 17975, 17930, 17945, 20562, 20577, 2, 13176, 13191, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 19994, 19949, 19964, 20682, 20697, 2, 14137, 14152, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 20486, 20456, 20468, 20770, 20779, 2, 14240, 14249, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 20607, 20562, 20577, 20821, 20836, 2, 14421, 14436, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 20721, 20682, 20697, 20885, 20897, 2, 14510, 14522, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 20791, 20770, 20779, 20928, 20934, 2, 14563, 14569, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 20860, 20821, 20836, 20950, 20962, 2, 14655, 14667, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 20915, 20885, 20897, 20978, 20987, 2, 14699, 14708, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 21803, 21758, 21773, 23430, 23445, 3, 17637, 17652, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 22267, 22237, 22249, 23647, 23656, 3, 17854, 17863, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 22403, 22358, 22373, 23708, 23723, 3, 19949, 19964, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 23469, 23430, 23445, 23828, 23840, 3, 20456, 20468, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 23668, 23647, 23656, 23896, 23902, 3, 20524, 20530, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 23753, 23708, 23723, 23924, 23939, 3, 20682, 20697, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 23858, 23828, 23840, 23987, 23996, 3, 20770, 20779, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 23963, 23924, 23939, 24032, 24044, 3, 20885, 20897, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 24008, 23987, 23996, 24059, 24065, 3, 20928, 20934, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 24388, 24352, 24364, 25072, 25084, 4, 22237, 22249, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 24538, 24517, 24526, 25140, 25146, 4, 22305, 22311, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 24638, 24593, 24608, 25178, 25193, 4, 23430, 23445, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25102, 25072, 25084, 25281, 25290, 4, 23647, 23656, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25217, 25178, 25193, 25332, 25344, 4, 23828, 23840, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25302, 25281, 25290, 25385, 25391, 4, 23896, 23902, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25362, 25332, 25344, 25406, 25415, 4, 23987, 23996, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25565, 25538, 25547, 25804, 25813, 5, 24517, 24526, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25675, 25639, 25651, 25854, 25866, 5, 25072, 25084, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25825, 25804, 25813, 25922, 25928, 5, 25140, 25146, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25884, 25854, 25866, 25943, 25952, 5, 25281, 25290, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 25964, 25943, 25952, 25983, 25989, 5, 25385, 25391, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 26022, 26004, 26010, 26069, 26075, 6, 25804, 25813, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 26087, 26069, 26075, 26103, 26109, 6, 25922, 25928, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 5020, 4984, 4996, 5493, 5976, 5994, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 5527, 5481, 5499, 5716, 6080, 6095, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 5908, 5872, 5884, 5988, 6236, 6254, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 6022, 5976, 5994, 6089, 6310, 6325, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 6117, 6080, 6095, 6151, 6364, 6376, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 6282, 6236, 6254, 6319, 6439, 6454, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 6347, 6310, 6325, 6370, 6473, 6485, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 9715, 9679, 9691, 10931, 13179, 13197, 1, 4984, 4996, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 10971, 10919, 10937, 11558, 14140, 14158, 1, 5481, 5499, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 11586, 11549, 11564, 11620, 14243, 14255, 1, 5710, 5722, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 11702, 11666, 11678, 13191, 14304, 14322, 1, 5872, 5884, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 13231, 13179, 13197, 14152, 14424, 14442, 1, 5976, 5994, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 14186, 14140, 14158, 14249, 14513, 14528, 1, 6080, 6095, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 14356, 14304, 14322, 14436, 14594, 14612, 1, 6236, 6254, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 14470, 14424, 14442, 14522, 14658, 14673, 1, 6310, 6325, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 14550, 14513, 14528, 14569, 14702, 14714, 1, 6364, 6376, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 14640, 14594, 14612, 14667, 14737, 14752, 1, 6439, 6454, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 16506, 16454, 16472, 17652, 19952, 19970, 2, 10919, 10937, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 17686, 17640, 17658, 17863, 20459, 20474, 2, 11549, 11564, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 17985, 17933, 17951, 19964, 20565, 20583, 2, 13179, 13197, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 20004, 19952, 19970, 20468, 20685, 20703, 2, 14140, 14158, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 20496, 20459, 20474, 20530, 20773, 20785, 2, 14243, 14255, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 20617, 20565, 20583, 20697, 20824, 20842, 2, 14424, 14442, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 20731, 20685, 20703, 20779, 20888, 20903, 2, 14513, 14528, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 20870, 20824, 20842, 20897, 20953, 20968, 2, 14658, 14673, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 21813, 21761, 21779, 22249, 23433, 23451, 3, 17640, 17658, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 22277, 22240, 22255, 22311, 23650, 23662, 3, 17857, 17869, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 22413, 22361, 22379, 23445, 23711, 23729, 3, 19952, 19970, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 23479, 23433, 23451, 23656, 23831, 23846, 3, 20459, 20474, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 23763, 23711, 23729, 23840, 23927, 23945, 3, 20685, 20703, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 23868, 23831, 23846, 23902, 23990, 24002, 3, 20773, 20785, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 23973, 23927, 23945, 23996, 24035, 24050, 3, 20888, 20903, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 24398, 24355, 24370, 24526, 25075, 25090, 4, 22240, 22255, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 24648, 24596, 24614, 25084, 25181, 25199, 4, 23433, 23451, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 25112, 25075, 25090, 25146, 25284, 25296, 4, 23650, 23662, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 25227, 25181, 25199, 25290, 25335, 25350, 4, 23831, 23846, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 25372, 25335, 25350, 25391, 25409, 25421, 4, 23990, 24002, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 25575, 25541, 25553, 25608, 25807, 25819, 5, 24520, 24532, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 25685, 25642, 25657, 25813, 25857, 25872, 5, 25075, 25090, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 25894, 25857, 25872, 25928, 25946, 25958, 5, 25284, 25296, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 26032, 26007, 26016, 26061, 26072, 26081, 6, 25807, 25819, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4094, 4067, 4076, 4996, 5875, 5890, 5973, 5988, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 5030, 4987, 5002, 5499, 5979, 6000, 6077, 6089, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 5537, 5484, 5505, 5722, 6083, 6101, 6145, 6151, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 5808, 5781, 5790, 5884, 6181, 6190, 6233, 6248, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 5918, 5875, 5890, 5994, 6239, 6260, 6307, 6319, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 6032, 5979, 6000, 6095, 6313, 6331, 6361, 6370, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 6208, 6181, 6190, 6254, 6410, 6419, 6436, 6448, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 6292, 6239, 6260, 6325, 6442, 6460, 6470, 6479, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 7741, 7714, 7723, 9691, 11669, 11684, 13176, 13191, 1, 4067, 4076, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 9725, 9682, 9697, 10937, 13182, 13203, 14137, 14152, 1, 4987, 5002, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 10981, 10922, 10943, 11564, 14143, 14164, 14240, 14249, 1, 5484, 5505, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 11712, 11669, 11684, 13197, 14307, 14328, 14421, 14436, 1, 5875, 5890, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 13241, 13182, 13203, 14158, 14427, 14448, 14510, 14522, 1, 5979, 6000, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 14196, 14143, 14164, 14255, 14516, 14534, 14563, 14569, 1, 6083, 6101, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 14366, 14307, 14328, 14442, 14597, 14618, 14655, 14667, 1, 6239, 6260, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 14480, 14427, 14448, 14528, 14661, 14679, 14699, 14708, 1, 6313, 6331, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 15459, 15432, 15441, 16472, 17936, 17957, 19949, 19964, 2, 9682, 9697, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 16516, 16457, 16478, 17658, 19955, 19976, 20456, 20468, 2, 10922, 10943, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 17696, 17643, 17664, 17869, 20462, 20480, 20524, 20530, 2, 11552, 11570, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 17995, 17936, 17957, 19970, 20568, 20589, 20682, 20697, 2, 13182, 13203, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 20014, 19955, 19976, 20474, 20688, 20709, 20770, 20779, 2, 14143, 14164, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 20627, 20568, 20589, 20703, 20827, 20848, 20885, 20897, 2, 14427, 14448, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 20741, 20688, 20709, 20785, 20891, 20909, 20928, 20934, 2, 14516, 14534, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 21333, 21306, 21315, 21779, 22364, 22385, 23430, 23445, 3, 16457, 16478, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 21823, 21764, 21785, 22255, 23436, 23457, 23647, 23656, 3, 17643, 17664, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 22423, 22364, 22385, 23451, 23714, 23735, 23828, 23840, 3, 19955, 19976, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 23489, 23436, 23457, 23662, 23834, 23852, 23896, 23902, 3, 20462, 20480, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 23773, 23714, 23735, 23846, 23930, 23951, 23987, 23996, 3, 20688, 20709, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 24222, 24195, 24204, 24370, 24599, 24620, 25072, 25084, 4, 21764, 21785, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 24408, 24358, 24376, 24532, 25078, 25096, 25140, 25146, 4, 22243, 22261, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 24658, 24599, 24620, 25090, 25184, 25205, 25281, 25290, 4, 23436, 23457, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 25237, 25184, 25205, 25296, 25338, 25356, 25385, 25391, 4, 23834, 23852, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 25500, 25482, 25488, 25553, 25645, 25663, 25804, 25813, 5, 24358, 24376, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 25695, 25645, 25663, 25819, 25860, 25878, 25922, 25928, 5, 25078, 25096, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 2630, 2612, 2618, 4076, 5784, 5796, 5872, 5884, 13191, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 4104, 4070, 4082, 5002, 5878, 5896, 5976, 5994, 14152, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 5040, 4990, 5008, 5505, 5982, 6006, 6080, 6095, 14249, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 5818, 5784, 5796, 5890, 6184, 6196, 6236, 6254, 14436, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 5928, 5878, 5896, 6000, 6242, 6266, 6310, 6325, 14522, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 6218, 6184, 6196, 6260, 6413, 6425, 6439, 6454, 14667, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 7751, 7717, 7729, 9697, 11672, 11690, 13179, 13197, 19964, 1, 4070, 4082, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 9735, 9685, 9703, 10943, 13185, 13209, 14140, 14158, 20468, 1, 4990, 5008, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 10991, 10925, 10949, 11570, 14146, 14170, 14243, 14255, 20530, 1, 5487, 5511, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 11722, 11672, 11690, 13203, 14310, 14334, 14424, 14442, 20697, 1, 5878, 5896, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 13251, 13185, 13209, 14164, 14430, 14454, 14513, 14528, 20779, 1, 5982, 6006, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 14376, 14310, 14334, 14448, 14600, 14624, 14658, 14673, 20897, 1, 6242, 6266, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 15469, 15435, 15447, 16478, 17939, 17963, 19952, 19970, 23445, 2, 9685, 9703, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 16526, 16460, 16484, 17664, 19958, 19982, 20459, 20474, 23656, 2, 10925, 10949, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 18005, 17939, 17963, 19976, 20571, 20595, 20685, 20703, 23840, 2, 13185, 13209, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 20024, 19958, 19982, 20480, 20691, 20715, 20773, 20785, 23902, 2, 14146, 14170, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 20637, 20571, 20595, 20709, 20830, 20854, 20888, 20903, 23996, 2, 14430, 14454, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 21343, 21309, 21321, 21785, 22367, 22391, 23433, 23451, 25084, 3, 16460, 16484, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 21833, 21767, 21791, 22261, 23439, 23463, 23650, 23662, 25146, 3, 17646, 17670, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 22433, 22367, 22391, 23457, 23717, 23741, 23831, 23846, 25290, 3, 19958, 19982, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 23783, 23717, 23741, 23852, 23933, 23957, 23990, 24002, 25391, 3, 20691, 20715, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 24232, 24198, 24210, 24376, 24602, 24626, 25075, 25090, 25813, 4, 21767, 21791, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 24668, 24602, 24626, 25096, 25187, 25211, 25284, 25296, 25928, 4, 23439, 23463, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 25510, 25485, 25494, 25559, 25648, 25669, 25807, 25819, 26061, 5, 24361, 24382, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 10, 1, 4, 2618, 5744, 5747, 5781, 5790, 11678, 14301, 14316, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 2640, 2615, 2624, 4082, 5787, 5802, 5875, 5890, 13197, 14421, 14436, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 4114, 4073, 4088, 5008, 5881, 5902, 5979, 6000, 14158, 14510, 14522, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 5753, 5744, 5747, 5796, 6159, 6162, 6181, 6190, 14322, 14591, 14606, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 5828, 5787, 5802, 5896, 6187, 6202, 6239, 6260, 14442, 14655, 14667, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 6168, 6159, 6162, 6196, 6398, 6401, 6410, 6419, 14612, 14734, 14746, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 6566, 6557, 6560, 7729, 11628, 11631, 11669, 11684, 17951, 20562, 20577, 1, 2615, 2624, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 7761, 7720, 7735, 9703, 11675, 11696, 13182, 13203, 19970, 20682, 20697, 1, 4073, 4088, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 9745, 9688, 9709, 10949, 13188, 13215, 14143, 14164, 20474, 20770, 20779, 1, 4993, 5014, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 11637, 11628, 11631, 11690, 14277, 14280, 14307, 14328, 20583, 20821, 20836, 1, 5787, 5802, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 11732, 11675, 11696, 13209, 14313, 14340, 14427, 14448, 20703, 20885, 20897, 1, 5881, 5902, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 14286, 14277, 14280, 14334, 14577, 14580, 14597, 14618, 20842, 20950, 20962, 1, 6187, 6202, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 14819, 14810, 14813, 15447, 17891, 17894, 17936, 17957, 22379, 23708, 23723, 2, 7720, 7735, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 15479, 15438, 15453, 16484, 17942, 17969, 19955, 19976, 23451, 23828, 23840, 2, 9688, 9709, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 16536, 16463, 16490, 17670, 19961, 19988, 20462, 20480, 23662, 23896, 23902, 2, 10928, 10955, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 17900, 17891, 17894, 17963, 20538, 20541, 20568, 20589, 23729, 23924, 23939, 2, 11675, 11696, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 18015, 17942, 17969, 19982, 20574, 20601, 20688, 20709, 23846, 23987, 23996, 2, 13188, 13215, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 20547, 20538, 20541, 20595, 20807, 20810, 20827, 20848, 23945, 24032, 24044, 2, 14313, 14340, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 21029, 21020, 21023, 21321, 22319, 22322, 22364, 22385, 24614, 25178, 25193, 3, 15438, 15453, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 21353, 21312, 21327, 21791, 22370, 22397, 23436, 23457, 25090, 25281, 25290, 3, 16463, 16490, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 22328, 22319, 22322, 22391, 23684, 23687, 23714, 23735, 25199, 25332, 25344, 3, 17942, 17969, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 22443, 22370, 22397, 23463, 23720, 23747, 23834, 23852, 25296, 25385, 25391, 3, 19961, 19988, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 23693, 23684, 23687, 23741, 23910, 23913, 23930, 23951, 25350, 25406, 25415, 3, 20574, 20601, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 24101, 24092, 24095, 24210, 24554, 24557, 24599, 24620, 25657, 25854, 25866, 4, 21312, 21327, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 24242, 24201, 24216, 24382, 24605, 24632, 25078, 25096, 25819, 25922, 25928, 4, 21770, 21797, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 24563, 24554, 24557, 24626, 25154, 25157, 25184, 25205, 25872, 25943, 25952, 4, 22370, 22397, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 25163, 25154, 25157, 25211, 25318, 25321, 25338, 25356, 25958, 25983, 25989, 4, 23720, 23747, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 25455, 25446, 25449, 25494, 25616, 25619, 25645, 25663, 26016, 26069, 26075, 5, 24201, 24216, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 25625, 25616, 25619, 25669, 25841, 25844, 25860, 25878, 26081, 26103, 26109, 5, 24605, 24632, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 5547, 5493, 5517, 6089, 6107, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 6042, 5988, 6012, 6319, 6337, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 6127, 6089, 6107, 6370, 6382, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 11001, 10931, 10961, 14152, 14176, 1, 5493, 5517, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 11596, 11558, 11576, 14249, 14261, 1, 5716, 5728, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 13261, 13191, 13221, 14436, 14460, 1, 5988, 6012, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 14206, 14152, 14176, 14522, 14540, 1, 6089, 6107, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 14386, 14316, 14346, 14606, 14630, 1, 6248, 6272, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 14490, 14436, 14460, 14667, 14685, 1, 6319, 6337, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 16546, 16466, 16496, 19964, 19994, 2, 10931, 10961, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 17706, 17652, 17676, 20468, 20486, 2, 11558, 11576, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 18025, 17945, 17975, 20577, 20607, 2, 13191, 13221, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 20034, 19964, 19994, 20697, 20721, 2, 14152, 14176, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 20506, 20468, 20486, 20779, 20791, 2, 14249, 14261, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 20647, 20577, 20607, 20836, 20860, 2, 14436, 14460, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 20751, 20697, 20721, 20897, 20915, 2, 14522, 14540, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 21843, 21773, 21803, 23445, 23469, 3, 17652, 17676, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 22287, 22249, 22267, 23656, 23668, 3, 17863, 17875, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 22453, 22373, 22403, 23723, 23753, 3, 19964, 19994, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 23499, 23445, 23469, 23840, 23858, 3, 20468, 20486, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 23793, 23723, 23753, 23939, 23963, 3, 20697, 20721, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 23878, 23840, 23858, 23996, 24008, 3, 20779, 20791, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 24418, 24364, 24388, 25084, 25102, 4, 22249, 22267, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 24678, 24608, 24638, 25193, 25217, 4, 23445, 23469, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 25122, 25084, 25102, 25290, 25302, 4, 23656, 23668, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 25247, 25193, 25217, 25344, 25362, 4, 23840, 23858, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 25585, 25547, 25565, 25813, 25825, 5, 24526, 24538, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 25705, 25651, 25675, 25866, 25884, 5, 25084, 25102, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 25904, 25866, 25884, 25952, 25964, 5, 25290, 25302, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 26042, 26010, 26022, 26075, 26087, 6, 25813, 25825, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 5050, 4996, 5020, 5517, 5994, 6022, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 5562, 5499, 5527, 5728, 6095, 6117, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 5938, 5884, 5908, 6012, 6254, 6282, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 6057, 5994, 6022, 6107, 6325, 6347, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 9755, 9691, 9715, 10961, 13197, 13231, 1, 4996, 5020, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 11016, 10937, 10971, 11576, 14158, 14186, 1, 5499, 5527, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 11742, 11678, 11702, 13221, 14322, 14356, 1, 5884, 5908, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 13276, 13197, 13231, 14176, 14442, 14470, 1, 5994, 6022, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 14221, 14158, 14186, 14261, 14528, 14550, 1, 6095, 6117, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 14401, 14322, 14356, 14460, 14612, 14640, 1, 6254, 6282, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 16561, 16472, 16506, 17676, 19970, 20004, 2, 10937, 10971, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 17721, 17658, 17686, 17875, 20474, 20496, 2, 11564, 11586, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 18040, 17951, 17985, 19994, 20583, 20617, 2, 13197, 13231, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 20049, 19970, 20004, 20486, 20703, 20731, 2, 14158, 14186, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 20662, 20583, 20617, 20721, 20842, 20870, 2, 14442, 14470, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 21858, 21779, 21813, 22267, 23451, 23479, 3, 17658, 17686, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 22468, 22379, 22413, 23469, 23729, 23763, 3, 19970, 20004, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 23514, 23451, 23479, 23668, 23846, 23868, 3, 20474, 20496, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 23808, 23729, 23763, 23858, 23945, 23973, 3, 20703, 20731, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 24433, 24370, 24398, 24538, 25090, 25112, 4, 22255, 22277, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 24693, 24614, 24648, 25102, 25199, 25227, 4, 23451, 23479, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 25262, 25199, 25227, 25302, 25350, 25372, 4, 23846, 23868, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 25720, 25657, 25685, 25825, 25872, 25894, 5, 25090, 25112, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 4124, 4076, 4094, 5020, 5890, 5918, 5988, 6012, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 5065, 5002, 5030, 5527, 6000, 6032, 6089, 6107, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 5838, 5790, 5808, 5908, 6190, 6208, 6248, 6272, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 5953, 5890, 5918, 6022, 6260, 6292, 6319, 6337, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 7771, 7723, 7741, 9715, 11684, 11712, 13191, 13221, 1, 4076, 4094, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 9770, 9697, 9725, 10971, 13203, 13241, 14152, 14176, 1, 5002, 5030, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 11031, 10943, 10981, 11586, 14164, 14196, 14249, 14261, 1, 5505, 5537, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 11757, 11684, 11712, 13231, 14328, 14366, 14436, 14460, 1, 5890, 5918, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 13291, 13203, 13241, 14186, 14448, 14480, 14522, 14540, 1, 6000, 6032, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 15489, 15441, 15459, 16506, 17957, 17995, 19964, 19994, 2, 9697, 9725, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 16576, 16478, 16516, 17686, 19976, 20014, 20468, 20486, 2, 10943, 10981, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 18055, 17957, 17995, 20004, 20589, 20627, 20697, 20721, 2, 13203, 13241, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 20064, 19976, 20014, 20496, 20709, 20741, 20779, 20791, 2, 14164, 14196, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 21363, 21315, 21333, 21813, 22385, 22423, 23445, 23469, 3, 16478, 16516, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 21873, 21785, 21823, 22277, 23457, 23489, 23656, 23668, 3, 17664, 17696, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 22483, 22385, 22423, 23479, 23735, 23773, 23840, 23858, 3, 19976, 20014, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 24252, 24204, 24222, 24398, 24620, 24658, 25084, 25102, 4, 21785, 21823, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 24708, 24620, 24658, 25112, 25205, 25237, 25290, 25302, 4, 23457, 23489, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 25520, 25488, 25500, 25575, 25663, 25695, 25813, 25825, 5, 24376, 24408, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 2650, 2618, 2630, 4094, 5796, 5818, 5884, 5908, 13221, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 4139, 4082, 4104, 5030, 5896, 5928, 5994, 6022, 14176, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 5853, 5796, 5818, 5918, 6196, 6218, 6254, 6282, 14460, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 7786, 7729, 7751, 9725, 11690, 11722, 13197, 13231, 19994, 1, 4082, 4104, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 9785, 9703, 9735, 10981, 13209, 13251, 14158, 14186, 20486, 1, 5008, 5040, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 11772, 11690, 11722, 13241, 14334, 14376, 14442, 14470, 20721, 1, 5896, 5928, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 15504, 15447, 15469, 16516, 17963, 18005, 19970, 20004, 23469, 2, 9703, 9735, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 16591, 16484, 16526, 17696, 19982, 20024, 20474, 20496, 23668, 2, 10949, 10991, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 18070, 17963, 18005, 20014, 20595, 20637, 20703, 20731, 23858, 2, 13209, 13251, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 21378, 21321, 21343, 21823, 22391, 22433, 23451, 23479, 25102, 3, 16484, 16526, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 22498, 22391, 22433, 23489, 23741, 23783, 23846, 23868, 25302, 3, 19982, 20024, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 24267, 24210, 24232, 24408, 24626, 24668, 25090, 25112, 25825, 4, 21791, 21833, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 20, 4, 10, 2630, 5747, 5753, 5790, 5808, 11702, 14316, 14346, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 2665, 2624, 2640, 4104, 5802, 5828, 5890, 5918, 13231, 14436, 14460, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 5763, 5747, 5753, 5818, 6162, 6168, 6190, 6208, 14356, 14606, 14630, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 6576, 6560, 6566, 7751, 11631, 11637, 11684, 11712, 17985, 20577, 20607, 1, 2624, 2640, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 7801, 7735, 7761, 9735, 11696, 11732, 13203, 13241, 20004, 20697, 20721, 1, 4088, 4114, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 11647, 11631, 11637, 11722, 14280, 14286, 14328, 14366, 20617, 20836, 20860, 1, 5802, 5828, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 14829, 14813, 14819, 15469, 17894, 17900, 17957, 17995, 22413, 23723, 23753, 2, 7735, 7761, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 15519, 15453, 15479, 16526, 17969, 18015, 19976, 20014, 23479, 23840, 23858, 2, 9709, 9745, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 17910, 17894, 17900, 18005, 20541, 20547, 20589, 20627, 23763, 23939, 23963, 2, 11696, 11732, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 21039, 21023, 21029, 21343, 22322, 22328, 22385, 22423, 24648, 25193, 25217, 3, 15453, 15479, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 21393, 21327, 21353, 21833, 22397, 22443, 23457, 23489, 25112, 25290, 25302, 3, 16490, 16536, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 22338, 22322, 22328, 22433, 23687, 23693, 23735, 23773, 25227, 25344, 25362, 3, 17969, 18015, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 24111, 24095, 24101, 24232, 24557, 24563, 24620, 24658, 25685, 25866, 25884, 4, 21327, 21353, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 24573, 24557, 24563, 24668, 25157, 25163, 25205, 25237, 25894, 25952, 25964, 4, 22397, 22443, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 25465, 25449, 25455, 25510, 25619, 25625, 25663, 25695, 26032, 26075, 26087, 5, 24216, 24242, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 5577, 5517, 5547, 6107, 6127, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 11046, 10961, 11001, 14176, 14206, 1, 5517, 5547, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 13306, 13221, 13261, 14460, 14490, 1, 6012, 6042, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 16606, 16496, 16546, 19994, 20034, 2, 10961, 11001, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 17736, 17676, 17706, 20486, 20506, 2, 11576, 11596, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 18085, 17975, 18025, 20607, 20647, 2, 13221, 13261, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 20079, 19994, 20034, 20721, 20751, 2, 14176, 14206, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 21888, 21803, 21843, 23469, 23499, 3, 17676, 17706, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 22513, 22403, 22453, 23753, 23793, 3, 19994, 20034, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 23529, 23469, 23499, 23858, 23878, 3, 20486, 20506, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 24448, 24388, 24418, 25102, 25122, 4, 22267, 22287, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 24723, 24638, 24678, 25217, 25247, 4, 23469, 23499, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 25735, 25675, 25705, 25884, 25904, 5, 25102, 25122, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 5080, 5020, 5050, 5547, 6022, 6057, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 9800, 9715, 9755, 11001, 13231, 13276, 1, 5020, 5050, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 11067, 10971, 11016, 11596, 14186, 14221, 1, 5527, 5562, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 11787, 11702, 11742, 13261, 14356, 14401, 1, 5908, 5938, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 16627, 16506, 16561, 17706, 20004, 20049, 2, 10971, 11016, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 18106, 17985, 18040, 20034, 20617, 20662, 2, 13231, 13276, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 21909, 21813, 21858, 22287, 23479, 23514, 3, 17686, 17721, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 22534, 22413, 22468, 23499, 23763, 23808, 3, 20004, 20049, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 24744, 24648, 24693, 25122, 25227, 25262, 4, 23479, 23514, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 4154, 4094, 4124, 5050, 5918, 5953, 6012, 6042, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 7816, 7741, 7771, 9755, 11712, 11757, 13221, 13261, 1, 4094, 4124, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 9821, 9725, 9770, 11016, 13241, 13291, 14176, 14206, 1, 5030, 5065, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 15534, 15459, 15489, 16561, 17995, 18055, 19994, 20034, 2, 9725, 9770, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 16648, 16516, 16576, 17721, 20014, 20064, 20486, 20506, 2, 10981, 11031, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 21408, 21333, 21363, 21858, 22423, 22483, 23469, 23499, 3, 16516, 16576, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 24282, 24222, 24252, 24433, 24658, 24708, 25102, 25122, 4, 21823, 21873, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 2680, 2630, 2650, 4124, 5818, 5853, 5908, 5938, 13261, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 7837, 7751, 7786, 9770, 11722, 11772, 13231, 13276, 20034, 1, 4104, 4139, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 15555, 15469, 15504, 16576, 18005, 18070, 20004, 20049, 23499, 2, 9735, 9785, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 21429, 21343, 21378, 21873, 22433, 22498, 23479, 23514, 25122, 3, 16526, 16591, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_g(pbuffer, 35, 10, 20, 2650, 5753, 5763, 5808, 5838, 11742, 14346, 14386, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_g(pbuffer, 6591, 6566, 6576, 7786, 11637, 11647, 11712, 11757, 18040, 20607, 20647, 1, 2640, 2665, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_g(pbuffer, 14844, 14819, 14829, 15504, 17900, 17910, 17995, 18055, 22468, 23753, 23793, 2, 7761, 7801, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_g(pbuffer, 21054, 21029, 21039, 21378, 22328, 22338, 22423, 22483, 24693, 25217, 25247, 3, 15479, 15519, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_g(pbuffer, 24126, 24101, 24111, 24267, 24563, 24573, 24658, 24708, 25720, 25884, 25904, 4, 21353, 21393, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_s(pbuffer, 5598, 5473, 11541, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_p(pbuffer, 5101, 4980, 5473, 10912, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 5601, 5478, 11546, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 11088, 10916, 17637, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 13327, 13176, 19949, 1, 10916, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 20100, 19949, 23430, 1, 17637, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 5104, 4984, 5473, 5478, 10919, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 11097, 10919, 11541, 11546, 17640, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_d(pbuffer, 4175, 4067, 4980, 4984, 9682, 10916, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_d(pbuffer, 9842, 9682, 10912, 10919, 16457, 17637, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 5610, 5493, 11558, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 11106, 10931, 17652, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 13336, 13191, 19964, 1, 10931, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 17757, 17652, 22249, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 20109, 19964, 23445, 1, 17652, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 23550, 23445, 25084, 1, 22249, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 5113, 4996, 5478, 5493, 10937, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 9851, 9691, 10916, 10931, 16472, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 11124, 10937, 11546, 11558, 17658, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 11808, 11678, 13176, 13191, 17951, 1, 9691, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 16669, 16472, 17637, 17652, 21779, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 18127, 17951, 19949, 19964, 22379, 1, 16472, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 22555, 22379, 23430, 23445, 24614, 1, 21779, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 4184, 4076, 4984, 4996, 9697, 10931, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 9869, 9697, 10919, 10937, 16478, 17652, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 16687, 16478, 17640, 17658, 21785, 22249, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 2701, 2618, 4067, 4076, 7729, 9691, 13176, 13191, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 7858, 7729, 9682, 9697, 15447, 16472, 19949, 19964, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 15576, 15447, 16457, 16478, 21321, 21779, 23430, 23445, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 5628, 5517, 11576, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 11142, 10961, 17676, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 13354, 13221, 19994, 1, 10961, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 16705, 16496, 21803, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 17775, 17676, 22267, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 18145, 17975, 22403, 1, 16496, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 20127, 19994, 23469, 1, 17676, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 21930, 21803, 24388, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 22573, 22403, 24638, 1, 21803, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 23568, 23469, 25102, 1, 22267, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 24765, 24638, 25675, 1, 24388, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 5131, 5020, 5493, 5517, 10971, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 9887, 9715, 10931, 10961, 16506, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 11172, 10971, 11558, 11576, 17686, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 11826, 11702, 13191, 13221, 17985, 1, 9715, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 16735, 16506, 17652, 17676, 21813, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 18175, 17985, 19964, 19994, 22413, 1, 16506, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 21960, 21813, 22249, 22267, 24398, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 22603, 22413, 23445, 23469, 24648, 1, 21813, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 24795, 24648, 25084, 25102, 25685, 1, 24398, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 4202, 4094, 4996, 5020, 9725, 10961, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 7876, 7741, 9691, 9715, 15459, 16496, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 9917, 9725, 10937, 10971, 16516, 17676, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 15594, 15459, 16472, 16506, 21333, 21803, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 16765, 16516, 17658, 17686, 21823, 22267, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 21450, 21333, 21779, 21813, 24222, 24388, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 2719, 2630, 4076, 4094, 7751, 9715, 13191, 13221, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 7906, 7751, 9697, 9725, 15469, 16506, 19964, 19994, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 15624, 15469, 16478, 16516, 21343, 21813, 23445, 23469, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 21480, 21343, 21785, 21823, 24232, 24398, 25084, 25102, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 56, 10, 2618, 2630, 6566, 7741, 11678, 11702, 17975, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 6612, 6566, 7729, 7751, 14819, 15459, 17951, 17985, 22403, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 14865, 14819, 15447, 15469, 21029, 21333, 22379, 22413, 24638, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 21075, 21029, 21321, 21343, 24101, 24222, 24614, 24648, 25675, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 5658, 5547, 11596, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 11202, 11001, 17706, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 13384, 13261, 20034, 1, 11001, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 16795, 16546, 21843, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 17805, 17706, 22287, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 18205, 18025, 22453, 1, 16546, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 20157, 20034, 23499, 1, 17706, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 21990, 21843, 24418, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 22633, 22453, 24678, 1, 21843, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 23598, 23499, 25122, 1, 22287, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 24469, 24418, 25585, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 24825, 24678, 25705, 1, 24418, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 25756, 25705, 26042, 1, 25585, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 5161, 5050, 5517, 5547, 11016, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 9947, 9755, 10961, 11001, 16561, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 11247, 11016, 11576, 11596, 17721, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 11856, 11742, 13221, 13261, 18040, 1, 9755, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 16840, 16561, 17676, 17706, 21858, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 18250, 18040, 19994, 20034, 22468, 1, 16561, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 22035, 21858, 22267, 22287, 24433, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 22678, 22468, 23469, 23499, 24693, 1, 21858, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 24870, 24693, 25102, 25122, 25720, 1, 24433, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 4232, 4124, 5020, 5050, 9770, 11001, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 7936, 7771, 9715, 9755, 15489, 16546, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 9992, 9770, 10971, 11016, 16576, 17706, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 15654, 15489, 16506, 16561, 21363, 21843, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 16885, 16576, 17686, 17721, 21873, 22287, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 21510, 21363, 21813, 21858, 24252, 24418, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 24303, 24252, 24398, 24433, 25520, 25585, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 2749, 2650, 4094, 4124, 7786, 9755, 13221, 13261, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 7981, 7786, 9725, 9770, 15504, 16561, 19994, 20034, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 15699, 15504, 16516, 16576, 21378, 21858, 23469, 23499, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 21555, 21378, 21823, 21873, 24267, 24433, 25102, 25122, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 86, 20, 2630, 2650, 6576, 7771, 11702, 11742, 18025, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 6642, 6576, 7751, 7786, 14829, 15489, 17985, 18040, 22453, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 14895, 14829, 15469, 15504, 21039, 21363, 22413, 22468, 24678, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 21105, 21039, 21343, 21378, 24111, 24252, 24648, 24693, 25705, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 24147, 24111, 24232, 24267, 25465, 25520, 25685, 25720, 26042, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 11292, 11046, 17736, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 13429, 13306, 20079, 1, 11046, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 16930, 16606, 21888, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 18295, 18085, 22513, 1, 16606, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 20202, 20079, 23529, 1, 17736, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 22080, 21888, 24448, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 22723, 22513, 24723, 1, 21888, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 24915, 24723, 25735, 1, 24448, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 5206, 5080, 5547, 5577, 11067, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 10037, 9800, 11001, 11046, 16627, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 11901, 11787, 13261, 13306, 18106, 1, 9800, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 16993, 16627, 17706, 17736, 21909, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 18358, 18106, 20034, 20079, 22534, 1, 16627, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 22786, 22534, 23499, 23529, 24744, 1, 21909, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 4277, 4154, 5050, 5080, 9821, 11046, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 8026, 7816, 9755, 9800, 15534, 16606, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 10100, 9821, 11016, 11067, 16648, 17736, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 15744, 15534, 16561, 16627, 21408, 21888, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 21600, 21408, 21858, 21909, 24282, 24448, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 2794, 2680, 4124, 4154, 7837, 9800, 13261, 13306, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 8089, 7837, 9770, 9821, 15555, 16627, 20034, 20079, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 15807, 15555, 16576, 16648, 21429, 21909, 23499, 23529, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_g(pbuffer, 131, 35, 2650, 2680, 6591, 7816, 11742, 11787, 18085, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_g(pbuffer, 6687, 6591, 7786, 7837, 14844, 15534, 18040, 18106, 22513, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_g(pbuffer, 14940, 14844, 15504, 15555, 21054, 21408, 22468, 22534, 24723, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_g(pbuffer, 21150, 21054, 21378, 21429, 24126, 24282, 24693, 24744, 25735, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_s(pbuffer, 13492, 13176, 13327, 19949, 20100, 1, 10916, 11088, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_p(pbuffer, 5269, 4984, 5104, 5598, 5601, 10919, 11097, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_d(pbuffer, 4340, 4067, 4175, 5101, 5104, 9682, 9842, 10916, 11088, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 11355, 10931, 11106, 17652, 17757, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 13510, 13191, 13336, 19964, 20109, 1, 10931, 11106, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 20265, 19964, 20109, 23445, 23550, 1, 17652, 17757, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 5287, 4996, 5113, 5601, 5610, 10937, 11124, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 10163, 9691, 9851, 11088, 11106, 16472, 16669, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 11964, 11678, 11808, 13327, 13336, 17951, 18127, 1, 9691, 9851, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 18421, 17951, 18127, 20100, 20109, 22379, 22555, 1, 16472, 16669, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_d(pbuffer, 4358, 4076, 4184, 5104, 5113, 9697, 9869, 10931, 11106, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_d(pbuffer, 10199, 9697, 9869, 11097, 11124, 16478, 16687, 17652, 17757, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_f(pbuffer, 2857, 2618, 2701, 4175, 4184, 7729, 7858, 9691, 9851, 13327, 13336, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_f(pbuffer, 8152, 7729, 7858, 9842, 9869, 15447, 15576, 16472, 16669, 20100, 20109, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 11391, 10961, 11142, 17676, 17775, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 13546, 13221, 13354, 19994, 20127, 1, 10961, 11142, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 17056, 16496, 16705, 21803, 21930, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 18457, 17975, 18145, 22403, 22573, 1, 16496, 16705, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 20301, 19994, 20127, 23469, 23568, 1, 17676, 17775, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 22849, 22403, 22573, 24638, 24765, 1, 21803, 21930, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 5323, 5020, 5131, 5610, 5628, 10971, 11172, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 10235, 9715, 9887, 11106, 11142, 16506, 16735, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 12000, 11702, 11826, 13336, 13354, 17985, 18175, 1, 9715, 9887, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 17116, 16506, 16735, 17757, 17775, 21813, 21960, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 18517, 17985, 18175, 20109, 20127, 22413, 22603, 1, 16506, 16735, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 22909, 22413, 22603, 23550, 23568, 24648, 24795, 1, 21813, 21960, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 4394, 4094, 4202, 5113, 5131, 9725, 9917, 10961, 11142, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 8188, 7741, 7876, 9851, 9887, 15459, 15594, 16496, 16705, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 10295, 9725, 9917, 11124, 11172, 16516, 16765, 17676, 17775, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 15870, 15459, 15594, 16669, 16735, 21333, 21450, 21803, 21930, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 2893, 2630, 2719, 4184, 4202, 7751, 7906, 9715, 9887, 13336, 13354, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 8248, 7751, 7906, 9869, 9917, 15469, 15624, 16506, 16735, 20109, 20127, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 15930, 15469, 15624, 16687, 16765, 21343, 21480, 21813, 21960, 23550, 23568, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_g(pbuffer, 194, 10, 56, 2701, 2719, 6566, 6612, 7741, 7876, 11808, 11826, 17975, 18145, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_g(pbuffer, 6750, 6566, 6612, 7858, 7906, 14819, 14865, 15459, 15594, 18127, 18175, 22403, 22573, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_g(pbuffer, 15003, 14819, 14865, 15576, 15624, 21029, 21075, 21333, 21450, 22555, 22603, 24638, 24765, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 11451, 11001, 11202, 17706, 17805, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 13606, 13261, 13384, 20034, 20157, 1, 11001, 11202, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 17176, 16546, 16795, 21843, 21990, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 18577, 18025, 18205, 22453, 22633, 1, 16546, 16795, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 20361, 20034, 20157, 23499, 23598, 1, 17706, 17805, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 22143, 21843, 21990, 24418, 24469, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 22969, 22453, 22633, 24678, 24825, 1, 21843, 21990, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 24978, 24678, 24825, 25705, 25756, 1, 24418, 24469, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 5383, 5050, 5161, 5628, 5658, 11016, 11247, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 10355, 9755, 9947, 11142, 11202, 16561, 16840, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 12060, 11742, 11856, 13354, 13384, 18040, 18250, 1, 9755, 9947, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 17266, 16561, 16840, 17775, 17805, 21858, 22035, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 18667, 18040, 18250, 20127, 20157, 22468, 22678, 1, 16561, 16840, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 23059, 22468, 22678, 23568, 23598, 24693, 24870, 1, 21858, 22035, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 4454, 4124, 4232, 5131, 5161, 9770, 9992, 11001, 11202, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 8308, 7771, 7936, 9887, 9947, 15489, 15654, 16546, 16795, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 10445, 9770, 9992, 11172, 11247, 16576, 16885, 17706, 17805, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 15990, 15489, 15654, 16735, 16840, 21363, 21510, 21843, 21990, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 21663, 21363, 21510, 21960, 22035, 24252, 24303, 24418, 24469, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_f(pbuffer, 2953, 2650, 2749, 4202, 4232, 7786, 7981, 9755, 9947, 13354, 13384, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_f(pbuffer, 8398, 7786, 7981, 9917, 9992, 15504, 15699, 16561, 16840, 20127, 20157, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_f(pbuffer, 16080, 15504, 15699, 16765, 16885, 21378, 21555, 21858, 22035, 23568, 23598, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_g(pbuffer, 254, 20, 86, 2719, 2749, 6576, 6642, 7771, 7936, 11826, 11856, 18025, 18205, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_g(pbuffer, 6810, 6576, 6642, 7906, 7981, 14829, 14895, 15489, 15654, 18175, 18250, 22453, 22633, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_g(pbuffer, 15063, 14829, 14895, 15624, 15699, 21039, 21105, 21363, 21510, 22603, 22678, 24678, 24825, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_g(pbuffer, 21213, 21039, 21105, 21480, 21555, 24111, 24147, 24252, 24303, 24795, 24870, 25705, 25756, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_s(pbuffer, 13696, 13306, 13429, 20079, 20202, 1, 11046, 11292, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_s(pbuffer, 17356, 16606, 16930, 21888, 22080, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_s(pbuffer, 18757, 18085, 18295, 22513, 22723, 1, 16606, 16930, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_s(pbuffer, 23149, 22513, 22723, 24723, 24915, 1, 21888, 22080, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_p(pbuffer, 10535, 9800, 10037, 11202, 11292, 16627, 16993, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_p(pbuffer, 12150, 11787, 11901, 13384, 13429, 18106, 18358, 1, 9800, 10037, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_p(pbuffer, 18883, 18106, 18358, 20157, 20202, 22534, 22786, 1, 16627, 16993, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_d(pbuffer, 4544, 4154, 4277, 5161, 5206, 9821, 10100, 11046, 11292, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_d(pbuffer, 8488, 7816, 8026, 9947, 10037, 15534, 15744, 16606, 16930, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_d(pbuffer, 16170, 15534, 15744, 16840, 16993, 21408, 21600, 21888, 22080, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_f(pbuffer, 3043, 2680, 2794, 4232, 4277, 7837, 8089, 9800, 10037, 13384, 13429, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_f(pbuffer, 8614, 7837, 8089, 9992, 10100, 15555, 15807, 16627, 16993, 20157, 20202, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_g(pbuffer, 344, 35, 131, 2749, 2794, 6591, 6687, 7816, 8026, 11856, 11901, 18085, 18295, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_g(pbuffer, 6900, 6591, 6687, 7981, 8089, 14844, 14940, 15534, 15744, 18250, 18358, 22513, 22723, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_g(pbuffer, 15153, 14844, 14940, 15699, 15807, 21054, 21150, 21408, 21600, 22678, 22786, 24723, 24915, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_s(pbuffer, 13822, 13336, 13510, 20109, 20265, 1, 11106, 11355, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_p(pbuffer, 12276, 11808, 11964, 13492, 13510, 18127, 18421, 1, 9851, 10163, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_d(pbuffer, 4670, 4184, 4358, 5269, 5287, 9869, 10199, 11106, 11355, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_f(pbuffer, 3169, 2701, 2857, 4340, 4358, 7858, 8152, 9851, 10163, 13492, 13510, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_s(pbuffer, 13882, 13354, 13546, 20127, 20301, 1, 11142, 11391, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_s(pbuffer, 19009, 18145, 18457, 22573, 22849, 1, 16705, 17056, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_p(pbuffer, 10661, 9887, 10235, 11355, 11391, 16735, 17116, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_p(pbuffer, 12336, 11826, 12000, 13510, 13546, 18175, 18517, 1, 9887, 10235, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_p(pbuffer, 19109, 18175, 18517, 20265, 20301, 22603, 22909, 1, 16735, 17116, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_d(pbuffer, 4730, 4202, 4394, 5287, 5323, 9917, 10295, 11142, 11391, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_d(pbuffer, 8740, 7876, 8188, 10163, 10235, 15594, 15870, 16705, 17056, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_f(pbuffer, 3229, 2719, 2893, 4358, 4394, 7906, 8248, 9887, 10235, 13510, 13546, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_f(pbuffer, 8840, 7906, 8248, 10199, 10295, 15624, 15930, 16735, 17116, 20265, 20301, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_g(pbuffer, 470, 56, 194, 2857, 2893, 6612, 6750, 7876, 8188, 11964, 12000, 18145, 18457, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_g(pbuffer, 7026, 6612, 6750, 8152, 8248, 14865, 15003, 15594, 15870, 18421, 18517, 22573, 22849, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_s(pbuffer, 13982, 13384, 13606, 20157, 20361, 1, 11202, 11451, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_s(pbuffer, 17482, 16795, 17176, 21990, 22143, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_s(pbuffer, 19209, 18205, 18577, 22633, 22969, 1, 16795, 17176, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_s(pbuffer, 23275, 22633, 22969, 24825, 24978, 1, 21990, 22143, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_p(pbuffer, 10761, 9947, 10355, 11391, 11451, 16840, 17266, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_p(pbuffer, 12436, 11856, 12060, 13546, 13606, 18250, 18667, 1, 9947, 10355, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_p(pbuffer, 19359, 18250, 18667, 20301, 20361, 22678, 23059, 1, 16840, 17266, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_d(pbuffer, 4830, 4232, 4454, 5323, 5383, 9992, 10445, 11202, 11451, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_d(pbuffer, 8940, 7936, 8308, 10235, 10355, 15654, 15990, 16795, 17176, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_d(pbuffer, 16296, 15654, 15990, 17116, 17266, 21510, 21663, 21990, 22143, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_f(pbuffer, 3329, 2749, 2953, 4394, 4454, 7981, 8398, 9947, 10355, 13546, 13606, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_f(pbuffer, 9090, 7981, 8398, 10295, 10445, 15699, 16080, 16840, 17266, 20301, 20361, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_g(pbuffer, 570, 86, 254, 2893, 2953, 6642, 6810, 7936, 8308, 12000, 12060, 18205, 18577, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_g(pbuffer, 7126, 6642, 6810, 8248, 8398, 14895, 15063, 15654, 15990, 18517, 18667, 22633, 22969, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_g(pbuffer, 15279, 14895, 15063, 15930, 16080, 21105, 21213, 21510, 21663, 22909, 23059, 24825, 24978, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_s(pbuffer, 19509, 18295, 18757, 22723, 23149, 1, 16930, 17356, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_p(pbuffer, 12586, 11901, 12150, 13606, 13696, 18358, 18883, 1, 10037, 10535, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_d(pbuffer, 9240, 8026, 8488, 10355, 10535, 15744, 16170, 16930, 17356, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_f(pbuffer, 3479, 2794, 3043, 4454, 4544, 8089, 8614, 10037, 10535, 13606, 13696, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_g(pbuffer, 720, 131, 344, 2953, 3043, 6687, 6900, 8026, 8488, 12060, 12150, 18295, 18757, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_g(pbuffer, 7276, 6687, 6900, 8398, 8614, 14940, 15153, 15744, 16170, 18667, 18883, 22723, 23149, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gf_p(pbuffer, 12796, 12000, 12336, 13822, 13882, 18517, 19109, 1, 10235, 10661, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gf_f(pbuffer, 3689, 2893, 3229, 4670, 4730, 8248, 8840, 10235, 10661, 13822, 13882, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gf_g(pbuffer, 930, 194, 470, 3169, 3229, 6750, 7026, 8188, 8740, 12276, 12336, 18457, 19009, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gg_s(pbuffer, 19719, 18577, 19209, 22969, 23275, 1, 17176, 17482, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gg_p(pbuffer, 12946, 12060, 12436, 13882, 13982, 18667, 19359, 1, 10355, 10761, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gg_d(pbuffer, 9450, 8308, 8940, 10661, 10761, 15990, 16296, 17176, 17482, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gg_f(pbuffer, 3839, 2953, 3329, 4730, 4830, 8398, 9090, 10355, 10761, 13882, 13982, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gg_g(pbuffer, 1080, 254, 570, 3229, 3329, 6810, 7126, 8308, 8940, 12336, 12436, 18577, 19209, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gg_g(pbuffer, 7486, 6810, 7126, 8840, 9090, 15063, 15279, 15990, 16296, 19109, 19359, 22969, 23275, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gh_g(pbuffer, 1980, 344, 720, 3329, 3479, 6900, 7276, 8488, 9240, 12436, 12586, 18757, 19509, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_hg_g(pbuffer, 2295, 570, 1080, 3689, 3839, 7126, 7486, 8940, 9450, 12796, 12946, 19209, 19719, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2cgeom::comp_prim_op_geom_010_gg(pbuffer, 1305, 570, 930, 1980, 2295, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, 0, pbuffer, 1305, 675, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2pecp namespace

#endif /* ProjectedCorePotentialGeom010GGForG_hpp */
