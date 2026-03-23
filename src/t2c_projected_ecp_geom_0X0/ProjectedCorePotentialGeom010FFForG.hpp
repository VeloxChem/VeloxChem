#ifndef ProjectedCorePotentialGeom010FFForG_hpp
#define ProjectedCorePotentialGeom010FFForG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "SimdArray.hpp"
#include "ProjectedCorePotentialPrimRecDDForD.hpp"
#include "ProjectedCorePotentialPrimRecDDForF.hpp"
#include "ProjectedCorePotentialPrimRecDDForG.hpp"
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
#include "ProjectedCorePotentialPrimRecDPForD.hpp"
#include "ProjectedCorePotentialPrimRecDPForF.hpp"
#include "ProjectedCorePotentialPrimRecDPForP.hpp"
#include "ProjectedCorePotentialPrimRecDPForS.hpp"
#include "ProjectedCorePotentialPrimRecFDForF.hpp"
#include "ProjectedCorePotentialPrimRecFDForG.hpp"
#include "ProjectedCorePotentialPrimRecFDForP.hpp"
#include "ProjectedCorePotentialPrimRecFFForD.hpp"
#include "ProjectedCorePotentialPrimRecFFForF.hpp"
#include "ProjectedCorePotentialPrimRecFFForG.hpp"
#include "ProjectedCorePotentialPrimRecFFForP.hpp"
#include "ProjectedCorePotentialPrimRecFFForS.hpp"
#include "ProjectedCorePotentialPrimRecFGForG.hpp"
#include "ProjectedCorePotentialPrimRecGFForG.hpp"
#include "ProjectedCorePotentialPrimRecPDForD.hpp"
#include "ProjectedCorePotentialPrimRecPDForF.hpp"
#include "ProjectedCorePotentialPrimRecPDForG.hpp"
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
#include "ProjectedCorePotentialPrimRecPPForD.hpp"
#include "ProjectedCorePotentialPrimRecPPForF.hpp"
#include "ProjectedCorePotentialPrimRecPPForP.hpp"
#include "ProjectedCorePotentialPrimRecPPForS.hpp"
#include "ProjectedCorePotentialPrimRecPSForD.hpp"
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
#include "ProjectedCorePotentialPrimRecSPForD.hpp"
#include "ProjectedCorePotentialPrimRecSPForF.hpp"
#include "ProjectedCorePotentialPrimRecSPForG.hpp"
#include "ProjectedCorePotentialPrimRecSPForP.hpp"
#include "ProjectedCorePotentialPrimRecSPForS.hpp"
#include "ProjectedCorePotentialPrimRecSS.hpp"
#include "GeometricalDerivatives010ForFF.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes (F|U_l|F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_projected_core_potential_geom_010_ff_for_g(T& distributor,
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

    CSimdArray<double> i_values(12, ket_npgtos);

    CSimdArray<double> l_values(5, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(8067, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(300, 1);

    CSimdArray<double> sbuffer(147, 1);

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

                    t2cfunc::comp_i_vals(i_values, 11, pfactors, 8);

                    t2cfunc::comp_l_vals(l_values, 4, pfactors, 8, 6);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 0, pbuffer, 0, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 1, pbuffer, 1074, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 1, pbuffer, 1075, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 2, pbuffer, 1595, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 2, pbuffer, 1596, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 2, pbuffer, 1597, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 3, pbuffer, 1879, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 3, pbuffer, 1880, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 3, pbuffer, 1881, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 3, pbuffer, 1882, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 4, pbuffer, 2008, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 4, pbuffer, 2009, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 4, pbuffer, 2010, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 4, pbuffer, 2011, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 5, pbuffer, 2043, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 5, pbuffer, 2044, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 0, pbuffer, 2048, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 1, pbuffer, 2068, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 1, pbuffer, 2069, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 1, pbuffer, 2070, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 2, pbuffer, 2118, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 2, pbuffer, 2119, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 2, pbuffer, 2120, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 2, pbuffer, 2121, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 3, pbuffer, 2172, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 3, pbuffer, 2173, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 3, pbuffer, 2174, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 3, pbuffer, 2175, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 3, pbuffer, 2176, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 4, pbuffer, 2227, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 4, pbuffer, 2228, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 4, pbuffer, 2229, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 5, pbuffer, 2242, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 0, pbuffer, 2243, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 1, pbuffer, 2253, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 1, pbuffer, 2254, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 1, pbuffer, 2255, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 2, pbuffer, 2277, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 2, pbuffer, 2278, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 2, pbuffer, 2279, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 2, pbuffer, 2280, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 2, pbuffer, 2281, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 3, pbuffer, 2312, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 3, pbuffer, 2313, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 3, pbuffer, 2314, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 3, pbuffer, 2315, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 4, pbuffer, 2337, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 4, pbuffer, 2338, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 3, 0, pbuffer, 2342, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 1, pbuffer, 2346, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 1, pbuffer, 2347, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 3, 1, pbuffer, 2348, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 2, pbuffer, 2355, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 2, pbuffer, 2356, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 2, pbuffer, 2357, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 2, pbuffer, 2358, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 3, pbuffer, 2368, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 3, pbuffer, 2369, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 3, pbuffer, 2370, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 4, pbuffer, 2377, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 4, 0, pbuffer, 2378, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 1, pbuffer, 2379, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 1, pbuffer, 2380, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 2, pbuffer, 2381, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 2, pbuffer, 2382, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 2, pbuffer, 2383, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 3, pbuffer, 2384, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 3, pbuffer, 2385, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 0, pbuffer, 2386, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 1, pbuffer, 2800, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 1, pbuffer, 2801, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 1, pbuffer, 2802, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 2, pbuffer, 3467, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 2, pbuffer, 3468, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 2, pbuffer, 3469, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 2, pbuffer, 3470, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 3, pbuffer, 3822, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 3, pbuffer, 3823, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 3, pbuffer, 3824, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 3, pbuffer, 3825, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 3, pbuffer, 3826, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 4, pbuffer, 3949, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 4, pbuffer, 3950, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 4, pbuffer, 3951, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 5, pbuffer, 3964, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 0, pbuffer, 3965, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 1, pbuffer, 3985, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 1, pbuffer, 3986, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 1, pbuffer, 3987, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 1, pbuffer, 3988, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 2, pbuffer, 4536, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 2, pbuffer, 4537, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 2, pbuffer, 4538, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 2, pbuffer, 4539, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 2, pbuffer, 4540, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 3, pbuffer, 4844, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 3, pbuffer, 4845, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 3, pbuffer, 4846, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 3, pbuffer, 4847, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 4, pbuffer, 4879, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 4, pbuffer, 4880, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 0, pbuffer, 4884, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 1, pbuffer, 4894, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 1, pbuffer, 4895, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 1, pbuffer, 4896, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 1, pbuffer, 4897, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 1, pbuffer, 4898, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 2, pbuffer, 4958, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 2, pbuffer, 4959, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 2, pbuffer, 4960, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 2, pbuffer, 4961, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 2, pbuffer, 4962, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 3, pbuffer, 5003, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 3, pbuffer, 5004, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 3, pbuffer, 5005, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 4, pbuffer, 5018, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 3, 0, pbuffer, 5019, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 1, pbuffer, 5023, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 1, pbuffer, 5024, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 1, pbuffer, 5025, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 1, pbuffer, 5026, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 3, 1, pbuffer, 5027, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 2, pbuffer, 5052, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 2, pbuffer, 5053, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 2, pbuffer, 5054, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 2, pbuffer, 5055, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 3, pbuffer, 5071, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 3, pbuffer, 5072, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 4, 0, pbuffer, 5076, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 1, pbuffer, 5077, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 1, pbuffer, 5078, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 1, pbuffer, 5079, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 4, 1, pbuffer, 5080, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 2, pbuffer, 5087, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 2, pbuffer, 5088, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 2, pbuffer, 5089, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 3, pbuffer, 5093, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 1, pbuffer, 5094, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 1, pbuffer, 5095, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 2, pbuffer, 5096, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 0, pbuffer, 5097, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 1, pbuffer, 5285, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 1, pbuffer, 5286, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 1, pbuffer, 5287, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 2, pbuffer, 5576, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 2, pbuffer, 5577, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 2, pbuffer, 5578, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 2, pbuffer, 5579, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 2, pbuffer, 5580, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 3, pbuffer, 5881, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 3, pbuffer, 5882, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 3, pbuffer, 5883, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 3, pbuffer, 5884, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 4, pbuffer, 5916, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 4, pbuffer, 5917, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 0, pbuffer, 5921, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 1, pbuffer, 5941, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 1, pbuffer, 5942, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 1, pbuffer, 5943, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 1, pbuffer, 5944, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 1, pbuffer, 5945, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 2, pbuffer, 6638, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 2, pbuffer, 6639, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 2, pbuffer, 6640, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 2, pbuffer, 6641, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 2, pbuffer, 6642, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 3, pbuffer, 6765, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 3, pbuffer, 6766, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 3, pbuffer, 6767, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 4, pbuffer, 6780, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 2, 0, pbuffer, 6781, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 1, pbuffer, 6791, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 1, pbuffer, 6792, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 1, pbuffer, 6793, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 1, pbuffer, 6794, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 2, 1, pbuffer, 6795, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 2, pbuffer, 6855, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 2, pbuffer, 6856, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 2, pbuffer, 6857, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 2, pbuffer, 6858, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 3, pbuffer, 6890, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 3, pbuffer, 6891, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 3, 0, pbuffer, 6895, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 1, pbuffer, 6899, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 1, pbuffer, 6900, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 1, pbuffer, 6901, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 1, pbuffer, 6902, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 3, 1, pbuffer, 6903, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 2, pbuffer, 6928, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 2, pbuffer, 6929, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 2, pbuffer, 6930, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 3, pbuffer, 6943, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 4, 0, pbuffer, 6944, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 1, pbuffer, 6945, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 1, pbuffer, 6946, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 4, 1, pbuffer, 6947, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 4, 1, pbuffer, 6948, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 2, pbuffer, 6955, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 2, pbuffer, 6956, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 1, pbuffer, 6960, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 5, 1, pbuffer, 6961, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 2, pbuffer, 6962, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 0, 0, pbuffer, 6963, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 1, pbuffer, 7028, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 1, pbuffer, 7029, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 0, 1, pbuffer, 7030, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 2, pbuffer, 7123, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 2, pbuffer, 7124, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 2, pbuffer, 7125, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 2, pbuffer, 7126, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 3, pbuffer, 7222, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 3, pbuffer, 7223, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 3, pbuffer, 7224, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 4, pbuffer, 7237, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 1, 0, pbuffer, 7238, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 1, pbuffer, 7258, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 1, pbuffer, 7259, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 1, pbuffer, 7260, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 1, pbuffer, 7261, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 1, 1, pbuffer, 7262, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 2, pbuffer, 7563, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 2, pbuffer, 7564, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 2, pbuffer, 7565, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 2, pbuffer, 7566, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 3, pbuffer, 7598, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 3, pbuffer, 7599, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 2, 0, pbuffer, 7603, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 1, pbuffer, 7613, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 1, pbuffer, 7614, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 1, pbuffer, 7615, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 2, 1, pbuffer, 7616, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 2, 1, pbuffer, 7617, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 2, pbuffer, 7668, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 2, pbuffer, 7669, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 2, pbuffer, 7670, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 3, pbuffer, 7683, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 3, 0, pbuffer, 7684, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 1, pbuffer, 7688, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 1, pbuffer, 7689, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 3, 1, pbuffer, 7690, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 3, 1, pbuffer, 7691, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 2, pbuffer, 7713, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 2, pbuffer, 7714, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 3, 4, 0, pbuffer, 7718, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 1, pbuffer, 7719, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 4, 1, pbuffer, 7720, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 4, 1, pbuffer, 7721, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 2, pbuffer, 7728, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 5, 1, pbuffer, 7729, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 5, 1, pbuffer, 7730, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 0, 0, pbuffer, 7731, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 0, 1, pbuffer, 7751, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 0, 1, pbuffer, 7752, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 2, pbuffer, 7781, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 0, 2, pbuffer, 7782, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 0, 2, pbuffer, 7783, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 3, pbuffer, 7815, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 0, 3, pbuffer, 7816, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 1, 0, pbuffer, 7820, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 1, pbuffer, 7830, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 1, 1, pbuffer, 7831, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 1, 1, pbuffer, 7832, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 1, 1, pbuffer, 7833, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 2, pbuffer, 7929, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 1, 2, pbuffer, 7930, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 1, 2, pbuffer, 7931, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 3, pbuffer, 7944, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 2, 0, pbuffer, 7945, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 1, pbuffer, 7949, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 2, 1, pbuffer, 7950, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 2, 1, pbuffer, 7951, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 2, 1, pbuffer, 7952, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 2, pbuffer, 7984, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 2, 2, pbuffer, 7985, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 4, 3, 0, pbuffer, 7989, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 1, pbuffer, 7990, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 3, 1, pbuffer, 7991, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 3, 1, pbuffer, 7992, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 2, pbuffer, 8005, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 4, 1, pbuffer, 8006, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 4, 1, pbuffer, 8007, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 5, 1, pbuffer, 8011, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 1, pbuffer, 8012, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 1, 1, pbuffer, 8013, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 2, pbuffer, 8042, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 2, 1, pbuffer, 8046, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 2, 1, pbuffer, 8047, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 2, 2, pbuffer, 8060, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 3, 1, pbuffer, 8061, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 5, 3, 1, pbuffer, 8062, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 4, 1, pbuffer, 8066, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2012, 2008, 2227, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2045, 2043, 2242, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2177, 2172, 2312, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2230, 2227, 2337, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2282, 2277, 2355, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2316, 2312, 2368, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2339, 2337, 2377, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2359, 2355, 2381, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2371, 2368, 2384, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3827, 3822, 4844, 1, 2008, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3952, 3949, 4879, 1, 2043, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4541, 4536, 4958, 1, 2172, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4848, 4844, 5003, 1, 2227, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4881, 4879, 5018, 1, 2242, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4899, 4894, 5023, 1, 2277, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4963, 4958, 5052, 1, 2312, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5006, 5003, 5071, 1, 2337, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5028, 5023, 5077, 1, 2355, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5056, 5052, 5087, 1, 2368, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5073, 5071, 5093, 1, 2377, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5081, 5077, 5094, 1, 2381, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5090, 5087, 5096, 1, 2384, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5581, 5576, 6638, 2, 3822, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5885, 5881, 6765, 2, 3949, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5918, 5916, 6780, 2, 3964, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5946, 5941, 6791, 2, 4536, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6643, 6638, 6855, 2, 4844, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6768, 6765, 6890, 2, 4879, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6796, 6791, 6899, 2, 4958, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6859, 6855, 6928, 2, 5003, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6892, 6890, 6943, 2, 5018, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6904, 6899, 6945, 2, 5052, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6931, 6928, 6955, 2, 5071, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6949, 6945, 6960, 2, 5087, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6957, 6955, 6962, 2, 5093, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7127, 7123, 7563, 3, 5881, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7225, 7222, 7598, 3, 5916, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7263, 7258, 7613, 3, 6638, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7567, 7563, 7668, 3, 6765, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7600, 7598, 7683, 3, 6780, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7618, 7613, 7688, 3, 6855, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7671, 7668, 7713, 3, 6890, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7692, 7688, 7719, 3, 6928, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7715, 7713, 7728, 3, 6943, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7722, 7719, 7729, 3, 6955, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7784, 7781, 7929, 4, 7222, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7817, 7815, 7944, 4, 7237, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7834, 7830, 7949, 4, 7563, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7932, 7929, 7984, 4, 7598, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7953, 7949, 7990, 4, 7668, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7986, 7984, 8005, 4, 7683, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7993, 7990, 8006, 4, 7713, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 8008, 8006, 8011, 4, 7728, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 8014, 8012, 8046, 5, 7929, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 8043, 8042, 8060, 5, 7944, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 8048, 8046, 8061, 5, 7984, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 8063, 8061, 8066, 5, 8005, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1883, 1879, 2008, 2173, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2015, 2009, 2043, 2228, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2122, 2118, 2172, 2278, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2180, 2173, 2227, 2313, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2233, 2228, 2242, 2338, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2285, 2278, 2312, 2356, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2319, 2313, 2337, 2369, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2362, 2356, 2368, 2382, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2374, 2369, 2377, 2385, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3471, 3467, 3822, 4537, 1, 1879, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3830, 3823, 3949, 4845, 1, 2009, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3955, 3950, 3964, 4880, 1, 2044, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3989, 3985, 4536, 4895, 1, 2118, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4544, 4537, 4844, 4959, 1, 2173, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4851, 4845, 4879, 5004, 1, 2228, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4902, 4895, 4958, 5024, 1, 2278, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4966, 4959, 5003, 5053, 1, 2313, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5009, 5004, 5018, 5072, 1, 2338, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5031, 5024, 5052, 5078, 1, 2356, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5059, 5053, 5071, 5088, 1, 2369, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5084, 5078, 5087, 5095, 1, 2382, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5584, 5577, 5881, 6639, 2, 3823, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5888, 5882, 5916, 6766, 2, 3950, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5949, 5942, 6638, 6792, 2, 4537, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6646, 6639, 6765, 6856, 2, 4845, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6771, 6766, 6780, 6891, 2, 4880, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6799, 6792, 6855, 6900, 2, 4959, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6862, 6856, 6890, 6929, 2, 5004, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6907, 6900, 6928, 6946, 2, 5053, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6934, 6929, 6943, 6956, 2, 5072, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6952, 6946, 6955, 6961, 2, 5088, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7130, 7124, 7222, 7564, 3, 5882, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7228, 7223, 7237, 7599, 3, 5917, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7266, 7259, 7563, 7614, 3, 6639, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7570, 7564, 7598, 7669, 3, 6766, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7621, 7614, 7668, 7689, 3, 6856, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7674, 7669, 7683, 7714, 3, 6891, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7695, 7689, 7713, 7720, 3, 6929, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7725, 7720, 7728, 7730, 3, 6956, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7787, 7782, 7815, 7930, 4, 7223, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7837, 7831, 7929, 7950, 4, 7564, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7935, 7930, 7944, 7985, 4, 7599, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7956, 7950, 7984, 7991, 4, 7669, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7996, 7991, 8005, 8007, 4, 7714, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 8017, 8013, 8042, 8047, 5, 7930, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 8051, 8047, 8060, 8062, 5, 7985, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1598, 1595, 1879, 2119, 2172, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1886, 1880, 2009, 2174, 2227, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2018, 2010, 2044, 2229, 2242, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2071, 2068, 2118, 2253, 2277, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2125, 2119, 2173, 2279, 2312, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2183, 2174, 2228, 2314, 2337, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2256, 2253, 2278, 2346, 2355, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2288, 2279, 2313, 2357, 2368, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2322, 2314, 2338, 2370, 2377, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2349, 2346, 2356, 2379, 2381, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2365, 2357, 2369, 2383, 2384, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2803, 2800, 3467, 3986, 4536, 1, 1595, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3474, 3468, 3823, 4538, 4844, 1, 1880, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3833, 3824, 3950, 4846, 4879, 1, 2010, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3992, 3986, 4537, 4896, 4958, 1, 2119, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4547, 4538, 4845, 4960, 5003, 1, 2174, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4854, 4846, 4880, 5005, 5018, 1, 2229, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4905, 4896, 4959, 5025, 5052, 1, 2279, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4969, 4960, 5004, 5054, 5071, 1, 2314, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5034, 5025, 5053, 5079, 5087, 1, 2357, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5062, 5054, 5072, 5089, 5093, 1, 2370, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5288, 5285, 5577, 5943, 6638, 2, 3468, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5587, 5578, 5882, 6640, 6765, 2, 3824, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5891, 5883, 5917, 6767, 6780, 2, 3951, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5952, 5943, 6639, 6793, 6855, 2, 4538, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6649, 6640, 6766, 6857, 6890, 2, 4846, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6802, 6793, 6856, 6901, 6928, 2, 4960, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6865, 6857, 6891, 6930, 6943, 2, 5005, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 6910, 6901, 6929, 6947, 6955, 2, 5054, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7031, 7028, 7124, 7260, 7563, 3, 5578, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7133, 7125, 7223, 7565, 7598, 3, 5883, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7269, 7260, 7564, 7615, 7668, 3, 6640, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7573, 7565, 7599, 7670, 7683, 3, 6767, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7624, 7615, 7669, 7690, 7713, 3, 6857, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7698, 7690, 7714, 7721, 7728, 3, 6930, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7753, 7751, 7782, 7832, 7929, 4, 7125, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7790, 7783, 7816, 7931, 7944, 4, 7224, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7840, 7832, 7930, 7951, 7984, 4, 7565, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7959, 7951, 7985, 7992, 8005, 4, 7670, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1076, 1074, 1595, 2069, 2118, 4536, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1601, 1596, 1880, 2120, 2173, 4844, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1889, 1881, 2010, 2175, 2228, 4879, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2074, 2069, 2119, 2254, 2278, 4958, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2128, 2120, 2174, 2280, 2313, 5003, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2186, 2175, 2229, 2315, 2338, 5018, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2259, 2254, 2279, 2347, 2356, 5052, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2291, 2280, 2314, 2358, 2369, 5071, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2352, 2347, 2357, 2380, 2382, 5087, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2806, 2801, 3468, 3987, 4537, 6638, 1, 1596, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3477, 3469, 3824, 4539, 4845, 6765, 1, 1881, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3836, 3825, 3951, 4847, 4880, 6780, 1, 2011, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3995, 3987, 4538, 4897, 4959, 6855, 1, 2120, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4550, 4539, 4846, 4961, 5004, 6890, 1, 2175, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4908, 4897, 4960, 5026, 5053, 6928, 1, 2280, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4972, 4961, 5005, 5055, 5072, 6943, 1, 2315, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5037, 5026, 5054, 5080, 5088, 6955, 1, 2358, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5291, 5286, 5578, 5944, 6639, 7563, 2, 3469, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5590, 5579, 5883, 6641, 6766, 7598, 2, 3825, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5955, 5944, 6640, 6794, 6856, 7668, 2, 4539, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6652, 6641, 6767, 6858, 6891, 7683, 2, 4847, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6805, 6794, 6857, 6902, 6929, 7713, 2, 4961, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 6913, 6902, 6930, 6948, 6956, 7728, 2, 5055, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7034, 7029, 7125, 7261, 7564, 7929, 3, 5579, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7136, 7126, 7224, 7566, 7599, 7944, 3, 5884, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7272, 7261, 7565, 7616, 7669, 7984, 3, 6641, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7627, 7616, 7670, 7691, 7714, 8005, 3, 6858, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7756, 7752, 7783, 7833, 7930, 8042, 4, 7126, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7843, 7833, 7931, 7952, 7985, 8060, 4, 7566, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1, 0, 1074, 2048, 2068, 3985, 4894, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1079, 1075, 1596, 2070, 2119, 4537, 4958, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1604, 1597, 1881, 2121, 2174, 4845, 5003, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1892, 1882, 2011, 2176, 2229, 4880, 5018, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2049, 2048, 2069, 2243, 2253, 4895, 5023, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2077, 2070, 2120, 2255, 2279, 4959, 5052, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2131, 2121, 2175, 2281, 2314, 5004, 5071, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2244, 2243, 2254, 2342, 2346, 5024, 5077, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2262, 2255, 2280, 2348, 2357, 5053, 5087, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2343, 2342, 2347, 2378, 2379, 5078, 5094, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2387, 2386, 2801, 3965, 3986, 5942, 6791, 1, 1075, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 2809, 2802, 3469, 3988, 4538, 6639, 6855, 1, 1597, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 3480, 3470, 3825, 4540, 4846, 6766, 6890, 1, 1882, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 3966, 3965, 3987, 4884, 4896, 6792, 6899, 1, 2070, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 3998, 3988, 4539, 4898, 4960, 6856, 6928, 1, 2121, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 4553, 4540, 4847, 4962, 5005, 6891, 6943, 1, 2176, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 4885, 4884, 4897, 5019, 5025, 6900, 6945, 1, 2255, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 4911, 4898, 4961, 5027, 5054, 6929, 6955, 1, 2281, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5020, 5019, 5026, 5076, 5079, 6946, 6960, 1, 2348, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5098, 5097, 5286, 5921, 5943, 7259, 7613, 2, 2802, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5294, 5287, 5579, 5945, 6640, 7564, 7668, 2, 3470, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5593, 5580, 5884, 6642, 6767, 7599, 7683, 2, 3826, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5922, 5921, 5944, 6781, 6793, 7614, 7688, 2, 3988, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 5958, 5945, 6641, 6795, 6857, 7669, 7713, 2, 4540, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6782, 6781, 6794, 6895, 6901, 7689, 7719, 2, 4898, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6808, 6795, 6858, 6903, 6930, 7714, 7728, 2, 4962, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6896, 6895, 6902, 6944, 6947, 7720, 7729, 2, 5027, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 6964, 6963, 7029, 7238, 7260, 7831, 7949, 3, 5287, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7037, 7030, 7126, 7262, 7565, 7930, 7984, 3, 5580, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7239, 7238, 7261, 7603, 7615, 7950, 7990, 3, 5945, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7275, 7262, 7566, 7617, 7670, 7985, 8005, 3, 6642, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7604, 7603, 7616, 7684, 7690, 7991, 8006, 3, 6795, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7685, 7684, 7691, 7718, 7721, 8007, 8011, 3, 6903, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7732, 7731, 7752, 7820, 7832, 8013, 8046, 4, 7030, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7821, 7820, 7833, 7945, 7951, 8047, 8061, 4, 7262, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 7946, 7945, 7952, 7989, 7992, 8062, 8066, 4, 7617, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2021, 2008, 2012, 2227, 2230, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2189, 2172, 2177, 2312, 2316, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2236, 2227, 2230, 2337, 2339, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2294, 2277, 2282, 2355, 2359, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2325, 2312, 2316, 2368, 2371, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3839, 3822, 3827, 4844, 4848, 1, 2008, 2012, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3958, 3949, 3952, 4879, 4881, 1, 2043, 2045, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4556, 4536, 4541, 4958, 4963, 1, 2172, 2177, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4857, 4844, 4848, 5003, 5006, 1, 2227, 2230, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4914, 4894, 4899, 5023, 5028, 1, 2277, 2282, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4975, 4958, 4963, 5052, 5056, 1, 2312, 2316, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5012, 5003, 5006, 5071, 5073, 1, 2337, 2339, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5040, 5023, 5028, 5077, 5081, 1, 2355, 2359, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5065, 5052, 5056, 5087, 5090, 1, 2368, 2371, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5596, 5576, 5581, 6638, 6643, 2, 3822, 3827, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5894, 5881, 5885, 6765, 6768, 2, 3949, 3952, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5961, 5941, 5946, 6791, 6796, 2, 4536, 4541, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6655, 6638, 6643, 6855, 6859, 2, 4844, 4848, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6774, 6765, 6768, 6890, 6892, 2, 4879, 4881, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6811, 6791, 6796, 6899, 6904, 2, 4958, 4963, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6868, 6855, 6859, 6928, 6931, 2, 5003, 5006, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6916, 6899, 6904, 6945, 6949, 2, 5052, 5056, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 6937, 6928, 6931, 6955, 6957, 2, 5071, 5073, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7139, 7123, 7127, 7563, 7567, 3, 5881, 5885, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7231, 7222, 7225, 7598, 7600, 3, 5916, 5918, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7278, 7258, 7263, 7613, 7618, 3, 6638, 6643, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7576, 7563, 7567, 7668, 7671, 3, 6765, 6768, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7630, 7613, 7618, 7688, 7692, 3, 6855, 6859, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7677, 7668, 7671, 7713, 7715, 3, 6890, 6892, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7701, 7688, 7692, 7719, 7722, 3, 6928, 6931, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7793, 7781, 7784, 7929, 7932, 4, 7222, 7225, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7846, 7830, 7834, 7949, 7953, 4, 7563, 7567, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7938, 7929, 7932, 7984, 7986, 4, 7598, 7600, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7962, 7949, 7953, 7990, 7993, 4, 7668, 7671, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7999, 7990, 7993, 8006, 8008, 4, 7713, 7715, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 8020, 8012, 8014, 8046, 8048, 5, 7929, 7932, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 8054, 8046, 8048, 8061, 8063, 5, 7984, 7986, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1895, 1879, 1883, 2012, 2173, 2180, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2027, 2009, 2015, 2045, 2228, 2233, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2134, 2118, 2122, 2177, 2278, 2285, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2195, 2173, 2180, 2230, 2313, 2319, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2300, 2278, 2285, 2316, 2356, 2362, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2331, 2313, 2319, 2339, 2369, 2374, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3483, 3467, 3471, 3827, 4537, 4544, 1, 1879, 1883, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3845, 3823, 3830, 3952, 4845, 4851, 1, 2009, 2015, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4001, 3985, 3989, 4541, 4895, 4902, 1, 2118, 2122, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4562, 4537, 4544, 4848, 4959, 4966, 1, 2173, 2180, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4863, 4845, 4851, 4881, 5004, 5009, 1, 2228, 2233, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4920, 4895, 4902, 4963, 5024, 5031, 1, 2278, 2285, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4981, 4959, 4966, 5006, 5053, 5059, 1, 2313, 2319, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5046, 5024, 5031, 5056, 5078, 5084, 1, 2356, 2362, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5602, 5577, 5584, 5885, 6639, 6646, 2, 3823, 3830, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5900, 5882, 5888, 5918, 6766, 6771, 2, 3950, 3955, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5967, 5942, 5949, 6643, 6792, 6799, 2, 4537, 4544, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6661, 6639, 6646, 6768, 6856, 6862, 2, 4845, 4851, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6817, 6792, 6799, 6859, 6900, 6907, 2, 4959, 4966, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6874, 6856, 6862, 6892, 6929, 6934, 2, 5004, 5009, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 6922, 6900, 6907, 6931, 6946, 6952, 2, 5053, 5059, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7145, 7124, 7130, 7225, 7564, 7570, 3, 5882, 5888, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7284, 7259, 7266, 7567, 7614, 7621, 3, 6639, 6646, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7582, 7564, 7570, 7600, 7669, 7674, 3, 6766, 6771, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7636, 7614, 7621, 7671, 7689, 7695, 3, 6856, 6862, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7707, 7689, 7695, 7715, 7720, 7725, 3, 6929, 6934, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7799, 7782, 7787, 7817, 7930, 7935, 4, 7223, 7228, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7852, 7831, 7837, 7932, 7950, 7956, 4, 7564, 7570, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7968, 7950, 7956, 7986, 7991, 7996, 4, 7669, 7674, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 8026, 8013, 8017, 8043, 8047, 8051, 5, 7930, 7935, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1607, 1595, 1598, 1883, 2119, 2125, 2172, 2177, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1901, 1880, 1886, 2015, 2174, 2183, 2227, 2230, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2080, 2068, 2071, 2122, 2253, 2256, 2277, 2282, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2140, 2119, 2125, 2180, 2279, 2288, 2312, 2316, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2201, 2174, 2183, 2233, 2314, 2322, 2337, 2339, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2265, 2253, 2256, 2285, 2346, 2349, 2355, 2359, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2306, 2279, 2288, 2319, 2357, 2365, 2368, 2371, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2812, 2800, 2803, 3471, 3986, 3992, 4536, 4541, 1, 1595, 1598, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3489, 3468, 3474, 3830, 4538, 4547, 4844, 4848, 1, 1880, 1886, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3851, 3824, 3833, 3955, 4846, 4854, 4879, 4881, 1, 2010, 2018, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4007, 3986, 3992, 4544, 4896, 4905, 4958, 4963, 1, 2119, 2125, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4568, 4538, 4547, 4851, 4960, 4969, 5003, 5006, 1, 2174, 2183, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4926, 4896, 4905, 4966, 5025, 5034, 5052, 5056, 1, 2279, 2288, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4987, 4960, 4969, 5009, 5054, 5062, 5071, 5073, 1, 2314, 2322, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5297, 5285, 5288, 5584, 5943, 5952, 6638, 6643, 2, 3468, 3474, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5608, 5578, 5587, 5888, 6640, 6649, 6765, 6768, 2, 3824, 3833, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5973, 5943, 5952, 6646, 6793, 6802, 6855, 6859, 2, 4538, 4547, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6667, 6640, 6649, 6771, 6857, 6865, 6890, 6892, 2, 4846, 4854, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 6823, 6793, 6802, 6862, 6901, 6910, 6928, 6931, 2, 4960, 4969, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7040, 7028, 7031, 7130, 7260, 7269, 7563, 7567, 3, 5578, 5587, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7151, 7125, 7133, 7228, 7565, 7573, 7598, 7600, 3, 5883, 5891, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7290, 7260, 7269, 7570, 7615, 7624, 7668, 7671, 3, 6640, 6649, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7642, 7615, 7624, 7674, 7690, 7698, 7713, 7715, 3, 6857, 6865, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7759, 7751, 7753, 7787, 7832, 7840, 7929, 7932, 4, 7125, 7133, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7858, 7832, 7840, 7935, 7951, 7959, 7984, 7986, 4, 7565, 7573, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1082, 1074, 1076, 1598, 2069, 2074, 2118, 2122, 4541, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1613, 1596, 1601, 1886, 2120, 2128, 2173, 2180, 4848, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1907, 1881, 1889, 2018, 2175, 2186, 2228, 2233, 4881, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2086, 2069, 2074, 2125, 2254, 2259, 2278, 2285, 4963, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2146, 2120, 2128, 2183, 2280, 2291, 2313, 2319, 5006, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2271, 2254, 2259, 2288, 2347, 2352, 2356, 2362, 5056, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2818, 2801, 2806, 3474, 3987, 3995, 4537, 4544, 6643, 1, 1596, 1601, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3495, 3469, 3477, 3833, 4539, 4550, 4845, 4851, 6768, 1, 1881, 1889, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4013, 3987, 3995, 4547, 4897, 4908, 4959, 4966, 6859, 1, 2120, 2128, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4574, 4539, 4550, 4854, 4961, 4972, 5004, 5009, 6892, 1, 2175, 2186, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4932, 4897, 4908, 4969, 5026, 5037, 5053, 5059, 6931, 1, 2280, 2291, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5303, 5286, 5291, 5587, 5944, 5955, 6639, 6646, 7567, 2, 3469, 3477, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5614, 5579, 5590, 5891, 6641, 6652, 6766, 6771, 7600, 2, 3825, 3836, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5979, 5944, 5955, 6649, 6794, 6805, 6856, 6862, 7671, 2, 4539, 4550, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 6829, 6794, 6805, 6865, 6902, 6913, 6929, 6934, 7715, 2, 4961, 4972, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 7046, 7029, 7034, 7133, 7261, 7272, 7564, 7570, 7932, 3, 5579, 5590, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 7296, 7261, 7272, 7573, 7616, 7627, 7669, 7674, 7986, 3, 6641, 6652, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 7765, 7752, 7756, 7790, 7833, 7843, 7930, 7935, 8043, 4, 7126, 7136, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 4, 0, 1, 1076, 2048, 2049, 2068, 2071, 3989, 4894, 4899, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 1088, 1075, 1079, 1601, 2070, 2077, 2119, 2125, 4544, 4958, 4963, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 1619, 1597, 1604, 1889, 2121, 2131, 2174, 2183, 4851, 5003, 5006, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 2052, 2048, 2049, 2074, 2243, 2244, 2253, 2256, 4902, 5023, 5028, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 2092, 2070, 2077, 2128, 2255, 2262, 2279, 2288, 4966, 5052, 5056, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 2247, 2243, 2244, 2259, 2342, 2343, 2346, 2349, 5031, 5077, 5081, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 2390, 2386, 2387, 2806, 3965, 3966, 3986, 3992, 5949, 6791, 6796, 1, 1075, 1079, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 2824, 2802, 2809, 3477, 3988, 3998, 4538, 4547, 6646, 6855, 6859, 1, 1597, 1604, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 3501, 3470, 3480, 3836, 4540, 4553, 4846, 4854, 6771, 6890, 6892, 1, 1882, 1892, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 3969, 3965, 3966, 3995, 4884, 4885, 4896, 4905, 6799, 6899, 6904, 1, 2070, 2077, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 4019, 3988, 3998, 4550, 4898, 4911, 4960, 4969, 6862, 6928, 6931, 1, 2121, 2131, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 4888, 4884, 4885, 4908, 5019, 5020, 5025, 5034, 6907, 6945, 6949, 1, 2255, 2262, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5101, 5097, 5098, 5291, 5921, 5922, 5943, 5952, 7266, 7613, 7618, 2, 2802, 2809, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5309, 5287, 5294, 5590, 5945, 5958, 6640, 6649, 7570, 7668, 7671, 2, 3470, 3480, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5925, 5921, 5922, 5955, 6781, 6782, 6793, 6802, 7621, 7688, 7692, 2, 3988, 3998, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 5985, 5945, 5958, 6652, 6795, 6808, 6857, 6865, 7674, 7713, 7715, 2, 4540, 4553, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 6785, 6781, 6782, 6805, 6895, 6896, 6901, 6910, 7695, 7719, 7722, 2, 4898, 4911, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 6967, 6963, 6964, 7034, 7238, 7239, 7260, 7269, 7837, 7949, 7953, 3, 5287, 5294, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 7052, 7030, 7037, 7136, 7262, 7275, 7565, 7573, 7935, 7984, 7986, 3, 5580, 5593, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 7242, 7238, 7239, 7272, 7603, 7604, 7615, 7624, 7956, 7990, 7993, 3, 5945, 5958, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 7607, 7603, 7604, 7627, 7684, 7685, 7690, 7698, 7996, 8006, 8008, 3, 6795, 6808, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 7735, 7731, 7732, 7756, 7820, 7821, 7832, 7840, 8017, 8046, 8048, 4, 7030, 7037, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 7824, 7820, 7821, 7843, 7945, 7946, 7951, 7959, 8051, 8061, 8063, 4, 7262, 7275, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2033, 2012, 2021, 2230, 2236, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2207, 2177, 2189, 2316, 2325, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3857, 3827, 3839, 4848, 4857, 1, 2012, 2021, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4580, 4541, 4556, 4963, 4975, 1, 2177, 2189, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4869, 4848, 4857, 5006, 5012, 1, 2230, 2236, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4938, 4899, 4914, 5028, 5040, 1, 2282, 2294, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4993, 4963, 4975, 5056, 5065, 1, 2316, 2325, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 5620, 5581, 5596, 6643, 6655, 2, 3827, 3839, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 5906, 5885, 5894, 6768, 6774, 2, 3952, 3958, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 5991, 5946, 5961, 6796, 6811, 2, 4541, 4556, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6673, 6643, 6655, 6859, 6868, 2, 4848, 4857, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6835, 6796, 6811, 6904, 6916, 2, 4963, 4975, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 6880, 6859, 6868, 6931, 6937, 2, 5006, 5012, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7157, 7127, 7139, 7567, 7576, 3, 5885, 5894, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7302, 7263, 7278, 7618, 7630, 3, 6643, 6655, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7588, 7567, 7576, 7671, 7677, 3, 6768, 6774, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7648, 7618, 7630, 7692, 7701, 3, 6859, 6868, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7805, 7784, 7793, 7932, 7938, 4, 7225, 7231, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7864, 7834, 7846, 7953, 7962, 4, 7567, 7576, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7974, 7953, 7962, 7993, 7999, 4, 7671, 7677, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 8032, 8014, 8020, 8048, 8054, 5, 7932, 7938, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1913, 1883, 1895, 2021, 2180, 2195, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 2152, 2122, 2134, 2189, 2285, 2300, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 2217, 2180, 2195, 2236, 2319, 2331, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 3507, 3471, 3483, 3839, 4544, 4562, 1, 1883, 1895, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 3867, 3830, 3845, 3958, 4851, 4863, 1, 2015, 2027, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4025, 3989, 4001, 4556, 4902, 4920, 1, 2122, 2134, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4590, 4544, 4562, 4857, 4966, 4981, 1, 2180, 2195, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4948, 4902, 4920, 4975, 5031, 5046, 1, 2285, 2300, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 5630, 5584, 5602, 5894, 6646, 6661, 2, 3830, 3845, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 6001, 5949, 5967, 6655, 6799, 6817, 2, 4544, 4562, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 6683, 6646, 6661, 6774, 6862, 6874, 2, 4851, 4863, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 6845, 6799, 6817, 6868, 6907, 6922, 2, 4966, 4981, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 7167, 7130, 7145, 7231, 7570, 7582, 3, 5888, 5900, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 7312, 7266, 7284, 7576, 7621, 7636, 3, 6646, 6661, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 7658, 7621, 7636, 7677, 7695, 7707, 3, 6862, 6874, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 7874, 7837, 7852, 7938, 7956, 7968, 4, 7570, 7582, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 1625, 1598, 1607, 1895, 2125, 2140, 2177, 2189, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 1923, 1886, 1901, 2027, 2183, 2201, 2230, 2236, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2098, 2071, 2080, 2134, 2256, 2265, 2282, 2294, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2162, 2125, 2140, 2195, 2288, 2306, 2316, 2325, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2830, 2803, 2812, 3483, 3992, 4007, 4541, 4556, 1, 1598, 1607, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 3517, 3474, 3489, 3845, 4547, 4568, 4848, 4857, 1, 1886, 1901, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4035, 3992, 4007, 4562, 4905, 4926, 4963, 4975, 1, 2125, 2140, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4600, 4547, 4568, 4863, 4969, 4987, 5006, 5012, 1, 2183, 2201, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 5315, 5288, 5297, 5602, 5952, 5973, 6643, 6655, 2, 3474, 3489, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 5640, 5587, 5608, 5900, 6649, 6667, 6768, 6774, 2, 3833, 3851, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 6011, 5952, 5973, 6661, 6802, 6823, 6859, 6868, 2, 4547, 4568, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 7058, 7031, 7040, 7145, 7269, 7290, 7567, 7576, 3, 5587, 5608, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 7322, 7269, 7290, 7582, 7624, 7642, 7671, 7677, 3, 6649, 6667, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 7771, 7753, 7759, 7799, 7840, 7858, 7932, 7938, 4, 7133, 7151, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1094, 1076, 1082, 1607, 2074, 2086, 2122, 2134, 4556, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1635, 1601, 1613, 1901, 2128, 2146, 2180, 2195, 4857, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 2108, 2074, 2086, 2140, 2259, 2271, 2285, 2300, 4975, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 2840, 2806, 2818, 3489, 3995, 4013, 4544, 4562, 6655, 1, 1601, 1613, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 3527, 3477, 3495, 3851, 4550, 4574, 4851, 4863, 6774, 1, 1889, 1907, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 4045, 3995, 4013, 4568, 4908, 4932, 4966, 4981, 6868, 1, 2128, 2146, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 5325, 5291, 5303, 5608, 5955, 5979, 6646, 6661, 7576, 2, 3477, 3495, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 6021, 5955, 5979, 6667, 6805, 6829, 6862, 6874, 7677, 2, 4550, 4574, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 7068, 7034, 7046, 7151, 7272, 7296, 7570, 7582, 7938, 3, 5590, 5614, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 10, 1, 4, 1082, 2049, 2052, 2071, 2080, 4001, 4899, 4914, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 1104, 1079, 1088, 1613, 2077, 2092, 2125, 2140, 4562, 4963, 4975, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 2058, 2049, 2052, 2086, 2244, 2247, 2256, 2265, 4920, 5028, 5040, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 2396, 2387, 2390, 2818, 3966, 3969, 3992, 4007, 5967, 6796, 6811, 1, 1079, 1088, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 2850, 2809, 2824, 3495, 3998, 4019, 4547, 4568, 6661, 6859, 6868, 1, 1604, 1619, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 3975, 3966, 3969, 4013, 4885, 4888, 4905, 4926, 6817, 6904, 6916, 1, 2077, 2092, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 5107, 5098, 5101, 5303, 5922, 5925, 5952, 5973, 7284, 7618, 7630, 2, 2809, 2824, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 5335, 5294, 5309, 5614, 5958, 5985, 6649, 6667, 7582, 7671, 7677, 2, 3480, 3501, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 5931, 5922, 5925, 5979, 6782, 6785, 6802, 6823, 7636, 7692, 7701, 2, 3998, 4019, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 6973, 6964, 6967, 7046, 7239, 7242, 7269, 7290, 7852, 7953, 7962, 3, 5294, 5309, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 7248, 7239, 7242, 7296, 7604, 7607, 7624, 7642, 7968, 7993, 7999, 3, 5958, 5985, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 7741, 7732, 7735, 7765, 7821, 7824, 7840, 7858, 8026, 8048, 8054, 4, 7037, 7052, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 3877, 3839, 3857, 4857, 4869, 1, 2021, 2033, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4610, 4556, 4580, 4975, 4993, 1, 2189, 2207, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 5650, 5596, 5620, 6655, 6673, 2, 3839, 3857, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 6031, 5961, 5991, 6811, 6835, 2, 4556, 4580, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 6693, 6655, 6673, 6868, 6880, 2, 4857, 4869, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 7177, 7139, 7157, 7576, 7588, 3, 5894, 5906, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 7332, 7278, 7302, 7630, 7648, 3, 6655, 6673, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 7884, 7846, 7864, 7962, 7974, 4, 7576, 7588, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 1933, 1895, 1913, 2033, 2195, 2217, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 3537, 3483, 3507, 3857, 4562, 4590, 1, 1895, 1913, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 4055, 4001, 4025, 4580, 4920, 4948, 1, 2134, 2152, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 5665, 5602, 5630, 5906, 6661, 6683, 2, 3845, 3867, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 6046, 5967, 6001, 6673, 6817, 6845, 2, 4562, 4590, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 7347, 7284, 7312, 7588, 7636, 7658, 3, 6661, 6683, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 1645, 1607, 1625, 1913, 2140, 2162, 2189, 2207, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 2860, 2812, 2830, 3507, 4007, 4035, 4556, 4580, 1, 1607, 1625, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 3552, 3489, 3517, 3867, 4568, 4600, 4857, 4869, 1, 1901, 1923, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 5345, 5297, 5315, 5630, 5973, 6011, 6655, 6673, 2, 3489, 3517, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 7078, 7040, 7058, 7167, 7290, 7322, 7576, 7588, 3, 5608, 5640, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 1114, 1082, 1094, 1625, 2086, 2108, 2134, 2152, 4580, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 2875, 2818, 2840, 3517, 4013, 4045, 4562, 4590, 6673, 1, 1613, 1635, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 5360, 5303, 5325, 5640, 5979, 6021, 6661, 6683, 7588, 2, 3495, 3527, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 20, 4, 10, 1094, 2052, 2058, 2080, 2098, 4025, 4914, 4938, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 2406, 2390, 2396, 2840, 3969, 3975, 4007, 4035, 6001, 6811, 6835, 1, 1088, 1104, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 5117, 5101, 5107, 5325, 5925, 5931, 5973, 6011, 7312, 7630, 7648, 2, 2824, 2850, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 6983, 6967, 6973, 7068, 7242, 7248, 7290, 7322, 7874, 7962, 7974, 3, 5309, 5335, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_s(pbuffer, 4625, 4536, 6638, 1, 3822, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_p(pbuffer, 1948, 1879, 2008, 3823, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_d(pbuffer, 1660, 1595, 1879, 3468, 3822, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 3892, 3827, 5885, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 4628, 4541, 6643, 1, 3827, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 6708, 6643, 7567, 1, 5885, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 1951, 1883, 2008, 2012, 3830, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 3567, 3471, 3822, 3827, 5584, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 4070, 3989, 4536, 4541, 5949, 1, 3471, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 6061, 5949, 6638, 6643, 7266, 1, 5584, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_d(pbuffer, 1663, 1598, 1879, 1883, 3474, 3827, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_d(pbuffer, 3576, 3474, 3823, 3830, 5587, 5885, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_f(pbuffer, 1129, 1076, 1595, 1598, 2806, 3471, 4536, 4541, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_f(pbuffer, 2890, 2806, 3468, 3474, 5291, 5584, 6638, 6643, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 3901, 3839, 5894, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 4637, 4556, 6655, 1, 3839, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 5680, 5596, 7139, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 6070, 5961, 7278, 1, 5596, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 6717, 6655, 7576, 1, 5894, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 7362, 7278, 7846, 1, 7139, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 1960, 1895, 2012, 2021, 3845, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 3585, 3483, 3827, 3839, 5602, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 4079, 4001, 4541, 4556, 5967, 1, 3483, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 5698, 5602, 5885, 5894, 7145, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 6088, 5967, 6643, 6655, 7284, 1, 5602, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 7380, 7284, 7567, 7576, 7852, 1, 7145, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 1672, 1607, 1883, 1895, 3489, 3839, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 2899, 2812, 3471, 3483, 5297, 5596, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 3603, 3489, 3830, 3845, 5608, 5894, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 5375, 5297, 5584, 5602, 7040, 7139, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 1138, 1082, 1598, 1607, 2818, 3483, 4541, 4556, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 2917, 2818, 3474, 3489, 5303, 5602, 6643, 6655, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 5393, 5303, 5587, 5608, 7046, 7145, 7567, 7576, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_g(pbuffer, 35, 4, 1076, 1082, 2390, 2812, 3989, 4001, 5961, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_g(pbuffer, 2421, 2390, 2806, 2818, 5101, 5297, 5949, 5967, 7278, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_g(pbuffer, 5132, 5101, 5291, 5303, 6967, 7040, 7266, 7284, 7846, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 3919, 3857, 5906, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 4655, 4580, 6673, 1, 3857, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 5716, 5620, 7157, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 6106, 5991, 7302, 1, 5620, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 6735, 6673, 7588, 1, 5906, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 7192, 7157, 7805, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 7398, 7302, 7864, 1, 7157, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 7899, 7864, 8032, 1, 7805, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 1978, 1913, 2021, 2033, 3867, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 3621, 3507, 3839, 3857, 5630, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 4097, 4025, 4556, 4580, 6001, 1, 3507, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 5746, 5630, 5894, 5906, 7167, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 6136, 6001, 6655, 6673, 7312, 1, 5630, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 7428, 7312, 7576, 7588, 7874, 1, 7167, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 1690, 1625, 1895, 1913, 3517, 3857, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 2935, 2830, 3483, 3507, 5315, 5620, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 3651, 3517, 3845, 3867, 5640, 5906, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 5411, 5315, 5602, 5630, 7058, 7157, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 7093, 7058, 7145, 7167, 7771, 7805, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 1156, 1094, 1607, 1625, 2840, 3507, 4556, 4580, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 2965, 2840, 3489, 3517, 5325, 5630, 6655, 6673, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 5441, 5325, 5608, 5640, 7068, 7167, 7576, 7588, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 53, 10, 1082, 1094, 2396, 2830, 4001, 4025, 5991, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 2439, 2396, 2818, 2840, 5107, 5315, 5967, 6001, 7302, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 5150, 5107, 5303, 5325, 6973, 7058, 7284, 7312, 7864, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 6998, 6973, 7046, 7068, 7741, 7771, 7852, 7874, 8032, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 4685, 4610, 6693, 1, 3877, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 5776, 5650, 7177, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 6166, 6031, 7332, 1, 5650, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 7458, 7332, 7884, 1, 7177, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 3681, 3537, 3857, 3877, 5665, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 4127, 4055, 4580, 4610, 6046, 1, 3537, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 6211, 6046, 6673, 6693, 7347, 1, 5665, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 1720, 1645, 1913, 1933, 3552, 3877, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 2995, 2860, 3507, 3537, 5345, 5650, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 5471, 5345, 5630, 5665, 7078, 7177, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 1186, 1114, 1625, 1645, 2875, 3537, 4580, 4610, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 3040, 2875, 3517, 3552, 5360, 5665, 6673, 6693, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 83, 20, 1094, 1114, 2406, 2860, 4025, 4055, 6031, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 2469, 2406, 2840, 2875, 5117, 5345, 6001, 6046, 7332, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 5180, 5117, 5325, 5360, 6983, 7078, 7312, 7347, 7884, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_s(pbuffer, 4730, 4541, 4628, 6643, 6708, 1, 3827, 3892, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_p(pbuffer, 4172, 3989, 4070, 4625, 4628, 5949, 6061, 1, 3471, 3567, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_d(pbuffer, 1765, 1598, 1663, 1948, 1951, 3474, 3576, 3827, 3892, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_f(pbuffer, 1231, 1076, 1129, 1660, 1663, 2806, 2890, 3471, 3567, 4625, 4628, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 4748, 4556, 4637, 6655, 6717, 1, 3839, 3901, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 6256, 5961, 6070, 7278, 7362, 1, 5596, 5680, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 3726, 3483, 3585, 3892, 3901, 5602, 5698, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 4190, 4001, 4079, 4628, 4637, 5967, 6088, 1, 3483, 3585, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 6292, 5967, 6088, 6708, 6717, 7284, 7380, 1, 5602, 5698, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_d(pbuffer, 1783, 1607, 1672, 1951, 1960, 3489, 3603, 3839, 3901, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_d(pbuffer, 3085, 2812, 2899, 3567, 3585, 5297, 5375, 5596, 5680, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_f(pbuffer, 1249, 1082, 1138, 1663, 1672, 2818, 2917, 3483, 3585, 4628, 4637, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_f(pbuffer, 3121, 2818, 2917, 3576, 3603, 5303, 5393, 5602, 5698, 6708, 6717, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_g(pbuffer, 128, 4, 35, 1129, 1138, 2390, 2421, 2812, 2899, 4070, 4079, 5961, 6070, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_g(pbuffer, 2514, 2390, 2421, 2890, 2917, 5101, 5132, 5297, 5375, 6061, 6088, 7278, 7362, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 4784, 4580, 4655, 6673, 6735, 1, 3857, 3919, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 5821, 5620, 5716, 7157, 7192, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 6328, 5991, 6106, 7302, 7398, 1, 5620, 5716, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 7503, 7302, 7398, 7864, 7899, 1, 7157, 7192, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 3762, 3507, 3621, 3901, 3919, 5630, 5746, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 4226, 4025, 4097, 4637, 4655, 6001, 6136, 1, 3507, 3621, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 6388, 6001, 6136, 6717, 6735, 7312, 7428, 1, 5630, 5746, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 1819, 1625, 1690, 1960, 1978, 3517, 3651, 3857, 3919, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 3157, 2830, 2935, 3585, 3621, 5315, 5411, 5620, 5716, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 5516, 5315, 5411, 5698, 5746, 7058, 7093, 7157, 7192, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 1285, 1094, 1156, 1672, 1690, 2840, 2965, 3507, 3621, 4637, 4655, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 3217, 2840, 2965, 3603, 3651, 5325, 5441, 5630, 5746, 6717, 6735, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_g(pbuffer, 164, 10, 53, 1138, 1156, 2396, 2439, 2830, 2935, 4079, 4097, 5991, 6106, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_g(pbuffer, 2550, 2396, 2439, 2917, 2965, 5107, 5150, 5315, 5411, 6088, 6136, 7302, 7398, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_g(pbuffer, 5225, 5107, 5150, 5393, 5441, 6973, 6998, 7058, 7093, 7380, 7428, 7864, 7899, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 6448, 6031, 6166, 7332, 7458, 1, 5650, 5776, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 4286, 4055, 4127, 4655, 4685, 6046, 6211, 1, 3537, 3681, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 3277, 2860, 2995, 3621, 3681, 5345, 5471, 5650, 5776, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_f(pbuffer, 1345, 1114, 1186, 1690, 1720, 2875, 3040, 3537, 3681, 4655, 4685, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_g(pbuffer, 224, 20, 83, 1156, 1186, 2406, 2469, 2860, 2995, 4097, 4127, 6031, 6166, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_g(pbuffer, 2610, 2406, 2469, 2965, 3040, 5117, 5180, 5345, 5471, 6136, 6211, 7332, 7458, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_p(pbuffer, 4376, 4079, 4190, 4730, 4748, 6088, 6292, 1, 3585, 3726, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_f(pbuffer, 1435, 1138, 1249, 1765, 1783, 2917, 3121, 3585, 3726, 4730, 4748, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_g(pbuffer, 314, 35, 128, 1231, 1249, 2421, 2514, 2899, 3085, 4172, 4190, 6070, 6256, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_s(pbuffer, 6538, 6106, 6328, 7398, 7503, 1, 5716, 5821, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_p(pbuffer, 4436, 4097, 4226, 4748, 4784, 6136, 6388, 1, 3621, 3762, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_d(pbuffer, 3367, 2935, 3157, 3726, 3762, 5411, 5516, 5716, 5821, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_f(pbuffer, 1495, 1156, 1285, 1783, 1819, 2965, 3217, 3621, 3762, 4748, 4784, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_g(pbuffer, 374, 53, 164, 1249, 1285, 2439, 2550, 2935, 3157, 4190, 4226, 6106, 6328, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_g(pbuffer, 2700, 2439, 2550, 3121, 3217, 5150, 5225, 5411, 5516, 6292, 6388, 7398, 7503, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_g(pbuffer, 774, 83, 224, 1285, 1345, 2469, 2610, 2995, 3277, 4226, 4286, 6166, 6448, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gf_g(pbuffer, 924, 164, 374, 1435, 1495, 2550, 2700, 3157, 3367, 4376, 4436, 6328, 6538, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2cgeom::comp_prim_op_geom_010_ff(pbuffer, 474, 164, 314, 774, 924, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, 0, pbuffer, 474, 300, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<3, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2pecp namespace

#endif /* ProjectedCorePotentialGeom010FFForG_hpp */
