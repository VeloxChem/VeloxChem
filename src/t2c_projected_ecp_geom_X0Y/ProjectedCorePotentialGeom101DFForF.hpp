#ifndef ProjectedCorePotentialGeom101DFForF_hpp
#define ProjectedCorePotentialGeom101DFForF_hpp

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
#include "ProjectedCorePotentialPrimRecDFForS.hpp"
#include "ProjectedCorePotentialPrimRecDGForD.hpp"
#include "ProjectedCorePotentialPrimRecDGForF.hpp"
#include "ProjectedCorePotentialPrimRecDGForP.hpp"
#include "ProjectedCorePotentialPrimRecDGForS.hpp"
#include "ProjectedCorePotentialPrimRecDPForD.hpp"
#include "ProjectedCorePotentialPrimRecDPForS.hpp"
#include "ProjectedCorePotentialPrimRecFDForF.hpp"
#include "ProjectedCorePotentialPrimRecFGForF.hpp"
#include "ProjectedCorePotentialPrimRecPDForD.hpp"
#include "ProjectedCorePotentialPrimRecPDForF.hpp"
#include "ProjectedCorePotentialPrimRecPDForP.hpp"
#include "ProjectedCorePotentialPrimRecPDForS.hpp"
#include "ProjectedCorePotentialPrimRecPFForD.hpp"
#include "ProjectedCorePotentialPrimRecPFForP.hpp"
#include "ProjectedCorePotentialPrimRecPFForS.hpp"
#include "ProjectedCorePotentialPrimRecPGForD.hpp"
#include "ProjectedCorePotentialPrimRecPGForF.hpp"
#include "ProjectedCorePotentialPrimRecPGForP.hpp"
#include "ProjectedCorePotentialPrimRecPGForS.hpp"
#include "ProjectedCorePotentialPrimRecPPForD.hpp"
#include "ProjectedCorePotentialPrimRecPPForP.hpp"
#include "ProjectedCorePotentialPrimRecPPForS.hpp"
#include "ProjectedCorePotentialPrimRecPSForP.hpp"
#include "ProjectedCorePotentialPrimRecSDForD.hpp"
#include "ProjectedCorePotentialPrimRecSDForF.hpp"
#include "ProjectedCorePotentialPrimRecSDForP.hpp"
#include "ProjectedCorePotentialPrimRecSDForS.hpp"
#include "ProjectedCorePotentialPrimRecSFForD.hpp"
#include "ProjectedCorePotentialPrimRecSFForF.hpp"
#include "ProjectedCorePotentialPrimRecSFForP.hpp"
#include "ProjectedCorePotentialPrimRecSFForS.hpp"
#include "ProjectedCorePotentialPrimRecSGForD.hpp"
#include "ProjectedCorePotentialPrimRecSGForF.hpp"
#include "ProjectedCorePotentialPrimRecSGForP.hpp"
#include "ProjectedCorePotentialPrimRecSGForS.hpp"
#include "ProjectedCorePotentialPrimRecSPForD.hpp"
#include "ProjectedCorePotentialPrimRecSPForF.hpp"
#include "ProjectedCorePotentialPrimRecSPForP.hpp"
#include "ProjectedCorePotentialPrimRecSPForS.hpp"
#include "ProjectedCorePotentialPrimRecSS.hpp"
#include "GeometricalDerivatives1X1ForDF.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes (d^(1)/dA^(1)D|U_l|d^(1)/dB^(1)F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_projected_core_potential_geom_101_df_for_f(T& distributor,
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

    CSimdArray<double> i_values(11, ket_npgtos);

    CSimdArray<double> l_values(4, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(4335, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(540, 1);

    CSimdArray<double> sbuffer(315, 1);

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

                    t2cfunc::comp_i_vals(i_values, 10, pfactors, 8);

                    t2cfunc::comp_l_vals(l_values, 3, pfactors, 8, 6);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 0, pbuffer, 0, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 1, pbuffer, 974, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 1, pbuffer, 975, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 2, pbuffer, 1335, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 2, pbuffer, 1336, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 2, pbuffer, 1337, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 3, pbuffer, 1505, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 3, pbuffer, 1506, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 3, pbuffer, 1507, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 3, pbuffer, 1508, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 4, pbuffer, 1574, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 4, pbuffer, 1575, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 4, pbuffer, 1576, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 5, pbuffer, 1589, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 0, pbuffer, 1590, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 1, pbuffer, 1610, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 1, pbuffer, 1611, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 1, pbuffer, 1612, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 2, pbuffer, 1660, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 2, pbuffer, 1661, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 2, pbuffer, 1662, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 2, pbuffer, 1663, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 3, pbuffer, 1714, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 3, pbuffer, 1715, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 3, pbuffer, 1716, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 3, pbuffer, 1717, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 4, pbuffer, 1749, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 4, pbuffer, 1750, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 0, pbuffer, 1754, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 1, pbuffer, 1764, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 1, pbuffer, 1765, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 1, pbuffer, 1766, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 2, pbuffer, 1788, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 2, pbuffer, 1789, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 2, pbuffer, 1790, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 2, pbuffer, 1791, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 3, pbuffer, 1813, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 3, pbuffer, 1814, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 3, pbuffer, 1815, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 4, pbuffer, 1828, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 0, pbuffer, 1829, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 1, pbuffer, 1833, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 1, pbuffer, 1834, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 1, pbuffer, 1835, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 2, pbuffer, 1842, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 2, pbuffer, 1843, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 2, pbuffer, 1844, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 3, pbuffer, 1851, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 3, pbuffer, 1852, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 0, pbuffer, 1856, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 1, pbuffer, 1857, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 1, pbuffer, 1858, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 2, pbuffer, 1859, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 2, pbuffer, 1860, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 3, pbuffer, 1861, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 0, pbuffer, 1862, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 1, pbuffer, 2086, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 1, pbuffer, 2087, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 1, pbuffer, 2088, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 2, pbuffer, 2467, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 2, pbuffer, 2468, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 2, pbuffer, 2469, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 2, pbuffer, 2470, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 3, pbuffer, 2669, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 3, pbuffer, 2670, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 3, pbuffer, 2671, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 4, pbuffer, 2703, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 4, pbuffer, 2704, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 0, pbuffer, 2708, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 1, pbuffer, 2728, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 1, pbuffer, 2729, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 1, pbuffer, 2730, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 1, pbuffer, 2731, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 2, pbuffer, 3119, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 2, pbuffer, 3120, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 2, pbuffer, 3121, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 2, pbuffer, 3122, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 3, pbuffer, 3173, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 3, pbuffer, 3174, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 3, pbuffer, 3175, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 4, pbuffer, 3188, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 0, pbuffer, 3189, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 1, pbuffer, 3199, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 1, pbuffer, 3200, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 1, pbuffer, 3201, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 1, pbuffer, 3202, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 2, pbuffer, 3243, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 2, pbuffer, 3244, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 2, pbuffer, 3245, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 2, pbuffer, 3246, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 3, pbuffer, 3268, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 3, pbuffer, 3269, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 0, pbuffer, 3273, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 1, pbuffer, 3277, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 1, pbuffer, 3278, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 1, pbuffer, 3279, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 1, pbuffer, 3280, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 2, pbuffer, 3296, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 2, pbuffer, 3297, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 2, pbuffer, 3298, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 3, pbuffer, 3305, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 4, 0, pbuffer, 3306, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 1, pbuffer, 3307, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 1, pbuffer, 3308, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 1, pbuffer, 3309, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 2, pbuffer, 3313, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 2, pbuffer, 3314, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 1, pbuffer, 3315, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 0, pbuffer, 3316, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 1, pbuffer, 3414, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 1, pbuffer, 3415, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 1, pbuffer, 3416, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 2, pbuffer, 3567, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 2, pbuffer, 3568, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 2, pbuffer, 3569, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 3, pbuffer, 3632, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 3, pbuffer, 3633, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 4, pbuffer, 3646, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 0, pbuffer, 3647, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 1, pbuffer, 3667, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 1, pbuffer, 3668, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 1, pbuffer, 3669, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 1, pbuffer, 3670, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 2, pbuffer, 3854, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 2, pbuffer, 3855, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 2, pbuffer, 3856, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 3, pbuffer, 3888, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 3, pbuffer, 3889, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 0, pbuffer, 3893, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 1, pbuffer, 3903, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 1, pbuffer, 3904, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 1, pbuffer, 3905, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 1, pbuffer, 3906, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 2, pbuffer, 3947, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 2, pbuffer, 3948, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 2, pbuffer, 3949, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 3, pbuffer, 3962, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 0, pbuffer, 3963, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 1, pbuffer, 3967, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 1, pbuffer, 3968, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 1, pbuffer, 3969, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 1, pbuffer, 3970, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 2, pbuffer, 3986, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 2, pbuffer, 3987, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 4, 0, pbuffer, 3991, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 1, pbuffer, 3992, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 1, pbuffer, 3993, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 4, 1, pbuffer, 3994, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 2, pbuffer, 3998, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 1, pbuffer, 3999, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 0, pbuffer, 4000, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 1, pbuffer, 4035, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 1, pbuffer, 4036, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 2, pbuffer, 4090, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 2, pbuffer, 4091, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 3, pbuffer, 4120, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 0, pbuffer, 4124, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 1, pbuffer, 4144, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 1, pbuffer, 4145, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 1, pbuffer, 4146, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 2, pbuffer, 4209, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 2, pbuffer, 4210, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 3, pbuffer, 4223, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 2, 0, pbuffer, 4224, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 1, pbuffer, 4234, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 1, pbuffer, 4235, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 1, pbuffer, 4236, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 2, pbuffer, 4268, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 2, pbuffer, 4269, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 3, 0, pbuffer, 4273, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 1, pbuffer, 4277, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 1, pbuffer, 4278, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 3, 1, pbuffer, 4279, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 2, pbuffer, 4292, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 4, 0, pbuffer, 4293, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 1, pbuffer, 4294, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 4, 1, pbuffer, 4295, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 5, 1, pbuffer, 4299, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 1, pbuffer, 4300, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 1, pbuffer, 4320, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 1, pbuffer, 4330, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 4, 1, pbuffer, 4334, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1509, 1505, 1714, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1577, 1574, 1749, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1664, 1660, 1788, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1718, 1714, 1813, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1751, 1749, 1828, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1792, 1788, 1842, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1816, 1813, 1851, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1845, 1842, 1859, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1853, 1851, 1861, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2471, 2467, 3119, 1, 1505, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2672, 2669, 3173, 1, 1574, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2705, 2703, 3188, 1, 1589, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2732, 2728, 3199, 1, 1660, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3123, 3119, 3243, 1, 1714, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3176, 3173, 3268, 1, 1749, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3203, 3199, 3277, 1, 1788, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3247, 3243, 3296, 1, 1813, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3270, 3268, 3305, 1, 1828, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3281, 3277, 3307, 1, 1842, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3299, 3296, 3313, 1, 1851, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3310, 3307, 3315, 1, 1859, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3570, 3567, 3854, 2, 2669, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3634, 3632, 3888, 2, 2703, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3671, 3667, 3903, 2, 3119, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3857, 3854, 3947, 2, 3173, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3890, 3888, 3962, 2, 3188, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3907, 3903, 3967, 2, 3243, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3950, 3947, 3986, 2, 3268, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3971, 3967, 3992, 2, 3296, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3988, 3986, 3998, 2, 3305, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3995, 3992, 3999, 2, 3313, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4092, 4090, 4209, 3, 3632, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4121, 4120, 4223, 3, 3646, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4147, 4144, 4234, 3, 3854, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4211, 4209, 4268, 3, 3888, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4237, 4234, 4277, 3, 3947, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4270, 4268, 4292, 3, 3962, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4280, 4277, 4294, 3, 3986, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4296, 4294, 4299, 3, 3998, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4301, 4300, 4320, 4, 4209, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4321, 4320, 4330, 4, 4268, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4331, 4330, 4334, 4, 4292, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1338, 1335, 1505, 1661, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1512, 1506, 1574, 1715, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1580, 1575, 1589, 1750, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1613, 1610, 1660, 1764, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1667, 1661, 1714, 1789, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1721, 1715, 1749, 1814, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1767, 1764, 1788, 1833, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1795, 1789, 1813, 1843, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1819, 1814, 1828, 1852, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1836, 1833, 1842, 1857, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1848, 1843, 1851, 1860, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2089, 2086, 2467, 2729, 1, 1335, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2474, 2468, 2669, 3120, 1, 1506, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2675, 2670, 2703, 3174, 1, 1575, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2735, 2729, 3119, 3200, 1, 1661, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3126, 3120, 3173, 3244, 1, 1715, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3179, 3174, 3188, 3269, 1, 1750, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3206, 3200, 3243, 3278, 1, 1789, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3250, 3244, 3268, 3297, 1, 1814, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3284, 3278, 3296, 3308, 1, 1843, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3302, 3297, 3305, 3314, 1, 1852, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3417, 3414, 3567, 3668, 2, 2468, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3573, 3568, 3632, 3855, 2, 2670, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3637, 3633, 3646, 3889, 2, 2704, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3674, 3668, 3854, 3904, 2, 3120, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3860, 3855, 3888, 3948, 2, 3174, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3910, 3904, 3947, 3968, 2, 3244, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3953, 3948, 3962, 3987, 2, 3269, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3974, 3968, 3986, 3993, 2, 3297, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4037, 4035, 4090, 4145, 3, 3568, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4095, 4091, 4120, 4210, 3, 3633, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4150, 4145, 4209, 4235, 3, 3855, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4214, 4210, 4223, 4269, 3, 3889, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4240, 4235, 4268, 4278, 3, 3948, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4283, 4278, 4292, 4295, 3, 3987, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 976, 974, 1335, 1611, 1660, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1341, 1336, 1506, 1662, 1714, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1515, 1507, 1575, 1716, 1749, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1616, 1611, 1661, 1765, 1788, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1670, 1662, 1715, 1790, 1813, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1724, 1716, 1750, 1815, 1828, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1770, 1765, 1789, 1834, 1842, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1798, 1790, 1814, 1844, 1851, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1839, 1834, 1843, 1858, 1859, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2092, 2087, 2468, 2730, 3119, 1, 1336, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2477, 2469, 2670, 3121, 3173, 1, 1507, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2678, 2671, 2704, 3175, 3188, 1, 1576, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2738, 2730, 3120, 3201, 3243, 1, 1662, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3129, 3121, 3174, 3245, 3268, 1, 1716, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3209, 3201, 3244, 3279, 3296, 1, 1790, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3253, 3245, 3269, 3298, 3305, 1, 1815, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3287, 3279, 3297, 3309, 3313, 1, 1844, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3420, 3415, 3568, 3669, 3854, 2, 2469, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3576, 3569, 3633, 3856, 3888, 2, 2671, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3677, 3669, 3855, 3905, 3947, 2, 3121, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3863, 3856, 3889, 3949, 3962, 2, 3175, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3913, 3905, 3948, 3969, 3986, 2, 3245, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3977, 3969, 3987, 3994, 3998, 2, 3298, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4040, 4036, 4091, 4146, 4209, 3, 3569, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4153, 4146, 4210, 4236, 4268, 3, 3856, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4243, 4236, 4269, 4279, 4292, 3, 3949, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1, 0, 974, 1590, 1610, 2728, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 979, 975, 1336, 1612, 1661, 3119, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1344, 1337, 1507, 1663, 1715, 3173, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1518, 1508, 1576, 1717, 1750, 3188, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1591, 1590, 1611, 1754, 1764, 3199, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1619, 1612, 1662, 1766, 1789, 3243, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1673, 1663, 1716, 1791, 1814, 3268, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1755, 1754, 1765, 1829, 1833, 3277, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1773, 1766, 1790, 1835, 1843, 3296, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1830, 1829, 1834, 1856, 1857, 3307, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1863, 1862, 2087, 2708, 2729, 3667, 1, 975, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2095, 2088, 2469, 2731, 3120, 3854, 1, 1337, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2480, 2470, 2671, 3122, 3174, 3888, 1, 1508, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2709, 2708, 2730, 3189, 3200, 3903, 1, 1612, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2741, 2731, 3121, 3202, 3244, 3947, 1, 1663, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3132, 3122, 3175, 3246, 3269, 3962, 1, 1717, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3190, 3189, 3201, 3273, 3278, 3967, 1, 1766, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3212, 3202, 3245, 3280, 3297, 3986, 1, 1791, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3274, 3273, 3279, 3306, 3308, 3992, 1, 1835, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3317, 3316, 3415, 3647, 3668, 4144, 2, 2088, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3423, 3416, 3569, 3670, 3855, 4209, 2, 2470, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3648, 3647, 3669, 3893, 3904, 4234, 2, 2731, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3680, 3670, 3856, 3906, 3948, 4268, 2, 3122, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3894, 3893, 3905, 3963, 3968, 4277, 2, 3202, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3916, 3906, 3949, 3970, 3987, 4292, 2, 3246, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3964, 3963, 3969, 3991, 3993, 4294, 2, 3280, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4001, 4000, 4036, 4124, 4145, 4300, 3, 3416, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4125, 4124, 4146, 4224, 4235, 4320, 3, 3670, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4225, 4224, 4236, 4273, 4278, 4330, 3, 3906, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4274, 4273, 4279, 4293, 4295, 4334, 3, 3970, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1521, 1505, 1509, 1714, 1718, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1583, 1574, 1577, 1749, 1751, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1676, 1660, 1664, 1788, 1792, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1727, 1714, 1718, 1813, 1816, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1801, 1788, 1792, 1842, 1845, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1822, 1813, 1816, 1851, 1853, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2483, 2467, 2471, 3119, 3123, 1, 1505, 1509, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2681, 2669, 2672, 3173, 3176, 1, 1574, 1577, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2744, 2728, 2732, 3199, 3203, 1, 1660, 1664, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3135, 3119, 3123, 3243, 3247, 1, 1714, 1718, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3182, 3173, 3176, 3268, 3270, 1, 1749, 1751, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3215, 3199, 3203, 3277, 3281, 1, 1788, 1792, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3256, 3243, 3247, 3296, 3299, 1, 1813, 1816, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3290, 3277, 3281, 3307, 3310, 1, 1842, 1845, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3579, 3567, 3570, 3854, 3857, 2, 2669, 2672, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3640, 3632, 3634, 3888, 3890, 2, 2703, 2705, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3683, 3667, 3671, 3903, 3907, 2, 3119, 3123, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3866, 3854, 3857, 3947, 3950, 2, 3173, 3176, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3919, 3903, 3907, 3967, 3971, 2, 3243, 3247, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3956, 3947, 3950, 3986, 3988, 2, 3268, 3270, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3980, 3967, 3971, 3992, 3995, 2, 3296, 3299, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4098, 4090, 4092, 4209, 4211, 3, 3632, 3634, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4156, 4144, 4147, 4234, 4237, 3, 3854, 3857, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4217, 4209, 4211, 4268, 4270, 3, 3888, 3890, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4246, 4234, 4237, 4277, 4280, 3, 3947, 3950, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4286, 4277, 4280, 4294, 4296, 3, 3986, 3988, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4304, 4300, 4301, 4320, 4321, 4, 4209, 4211, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4324, 4320, 4321, 4330, 4331, 4, 4268, 4270, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1347, 1335, 1338, 1509, 1661, 1667, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1527, 1506, 1512, 1577, 1715, 1721, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1622, 1610, 1613, 1664, 1764, 1767, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1682, 1661, 1667, 1718, 1789, 1795, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1733, 1715, 1721, 1751, 1814, 1819, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1776, 1764, 1767, 1792, 1833, 1836, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1807, 1789, 1795, 1816, 1843, 1848, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2098, 2086, 2089, 2471, 2729, 2735, 1, 1335, 1338, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2489, 2468, 2474, 2672, 3120, 3126, 1, 1506, 1512, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2687, 2670, 2675, 2705, 3174, 3179, 1, 1575, 1580, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2750, 2729, 2735, 3123, 3200, 3206, 1, 1661, 1667, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3141, 3120, 3126, 3176, 3244, 3250, 1, 1715, 1721, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3221, 3200, 3206, 3247, 3278, 3284, 1, 1789, 1795, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3262, 3244, 3250, 3270, 3297, 3302, 1, 1814, 1819, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3426, 3414, 3417, 3570, 3668, 3674, 2, 2468, 2474, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3585, 3568, 3573, 3634, 3855, 3860, 2, 2670, 2675, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3689, 3668, 3674, 3857, 3904, 3910, 2, 3120, 3126, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3872, 3855, 3860, 3890, 3948, 3953, 2, 3174, 3179, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3925, 3904, 3910, 3950, 3968, 3974, 2, 3244, 3250, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4043, 4035, 4037, 4092, 4145, 4150, 3, 3568, 3573, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4104, 4091, 4095, 4121, 4210, 4214, 3, 3633, 3637, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4162, 4145, 4150, 4211, 4235, 4240, 3, 3855, 3860, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4252, 4235, 4240, 4270, 4278, 4283, 3, 3948, 3953, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 982, 974, 976, 1338, 1611, 1616, 1660, 1664, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1353, 1336, 1341, 1512, 1662, 1670, 1714, 1718, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1533, 1507, 1515, 1580, 1716, 1724, 1749, 1751, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1628, 1611, 1616, 1667, 1765, 1770, 1788, 1792, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1688, 1662, 1670, 1721, 1790, 1798, 1813, 1816, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1782, 1765, 1770, 1795, 1834, 1839, 1842, 1845, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2104, 2087, 2092, 2474, 2730, 2738, 3119, 3123, 1, 1336, 1341, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2495, 2469, 2477, 2675, 3121, 3129, 3173, 3176, 1, 1507, 1515, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2756, 2730, 2738, 3126, 3201, 3209, 3243, 3247, 1, 1662, 1670, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3147, 3121, 3129, 3179, 3245, 3253, 3268, 3270, 1, 1716, 1724, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3227, 3201, 3209, 3250, 3279, 3287, 3296, 3299, 1, 1790, 1798, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3432, 3415, 3420, 3573, 3669, 3677, 3854, 3857, 2, 2469, 2477, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3591, 3569, 3576, 3637, 3856, 3863, 3888, 3890, 2, 2671, 2678, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3695, 3669, 3677, 3860, 3905, 3913, 3947, 3950, 2, 3121, 3129, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3931, 3905, 3913, 3953, 3969, 3977, 3986, 3988, 2, 3245, 3253, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4049, 4036, 4040, 4095, 4146, 4153, 4209, 4211, 3, 3569, 3576, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4168, 4146, 4153, 4214, 4236, 4243, 4268, 4270, 3, 3856, 3863, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4, 0, 1, 976, 1590, 1591, 1610, 1613, 2732, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 988, 975, 979, 1341, 1612, 1619, 1661, 1667, 3123, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1359, 1337, 1344, 1515, 1663, 1673, 1715, 1721, 3176, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1594, 1590, 1591, 1616, 1754, 1755, 1764, 1767, 3203, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1634, 1612, 1619, 1670, 1766, 1773, 1789, 1795, 3247, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1758, 1754, 1755, 1770, 1829, 1830, 1833, 1836, 3281, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1866, 1862, 1863, 2092, 2708, 2709, 2729, 2735, 3671, 1, 975, 979, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2110, 2088, 2095, 2477, 2731, 2741, 3120, 3126, 3857, 1, 1337, 1344, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2501, 2470, 2480, 2678, 3122, 3132, 3174, 3179, 3890, 1, 1508, 1518, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2712, 2708, 2709, 2738, 3189, 3190, 3200, 3206, 3907, 1, 1612, 1619, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2762, 2731, 2741, 3129, 3202, 3212, 3244, 3250, 3950, 1, 1663, 1673, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3193, 3189, 3190, 3209, 3273, 3274, 3278, 3284, 3971, 1, 1766, 1773, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3320, 3316, 3317, 3420, 3647, 3648, 3668, 3674, 4147, 2, 2088, 2095, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3438, 3416, 3423, 3576, 3670, 3680, 3855, 3860, 4211, 2, 2470, 2480, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3651, 3647, 3648, 3677, 3893, 3894, 3904, 3910, 4237, 2, 2731, 2741, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3701, 3670, 3680, 3863, 3906, 3916, 3948, 3953, 4270, 2, 3122, 3132, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3897, 3893, 3894, 3913, 3963, 3964, 3968, 3974, 4280, 2, 3202, 3212, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4004, 4000, 4001, 4040, 4124, 4125, 4145, 4150, 4301, 3, 3416, 3423, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4128, 4124, 4125, 4153, 4224, 4225, 4235, 4240, 4321, 3, 3670, 3680, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4228, 4224, 4225, 4243, 4273, 4274, 4278, 4283, 4331, 3, 3906, 3916, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1539, 1509, 1521, 1718, 1727, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1694, 1664, 1676, 1792, 1801, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1739, 1718, 1727, 1816, 1822, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2507, 2471, 2483, 3123, 3135, 1, 1509, 1521, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2693, 2672, 2681, 3176, 3182, 1, 1577, 1583, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2768, 2732, 2744, 3203, 3215, 1, 1664, 1676, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3153, 3123, 3135, 3247, 3256, 1, 1718, 1727, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3233, 3203, 3215, 3281, 3290, 1, 1792, 1801, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3597, 3570, 3579, 3857, 3866, 2, 2672, 2681, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3707, 3671, 3683, 3907, 3919, 2, 3123, 3135, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3878, 3857, 3866, 3950, 3956, 2, 3176, 3182, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3937, 3907, 3919, 3971, 3980, 2, 3247, 3256, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4110, 4092, 4098, 4211, 4217, 3, 3634, 3640, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4174, 4147, 4156, 4237, 4246, 3, 3857, 3866, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4258, 4237, 4246, 4280, 4286, 3, 3950, 3956, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4310, 4301, 4304, 4321, 4324, 4, 4211, 4217, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1365, 1338, 1347, 1521, 1667, 1682, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1549, 1512, 1527, 1583, 1721, 1733, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1640, 1613, 1622, 1676, 1767, 1776, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1704, 1667, 1682, 1727, 1795, 1807, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 2116, 2089, 2098, 2483, 2735, 2750, 1, 1338, 1347, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 2517, 2474, 2489, 2681, 3126, 3141, 1, 1512, 1527, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 2778, 2735, 2750, 3135, 3206, 3221, 1, 1667, 1682, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 3163, 3126, 3141, 3182, 3250, 3262, 1, 1721, 1733, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 3444, 3417, 3426, 3579, 3674, 3689, 2, 2474, 2489, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 3607, 3573, 3585, 3640, 3860, 3872, 2, 2675, 2687, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 3717, 3674, 3689, 3866, 3910, 3925, 2, 3126, 3141, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4055, 4037, 4043, 4098, 4150, 4162, 3, 3573, 3585, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4184, 4150, 4162, 4217, 4240, 4252, 3, 3860, 3872, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 994, 976, 982, 1347, 1616, 1628, 1664, 1676, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 1375, 1341, 1353, 1527, 1670, 1688, 1718, 1727, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 1650, 1616, 1628, 1682, 1770, 1782, 1792, 1801, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2126, 2092, 2104, 2489, 2738, 2756, 3123, 3135, 1, 1341, 1353, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2527, 2477, 2495, 2687, 3129, 3147, 3176, 3182, 1, 1515, 1533, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2788, 2738, 2756, 3141, 3209, 3227, 3247, 3256, 1, 1670, 1688, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 3454, 3420, 3432, 3585, 3677, 3695, 3857, 3866, 2, 2477, 2495, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 3727, 3677, 3695, 3872, 3913, 3931, 3950, 3956, 2, 3129, 3147, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4065, 4040, 4049, 4104, 4153, 4168, 4211, 4217, 3, 3576, 3591, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 10, 1, 4, 982, 1591, 1594, 1613, 1622, 2744, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1004, 979, 988, 1353, 1619, 1634, 1667, 1682, 3135, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1600, 1591, 1594, 1628, 1755, 1758, 1767, 1776, 3215, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1872, 1863, 1866, 2104, 2709, 2712, 2735, 2750, 3683, 1, 979, 988, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 2136, 2095, 2110, 2495, 2741, 2762, 3126, 3141, 3866, 1, 1344, 1359, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 2718, 2709, 2712, 2756, 3190, 3193, 3206, 3221, 3919, 1, 1619, 1634, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 3326, 3317, 3320, 3432, 3648, 3651, 3674, 3689, 4156, 2, 2095, 2110, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 3464, 3423, 3438, 3591, 3680, 3701, 3860, 3872, 4217, 2, 2480, 2501, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 3657, 3648, 3651, 3695, 3894, 3897, 3910, 3925, 4246, 2, 2741, 2762, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 4010, 4001, 4004, 4049, 4125, 4128, 4150, 4162, 4304, 3, 3423, 3438, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 4134, 4125, 4128, 4168, 4225, 4228, 4240, 4252, 4324, 3, 3680, 3701, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 1559, 1521, 1539, 1727, 1739, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 2537, 2483, 2507, 3135, 3153, 1, 1521, 1539, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 2798, 2744, 2768, 3215, 3233, 1, 1676, 1694, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 3617, 3579, 3597, 3866, 3878, 2, 2681, 2693, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 3737, 3683, 3707, 3919, 3937, 2, 3135, 3153, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4194, 4156, 4174, 4246, 4258, 3, 3866, 3878, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 1385, 1347, 1365, 1539, 1682, 1704, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 2146, 2098, 2116, 2507, 2750, 2778, 1, 1347, 1365, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 2552, 2489, 2517, 2693, 3141, 3163, 1, 1527, 1549, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 3474, 3426, 3444, 3597, 3689, 3717, 2, 2489, 2517, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 4075, 4043, 4055, 4110, 4162, 4184, 3, 3585, 3607, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 1014, 982, 994, 1365, 1628, 1650, 1676, 1694, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 2161, 2104, 2126, 2517, 2756, 2788, 3135, 3153, 1, 1353, 1375, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 3489, 3432, 3454, 3607, 3695, 3727, 3866, 3878, 2, 2495, 2527, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 20, 4, 10, 994, 1594, 1600, 1622, 1640, 2768, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 1882, 1866, 1872, 2126, 2712, 2718, 2750, 2778, 3707, 1, 988, 1004, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 3336, 3320, 3326, 3454, 3651, 3657, 3689, 3717, 4174, 2, 2110, 2136, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 4020, 4004, 4010, 4065, 4128, 4134, 4162, 4184, 4310, 3, 3438, 3464, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_p(pbuffer, 1400, 1335, 1505, 2468, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 2567, 2471, 3570, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 2813, 2732, 3671, 1, 2471, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 3752, 3671, 4147, 1, 3570, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 1403, 1338, 1505, 1509, 2474, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_d(pbuffer, 1029, 976, 1335, 1338, 2092, 2471, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_d(pbuffer, 2176, 2092, 2468, 2474, 3420, 3570, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 2576, 2483, 3579, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 2822, 2744, 3683, 1, 2483, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 3761, 3683, 4156, 1, 3579, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 1412, 1347, 1509, 1521, 2489, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 2185, 2098, 2471, 2483, 3426, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 3504, 3426, 3570, 3579, 4043, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 1038, 982, 1338, 1347, 2104, 2483, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 2203, 2104, 2474, 2489, 3432, 3579, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 35, 4, 976, 982, 1866, 2098, 2732, 2744, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 1897, 1866, 2092, 2104, 3320, 3426, 3671, 3683, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 3351, 3320, 3420, 3432, 4004, 4043, 4147, 4156, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 2594, 2507, 3597, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 2840, 2768, 3707, 1, 2507, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 3779, 3707, 4174, 1, 3597, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 1430, 1365, 1521, 1539, 2517, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 1056, 994, 1347, 1365, 2126, 2507, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 2221, 2126, 2489, 2517, 3454, 3597, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 2624, 2537, 3617, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 2870, 2798, 3737, 1, 2537, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 3809, 3737, 4194, 1, 3617, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 1460, 1385, 1539, 1559, 2552, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 2251, 2146, 2507, 2537, 3474, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 3522, 3474, 3597, 3617, 4075, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 1086, 1014, 1365, 1385, 2161, 2537, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 2296, 2161, 2517, 2552, 3489, 3617, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 53, 20, 994, 1014, 1882, 2146, 2768, 2798, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 1915, 1882, 2126, 2161, 3336, 3474, 3707, 3737, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 3369, 3336, 3454, 3489, 4020, 4075, 4174, 4194, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_s(pbuffer, 2915, 2732, 2813, 3671, 3752, 1, 2471, 2567, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_d(pbuffer, 1131, 976, 1029, 1400, 1403, 2092, 2176, 2471, 2567, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 2933, 2744, 2822, 3683, 3761, 1, 2483, 2576, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 2341, 2098, 2185, 2567, 2576, 3426, 3504, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_d(pbuffer, 1149, 982, 1038, 1403, 1412, 2104, 2203, 2483, 2576, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_f(pbuffer, 98, 4, 35, 1029, 1038, 1866, 1897, 2098, 2185, 2813, 2822, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_f(pbuffer, 1960, 1866, 1897, 2176, 2203, 3320, 3351, 3426, 3504, 3752, 3761, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 2969, 2768, 2840, 3707, 3779, 1, 2507, 2594, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 1185, 994, 1056, 1412, 1430, 2126, 2221, 2507, 2594, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 3029, 2798, 2870, 3737, 3809, 1, 2537, 2624, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 2377, 2146, 2251, 2594, 2624, 3474, 3522, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 1245, 1014, 1086, 1430, 1460, 2161, 2296, 2537, 2624, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_f(pbuffer, 674, 20, 53, 1056, 1086, 1882, 1915, 2146, 2251, 2840, 2870, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_f(pbuffer, 1996, 1882, 1915, 2221, 2296, 3336, 3369, 3474, 3522, 3779, 3809, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_f(pbuffer, 764, 35, 98, 1131, 1149, 1897, 1960, 2185, 2341, 2915, 2933, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_f(pbuffer, 824, 53, 674, 1185, 1245, 1915, 1996, 2251, 2377, 2969, 3029, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2cgeom::comp_prim_op_geom_11_df(pbuffer, 134, 35, 53, 764, 824, 1, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, 0, pbuffer, 134, 540, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<2, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2pecp namespace

#endif /* ProjectedCorePotentialGeom101DFForF_hpp */
