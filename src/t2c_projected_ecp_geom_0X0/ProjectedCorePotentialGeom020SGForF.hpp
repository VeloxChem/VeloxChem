#ifndef ProjectedCorePotentialGeom020SGForF_hpp
#define ProjectedCorePotentialGeom020SGForF_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "SimdArray.hpp"
#include "ProjectedCorePotentialPrimRecDGForF.hpp"
#include "ProjectedCorePotentialPrimRecPFForD.hpp"
#include "ProjectedCorePotentialPrimRecPFForF.hpp"
#include "ProjectedCorePotentialPrimRecPFForS.hpp"
#include "ProjectedCorePotentialPrimRecPGForD.hpp"
#include "ProjectedCorePotentialPrimRecPGForF.hpp"
#include "ProjectedCorePotentialPrimRecPGForP.hpp"
#include "ProjectedCorePotentialPrimRecPGForS.hpp"
#include "ProjectedCorePotentialPrimRecPHForF.hpp"
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
#include "ProjectedCorePotentialPrimRecSHForD.hpp"
#include "ProjectedCorePotentialPrimRecSHForF.hpp"
#include "ProjectedCorePotentialPrimRecSHForP.hpp"
#include "ProjectedCorePotentialPrimRecSHForS.hpp"
#include "ProjectedCorePotentialPrimRecSIForF.hpp"
#include "ProjectedCorePotentialPrimRecSPForD.hpp"
#include "ProjectedCorePotentialPrimRecSPForF.hpp"
#include "ProjectedCorePotentialPrimRecSPForP.hpp"
#include "ProjectedCorePotentialPrimRecSPForS.hpp"
#include "ProjectedCorePotentialPrimRecSS.hpp"
#include "GeometricalDerivatives020ForSG.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes (S|U_l|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_projected_core_potential_geom_020_sg_for_f(T& distributor,
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

    CSimdArray<double> i_values(10, ket_npgtos);

    CSimdArray<double> l_values(4, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(2300, ket_npgtos);

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

                    t2cfunc::comp_i_vals(i_values, 9, pfactors, 8);

                    t2cfunc::comp_l_vals(l_values, 3, pfactors, 8, 6);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 0, pbuffer, 0, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 1, pbuffer, 402, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 1, pbuffer, 403, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 2, pbuffer, 568, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 2, pbuffer, 569, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 2, pbuffer, 570, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 3, pbuffer, 633, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 3, pbuffer, 634, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 3, pbuffer, 635, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 4, pbuffer, 667, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 4, pbuffer, 668, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 0, pbuffer, 672, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 1, pbuffer, 728, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 1, pbuffer, 729, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 1, pbuffer, 730, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 2, pbuffer, 839, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 2, pbuffer, 840, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 2, pbuffer, 841, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 2, pbuffer, 842, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 3, pbuffer, 908, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 3, pbuffer, 909, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 3, pbuffer, 910, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 4, pbuffer, 923, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 0, pbuffer, 924, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 1, pbuffer, 959, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 1, pbuffer, 960, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 1, pbuffer, 961, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 2, pbuffer, 1024, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 2, pbuffer, 1025, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 2, pbuffer, 1026, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 2, pbuffer, 1027, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 3, pbuffer, 1059, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 3, pbuffer, 1060, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 0, pbuffer, 1064, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 1, pbuffer, 1084, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 1, pbuffer, 1085, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 1, pbuffer, 1086, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 2, pbuffer, 1118, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 2, pbuffer, 1119, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 2, pbuffer, 1120, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 3, pbuffer, 1133, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 0, pbuffer, 1134, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 1, pbuffer, 1144, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 1, pbuffer, 1145, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 1, pbuffer, 1146, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 2, pbuffer, 1159, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 2, pbuffer, 1160, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 5, 0, pbuffer, 1164, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 5, 1, pbuffer, 1168, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 5, 1, pbuffer, 1169, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 5, 2, pbuffer, 1173, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 6, 0, pbuffer, 1174, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 6, 1, pbuffer, 1175, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 0, pbuffer, 1176, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 1, pbuffer, 1277, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 1, pbuffer, 1278, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 1, pbuffer, 1279, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 2, pbuffer, 1433, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 2, pbuffer, 1434, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 2, pbuffer, 1435, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 3, pbuffer, 1498, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 3, pbuffer, 1499, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 4, pbuffer, 1512, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 0, pbuffer, 1513, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 1, pbuffer, 1548, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 1, pbuffer, 1549, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 1, pbuffer, 1550, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 1, pbuffer, 1551, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 2, pbuffer, 1744, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 2, pbuffer, 1745, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 2, pbuffer, 1746, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 3, pbuffer, 1778, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 3, pbuffer, 1779, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 0, pbuffer, 1783, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 1, pbuffer, 1803, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 1, pbuffer, 1804, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 1, pbuffer, 1805, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 1, pbuffer, 1806, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 2, pbuffer, 1872, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 2, pbuffer, 1873, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 2, pbuffer, 1874, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 3, pbuffer, 1887, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 0, pbuffer, 1888, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 1, pbuffer, 1898, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 1, pbuffer, 1899, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 1, pbuffer, 1900, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 1, pbuffer, 1901, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 2, pbuffer, 1933, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 2, pbuffer, 1934, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 4, 0, pbuffer, 1938, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 1, pbuffer, 1942, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 1, pbuffer, 1943, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 1, pbuffer, 1944, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 2, pbuffer, 1957, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 5, 0, pbuffer, 1958, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 1, pbuffer, 1959, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 1, pbuffer, 1960, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 6, 1, pbuffer, 1964, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 0, pbuffer, 1965, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 1, pbuffer, 2000, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 1, pbuffer, 2001, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 2, pbuffer, 2055, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 2, pbuffer, 2056, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 3, pbuffer, 2085, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 0, pbuffer, 2089, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 1, pbuffer, 2109, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 1, pbuffer, 2110, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 1, pbuffer, 2111, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 2, pbuffer, 2174, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 2, pbuffer, 2175, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 3, pbuffer, 2188, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 0, pbuffer, 2189, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 1, pbuffer, 2199, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 1, pbuffer, 2200, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 1, pbuffer, 2201, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 2, pbuffer, 2233, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 2, pbuffer, 2234, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 0, pbuffer, 2238, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 1, pbuffer, 2242, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 1, pbuffer, 2243, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 1, pbuffer, 2244, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 2, pbuffer, 2257, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 4, 0, pbuffer, 2258, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 1, pbuffer, 2259, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 1, pbuffer, 2260, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 1, pbuffer, 2264, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 1, pbuffer, 2265, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 1, pbuffer, 2285, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 1, pbuffer, 2295, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 1, pbuffer, 2299, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 636, 633, 908, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 669, 667, 923, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 843, 839, 1024, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 911, 908, 1059, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1028, 1024, 1118, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1061, 1059, 1133, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1121, 1118, 1159, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1161, 1159, 1173, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1436, 1433, 1744, 1, 633, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1500, 1498, 1778, 1, 667, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1552, 1548, 1803, 1, 839, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1747, 1744, 1872, 1, 908, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1780, 1778, 1887, 1, 923, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1807, 1803, 1898, 1, 1024, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1875, 1872, 1933, 1, 1059, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1902, 1898, 1942, 1, 1118, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1935, 1933, 1957, 1, 1133, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1945, 1942, 1959, 1, 1159, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1961, 1959, 1964, 1, 1173, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2057, 2055, 2174, 2, 1498, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2086, 2085, 2188, 2, 1512, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2112, 2109, 2199, 2, 1744, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2176, 2174, 2233, 2, 1778, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2202, 2199, 2242, 2, 1872, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2235, 2233, 2257, 2, 1887, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2245, 2242, 2259, 2, 1933, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2261, 2259, 2264, 2, 1957, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2266, 2265, 2285, 3, 2174, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2286, 2285, 2295, 3, 2233, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2296, 2295, 2299, 3, 2257, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 571, 568, 633, 840, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 639, 634, 667, 909, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 731, 728, 839, 959, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 846, 840, 908, 1025, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 914, 909, 923, 1060, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 962, 959, 1024, 1084, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1031, 1025, 1059, 1119, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1087, 1084, 1118, 1144, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1124, 1119, 1133, 1160, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1147, 1144, 1159, 1168, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1170, 1168, 1173, 1175, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1280, 1277, 1433, 1549, 1, 568, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1439, 1434, 1498, 1745, 1, 634, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1503, 1499, 1512, 1779, 1, 668, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1555, 1549, 1744, 1804, 1, 840, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1750, 1745, 1778, 1873, 1, 909, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1810, 1804, 1872, 1899, 1, 1025, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1878, 1873, 1887, 1934, 1, 1060, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1905, 1899, 1933, 1943, 1, 1119, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1948, 1943, 1957, 1960, 1, 1160, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2002, 2000, 2055, 2110, 2, 1434, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2060, 2056, 2085, 2175, 2, 1499, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2115, 2110, 2174, 2200, 2, 1745, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2179, 2175, 2188, 2234, 2, 1779, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2205, 2200, 2233, 2243, 2, 1873, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 2248, 2243, 2257, 2260, 2, 1934, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 404, 402, 568, 729, 839, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 574, 569, 634, 841, 908, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 642, 635, 668, 910, 923, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 734, 729, 840, 960, 1024, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 849, 841, 909, 1026, 1059, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 965, 960, 1025, 1085, 1118, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1034, 1026, 1060, 1120, 1133, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1090, 1085, 1119, 1145, 1159, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1150, 1145, 1160, 1169, 1173, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1283, 1278, 1434, 1550, 1744, 1, 569, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1442, 1435, 1499, 1746, 1778, 1, 635, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1558, 1550, 1745, 1805, 1872, 1, 841, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1753, 1746, 1779, 1874, 1887, 1, 910, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1813, 1805, 1873, 1900, 1933, 1, 1026, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1908, 1900, 1934, 1944, 1957, 1, 1120, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2005, 2001, 2056, 2111, 2174, 2, 1435, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2118, 2111, 2175, 2201, 2233, 2, 1746, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2208, 2201, 2234, 2244, 2257, 2, 1874, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1, 0, 402, 672, 728, 1548, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 407, 403, 569, 730, 840, 1744, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 577, 570, 635, 842, 909, 1778, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 673, 672, 729, 924, 959, 1803, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 737, 730, 841, 961, 1025, 1872, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 852, 842, 910, 1027, 1060, 1887, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 925, 924, 960, 1064, 1084, 1898, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 968, 961, 1026, 1086, 1119, 1933, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1065, 1064, 1085, 1134, 1144, 1942, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1093, 1086, 1120, 1146, 1160, 1957, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1135, 1134, 1145, 1164, 1168, 1959, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1165, 1164, 1169, 1174, 1175, 1964, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1177, 1176, 1278, 1513, 1549, 2109, 1, 403, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1286, 1279, 1435, 1551, 1745, 2174, 1, 570, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1514, 1513, 1550, 1783, 1804, 2199, 1, 730, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1561, 1551, 1746, 1806, 1873, 2233, 1, 842, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1784, 1783, 1805, 1888, 1899, 2242, 1, 961, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1816, 1806, 1874, 1901, 1934, 2257, 1, 1027, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1889, 1888, 1900, 1938, 1943, 2259, 1, 1086, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1939, 1938, 1944, 1958, 1960, 2264, 1, 1146, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1966, 1965, 2001, 2089, 2110, 2265, 2, 1279, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2090, 2089, 2111, 2189, 2200, 2285, 2, 1551, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2190, 2189, 2201, 2238, 2243, 2295, 2, 1806, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2239, 2238, 2244, 2258, 2260, 2299, 2, 1901, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 645, 633, 636, 908, 911, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 855, 839, 843, 1024, 1028, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 917, 908, 911, 1059, 1061, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1037, 1024, 1028, 1118, 1121, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1127, 1118, 1121, 1159, 1161, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1445, 1433, 1436, 1744, 1747, 1, 633, 636, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1506, 1498, 1500, 1778, 1780, 1, 667, 669, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1564, 1548, 1552, 1803, 1807, 1, 839, 843, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1756, 1744, 1747, 1872, 1875, 1, 908, 911, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1819, 1803, 1807, 1898, 1902, 1, 1024, 1028, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1881, 1872, 1875, 1933, 1935, 1, 1059, 1061, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1911, 1898, 1902, 1942, 1945, 1, 1118, 1121, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1951, 1942, 1945, 1959, 1961, 1, 1159, 1161, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2063, 2055, 2057, 2174, 2176, 2, 1498, 1500, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2121, 2109, 2112, 2199, 2202, 2, 1744, 1747, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2182, 2174, 2176, 2233, 2235, 2, 1778, 1780, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2211, 2199, 2202, 2242, 2245, 2, 1872, 1875, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2251, 2242, 2245, 2259, 2261, 2, 1933, 1935, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2269, 2265, 2266, 2285, 2286, 3, 2174, 2176, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2289, 2285, 2286, 2295, 2296, 3, 2233, 2235, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 580, 568, 571, 636, 840, 846, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 651, 634, 639, 669, 909, 914, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 740, 728, 731, 843, 959, 962, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 861, 840, 846, 911, 1025, 1031, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 971, 959, 962, 1028, 1084, 1087, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1043, 1025, 1031, 1061, 1119, 1124, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1096, 1084, 1087, 1121, 1144, 1147, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1153, 1144, 1147, 1161, 1168, 1170, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1289, 1277, 1280, 1436, 1549, 1555, 1, 568, 571, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1451, 1434, 1439, 1500, 1745, 1750, 1, 634, 639, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1570, 1549, 1555, 1747, 1804, 1810, 1, 840, 846, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1762, 1745, 1750, 1780, 1873, 1878, 1, 909, 914, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1825, 1804, 1810, 1875, 1899, 1905, 1, 1025, 1031, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1917, 1899, 1905, 1935, 1943, 1948, 1, 1119, 1124, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2008, 2000, 2002, 2057, 2110, 2115, 2, 1434, 1439, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2069, 2056, 2060, 2086, 2175, 2179, 2, 1499, 1503, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2127, 2110, 2115, 2176, 2200, 2205, 2, 1745, 1750, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 2217, 2200, 2205, 2235, 2243, 2248, 2, 1873, 1878, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 410, 402, 404, 571, 729, 734, 839, 843, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 586, 569, 574, 639, 841, 849, 908, 911, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 746, 729, 734, 846, 960, 965, 1024, 1028, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 867, 841, 849, 914, 1026, 1034, 1059, 1061, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 977, 960, 965, 1031, 1085, 1090, 1118, 1121, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1102, 1085, 1090, 1124, 1145, 1150, 1159, 1161, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1295, 1278, 1283, 1439, 1550, 1558, 1744, 1747, 1, 569, 574, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1457, 1435, 1442, 1503, 1746, 1753, 1778, 1780, 1, 635, 642, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1576, 1550, 1558, 1750, 1805, 1813, 1872, 1875, 1, 841, 849, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1831, 1805, 1813, 1878, 1900, 1908, 1933, 1935, 1, 1026, 1034, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2014, 2001, 2005, 2060, 2111, 2118, 2174, 2176, 2, 1435, 1442, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2133, 2111, 2118, 2179, 2201, 2208, 2233, 2235, 2, 1746, 1753, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4, 0, 1, 404, 672, 673, 728, 731, 1552, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 416, 403, 407, 574, 730, 737, 840, 846, 1747, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 592, 570, 577, 642, 842, 852, 909, 914, 1780, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 676, 672, 673, 734, 924, 925, 959, 962, 1807, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 752, 730, 737, 849, 961, 968, 1025, 1031, 1875, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 928, 924, 925, 965, 1064, 1065, 1084, 1087, 1902, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 983, 961, 968, 1034, 1086, 1093, 1119, 1124, 1935, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1068, 1064, 1065, 1090, 1134, 1135, 1144, 1147, 1945, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1138, 1134, 1135, 1150, 1164, 1165, 1168, 1170, 1961, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1180, 1176, 1177, 1283, 1513, 1514, 1549, 1555, 2112, 1, 403, 407, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1301, 1279, 1286, 1442, 1551, 1561, 1745, 1750, 2176, 1, 570, 577, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1517, 1513, 1514, 1558, 1783, 1784, 1804, 1810, 2202, 1, 730, 737, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1582, 1551, 1561, 1753, 1806, 1816, 1873, 1878, 2235, 1, 842, 852, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1787, 1783, 1784, 1813, 1888, 1889, 1899, 1905, 2245, 1, 961, 968, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1892, 1888, 1889, 1908, 1938, 1939, 1943, 1948, 2261, 1, 1086, 1093, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1969, 1965, 1966, 2005, 2089, 2090, 2110, 2115, 2266, 2, 1279, 1286, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2093, 2089, 2090, 2118, 2189, 2190, 2200, 2205, 2286, 2, 1551, 1561, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2193, 2189, 2190, 2208, 2238, 2239, 2243, 2248, 2296, 2, 1806, 1816, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 657, 636, 645, 911, 917, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 873, 843, 855, 1028, 1037, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1049, 1028, 1037, 1121, 1127, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1463, 1436, 1445, 1747, 1756, 1, 636, 645, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1588, 1552, 1564, 1807, 1819, 1, 843, 855, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1768, 1747, 1756, 1875, 1881, 1, 911, 917, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1837, 1807, 1819, 1902, 1911, 1, 1028, 1037, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1923, 1902, 1911, 1945, 1951, 1, 1121, 1127, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2075, 2057, 2063, 2176, 2182, 2, 1500, 1506, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2139, 2112, 2121, 2202, 2211, 2, 1747, 1756, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2223, 2202, 2211, 2245, 2251, 2, 1875, 1881, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2275, 2266, 2269, 2286, 2289, 3, 2176, 2182, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 598, 571, 580, 645, 846, 861, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 758, 731, 740, 855, 962, 971, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 883, 846, 861, 917, 1031, 1043, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 989, 962, 971, 1037, 1087, 1096, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1108, 1087, 1096, 1127, 1147, 1153, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1307, 1280, 1289, 1445, 1555, 1570, 1, 571, 580, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1473, 1439, 1451, 1506, 1750, 1762, 1, 639, 651, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1598, 1555, 1570, 1756, 1810, 1825, 1, 846, 861, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1847, 1810, 1825, 1881, 1905, 1917, 1, 1031, 1043, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 2020, 2002, 2008, 2063, 2115, 2127, 2, 1439, 1451, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 2149, 2115, 2127, 2182, 2205, 2217, 2, 1750, 1762, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 422, 404, 410, 580, 734, 746, 843, 855, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 608, 574, 586, 651, 849, 867, 911, 917, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 768, 734, 746, 861, 965, 977, 1028, 1037, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 999, 965, 977, 1043, 1090, 1102, 1121, 1127, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 1317, 1283, 1295, 1451, 1558, 1576, 1747, 1756, 1, 574, 586, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 1608, 1558, 1576, 1762, 1813, 1831, 1875, 1881, 1, 849, 867, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2030, 2005, 2014, 2069, 2118, 2133, 2176, 2182, 2, 1442, 1457, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 10, 1, 4, 410, 673, 676, 731, 740, 1564, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 432, 407, 416, 586, 737, 752, 846, 861, 1756, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 682, 673, 676, 746, 925, 928, 962, 971, 1819, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 778, 737, 752, 867, 968, 983, 1031, 1043, 1881, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 934, 925, 928, 977, 1065, 1068, 1087, 1096, 1911, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1074, 1065, 1068, 1102, 1135, 1138, 1147, 1153, 1951, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1186, 1177, 1180, 1295, 1514, 1517, 1555, 1570, 2121, 1, 407, 416, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1327, 1286, 1301, 1457, 1561, 1582, 1750, 1762, 2182, 1, 577, 592, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1523, 1514, 1517, 1576, 1784, 1787, 1810, 1825, 2211, 1, 737, 752, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1793, 1784, 1787, 1831, 1889, 1892, 1905, 1917, 2251, 1, 968, 983, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 1975, 1966, 1969, 2014, 2090, 2093, 2115, 2127, 2269, 2, 1286, 1301, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 2099, 2090, 2093, 2133, 2190, 2193, 2205, 2217, 2289, 2, 1561, 1582, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 893, 855, 873, 1037, 1049, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 1483, 1445, 1463, 1756, 1768, 1, 645, 657, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 1618, 1564, 1588, 1819, 1837, 1, 855, 873, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 1857, 1819, 1837, 1911, 1923, 1, 1037, 1049, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 2159, 2121, 2139, 2211, 2223, 2, 1756, 1768, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 618, 580, 598, 657, 861, 883, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 788, 740, 758, 873, 971, 989, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 1009, 971, 989, 1049, 1096, 1108, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 1337, 1289, 1307, 1463, 1570, 1598, 1, 580, 598, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 1633, 1570, 1598, 1768, 1825, 1847, 1, 861, 883, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 2040, 2008, 2020, 2075, 2127, 2149, 2, 1451, 1473, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 442, 410, 422, 598, 746, 768, 855, 873, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 803, 746, 768, 883, 977, 999, 1037, 1049, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 1352, 1295, 1317, 1473, 1576, 1608, 1756, 1768, 1, 586, 608, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 20, 4, 10, 422, 676, 682, 740, 758, 1588, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 457, 416, 432, 608, 752, 778, 861, 883, 1768, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 692, 676, 682, 768, 928, 934, 971, 989, 1837, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 944, 928, 934, 999, 1068, 1074, 1096, 1108, 1923, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 1196, 1180, 1186, 1317, 1517, 1523, 1570, 1598, 2139, 1, 416, 432, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 1533, 1517, 1523, 1608, 1787, 1793, 1825, 1847, 2223, 1, 752, 778, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 1985, 1969, 1975, 2030, 2093, 2099, 2127, 2149, 2275, 2, 1301, 1327, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 1648, 1588, 1618, 1837, 1857, 1, 873, 893, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 818, 758, 788, 893, 989, 1009, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 1367, 1307, 1337, 1483, 1598, 1633, 1, 598, 618, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 472, 422, 442, 618, 768, 803, 873, 893, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 125, 10, 20, 442, 682, 692, 758, 788, 1618, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 707, 682, 692, 803, 934, 944, 989, 1009, 1857, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 1211, 1186, 1196, 1352, 1523, 1533, 1598, 1633, 2159, 1, 432, 457, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_si_f(pbuffer, 146, 20, 125, 472, 692, 707, 788, 818, 1648, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 1669, 1588, 2139, 1, 1463, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 493, 422, 580, 598, 1317, 1463, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 174, 10, 410, 422, 1186, 1307, 1564, 1588, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 1699, 1618, 2159, 1, 1483, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 1388, 1337, 1463, 1483, 2040, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 523, 442, 598, 618, 1352, 1483, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 204, 20, 422, 442, 1196, 1337, 1588, 1618, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_f(pbuffer, 1232, 1196, 1317, 1352, 1985, 2040, 2139, 2159, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 249, 125, 442, 472, 1211, 1367, 1618, 1648, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_f(pbuffer, 312, 20, 204, 493, 523, 1196, 1232, 1337, 1388, 1669, 1699, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2cgeom::comp_prim_op_geom_020_sg(pbuffer, 35, 4, 20, 146, 174, 249, 312, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, 0, pbuffer, 35, 90, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<0, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 0, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2pecp namespace

#endif /* ProjectedCorePotentialGeom020SGForF_hpp */
