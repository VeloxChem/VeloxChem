#ifndef ProjectedCorePotentialGeom010PFForG_hpp
#define ProjectedCorePotentialGeom010PFForG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "SimdArray.hpp"
#include "ProjectedCorePotentialPrimRecDFForG.hpp"
#include "ProjectedCorePotentialPrimRecPDForF.hpp"
#include "ProjectedCorePotentialPrimRecPDForG.hpp"
#include "ProjectedCorePotentialPrimRecPDForP.hpp"
#include "ProjectedCorePotentialPrimRecPFForD.hpp"
#include "ProjectedCorePotentialPrimRecPFForF.hpp"
#include "ProjectedCorePotentialPrimRecPFForG.hpp"
#include "ProjectedCorePotentialPrimRecPFForP.hpp"
#include "ProjectedCorePotentialPrimRecPFForS.hpp"
#include "ProjectedCorePotentialPrimRecPGForG.hpp"
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
#include "GeometricalDerivatives010ForPF.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes (P|U_l|F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_projected_core_potential_geom_010_pf_for_g(T& distributor,
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

    CSimdArray<double> l_values(5, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(1499, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(90, 1);

    CSimdArray<double> sbuffer(63, 1);

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

                    t2cfunc::comp_l_vals(l_values, 4, pfactors, 8, 6);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 0, pbuffer, 0, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 1, pbuffer, 278, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 1, pbuffer, 279, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 2, pbuffer, 381, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 2, pbuffer, 382, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 0, 2, pbuffer, 383, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 3, pbuffer, 415, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 3, pbuffer, 416, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 3, pbuffer, 417, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 4, pbuffer, 430, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 4, pbuffer, 431, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 0, pbuffer, 435, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 1, pbuffer, 455, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 1, pbuffer, 456, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 1, pbuffer, 457, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 2, pbuffer, 505, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 2, pbuffer, 506, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 2, pbuffer, 507, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 1, 2, pbuffer, 508, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 3, pbuffer, 540, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 3, pbuffer, 541, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 3, pbuffer, 542, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 4, pbuffer, 555, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 0, pbuffer, 556, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 1, pbuffer, 566, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 1, pbuffer, 567, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 2, 1, pbuffer, 568, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 2, pbuffer, 590, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 2, pbuffer, 591, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 2, pbuffer, 592, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 2, pbuffer, 593, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 3, pbuffer, 615, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 3, pbuffer, 616, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 3, 0, pbuffer, 620, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 1, pbuffer, 624, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 1, pbuffer, 625, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 3, 1, pbuffer, 626, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 2, pbuffer, 633, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 2, pbuffer, 634, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 2, pbuffer, 635, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 3, pbuffer, 642, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 0, 4, 0, pbuffer, 643, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 1, pbuffer, 644, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 1, pbuffer, 645, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 2, pbuffer, 646, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 2, pbuffer, 647, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 0, pbuffer, 648, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 1, pbuffer, 713, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 1, pbuffer, 714, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 0, 1, pbuffer, 715, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 2, pbuffer, 808, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 2, pbuffer, 809, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 2, pbuffer, 810, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 3, pbuffer, 842, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 3, pbuffer, 843, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 3, pbuffer, 844, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 4, pbuffer, 857, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 0, pbuffer, 858, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 1, pbuffer, 878, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 1, pbuffer, 879, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 1, pbuffer, 880, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 1, 1, pbuffer, 881, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 2, pbuffer, 995, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 2, pbuffer, 996, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 2, pbuffer, 997, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 2, pbuffer, 998, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 3, pbuffer, 1030, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 3, pbuffer, 1031, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 0, pbuffer, 1035, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 1, pbuffer, 1045, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 1, pbuffer, 1046, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 1, pbuffer, 1047, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 1, pbuffer, 1048, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 2, 1, pbuffer, 1049, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 2, pbuffer, 1100, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 2, pbuffer, 1101, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 2, pbuffer, 1102, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 3, pbuffer, 1115, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 3, 0, pbuffer, 1116, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 1, pbuffer, 1120, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 1, pbuffer, 1121, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 1, pbuffer, 1122, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 1, pbuffer, 1123, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 2, pbuffer, 1145, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 2, pbuffer, 1146, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 1, 4, 0, pbuffer, 1150, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 1, pbuffer, 1151, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 1, pbuffer, 1152, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 1, pbuffer, 1153, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 2, pbuffer, 1160, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 1, pbuffer, 1161, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 1, pbuffer, 1162, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 0, 0, pbuffer, 1163, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 1, pbuffer, 1183, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 1, pbuffer, 1184, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 2, pbuffer, 1213, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 2, pbuffer, 1214, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 2, pbuffer, 1215, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 3, pbuffer, 1247, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 3, pbuffer, 1248, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 1, 0, pbuffer, 1252, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 1, pbuffer, 1262, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 1, pbuffer, 1263, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 1, pbuffer, 1264, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 1, pbuffer, 1265, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 2, pbuffer, 1361, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 2, pbuffer, 1362, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 2, pbuffer, 1363, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 3, pbuffer, 1376, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 2, 0, pbuffer, 1377, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 1, pbuffer, 1381, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 1, pbuffer, 1382, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 1, pbuffer, 1383, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 1, pbuffer, 1384, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 2, pbuffer, 1416, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 2, pbuffer, 1417, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(4, 2, 3, 0, pbuffer, 1421, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 1, pbuffer, 1422, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 1, pbuffer, 1423, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 1, pbuffer, 1424, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 2, pbuffer, 1437, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 1, pbuffer, 1438, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 1, pbuffer, 1439, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 1, pbuffer, 1443, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 1, pbuffer, 1444, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 1, pbuffer, 1445, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 2, pbuffer, 1474, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 1, pbuffer, 1478, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 1, pbuffer, 1479, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 2, pbuffer, 1492, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 1, pbuffer, 1493, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 1, pbuffer, 1494, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 1, pbuffer, 1498, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 432, 430, 555, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 543, 540, 615, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 594, 590, 633, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 617, 615, 642, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 636, 633, 646, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 845, 842, 1030, 1, 430, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 999, 995, 1100, 1, 540, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1032, 1030, 1115, 1, 555, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1050, 1045, 1120, 1, 590, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1103, 1100, 1145, 1, 615, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1124, 1120, 1151, 1, 633, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1147, 1145, 1160, 1, 642, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1154, 1151, 1161, 1, 646, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1216, 1213, 1361, 2, 842, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1249, 1247, 1376, 2, 857, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1266, 1262, 1381, 2, 995, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1364, 1361, 1416, 2, 1030, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1385, 1381, 1422, 2, 1100, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1418, 1416, 1437, 2, 1115, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1425, 1422, 1438, 2, 1145, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1440, 1438, 1443, 2, 1160, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1446, 1444, 1478, 3, 1361, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1475, 1474, 1492, 3, 1376, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1480, 1478, 1493, 3, 1416, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1495, 1493, 1498, 3, 1437, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 418, 415, 430, 541, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 509, 505, 540, 591, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 546, 541, 555, 616, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 597, 591, 615, 634, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 639, 634, 642, 647, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 811, 808, 842, 996, 1, 415, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 848, 843, 857, 1031, 1, 431, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 882, 878, 995, 1046, 1, 505, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1002, 996, 1030, 1101, 1, 541, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1053, 1046, 1100, 1121, 1, 591, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1106, 1101, 1115, 1146, 1, 616, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1127, 1121, 1145, 1152, 1, 634, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1157, 1152, 1160, 1162, 1, 647, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1219, 1214, 1247, 1362, 2, 843, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1269, 1263, 1361, 1382, 2, 996, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1367, 1362, 1376, 1417, 2, 1031, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1388, 1382, 1416, 1423, 2, 1101, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1428, 1423, 1437, 1439, 2, 1146, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1449, 1445, 1474, 1479, 3, 1362, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 1483, 1479, 1492, 1494, 3, 1417, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 384, 381, 415, 506, 540, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 421, 416, 431, 542, 555, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 458, 455, 505, 566, 590, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 512, 506, 541, 592, 615, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 569, 566, 591, 624, 633, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 600, 592, 616, 635, 642, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 627, 624, 634, 644, 646, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 716, 713, 808, 879, 995, 1, 381, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 814, 809, 843, 997, 1030, 1, 416, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 885, 879, 996, 1047, 1100, 1, 506, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1005, 997, 1031, 1102, 1115, 1, 542, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1056, 1047, 1101, 1122, 1145, 1, 592, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1130, 1122, 1146, 1153, 1160, 1, 635, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1185, 1183, 1214, 1264, 1361, 2, 809, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1222, 1215, 1248, 1363, 1376, 2, 844, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1272, 1264, 1362, 1383, 1416, 2, 997, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 1391, 1383, 1417, 1424, 1437, 2, 1102, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 280, 278, 381, 456, 505, 995, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 387, 382, 416, 507, 541, 1030, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 461, 456, 506, 567, 591, 1100, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 515, 507, 542, 593, 616, 1115, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 572, 567, 592, 625, 634, 1145, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 630, 625, 635, 645, 647, 1160, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 719, 714, 809, 880, 996, 1361, 1, 382, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 817, 810, 844, 998, 1031, 1376, 1, 417, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 888, 880, 997, 1048, 1101, 1416, 1, 507, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1059, 1048, 1102, 1123, 1146, 1437, 1, 593, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1188, 1184, 1215, 1265, 1362, 1474, 2, 810, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1275, 1265, 1363, 1384, 1417, 1492, 2, 998, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1, 0, 278, 435, 455, 878, 1045, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 283, 279, 382, 457, 506, 996, 1100, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 390, 383, 417, 508, 542, 1031, 1115, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 436, 435, 456, 556, 566, 1046, 1120, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 464, 457, 507, 568, 592, 1101, 1145, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 557, 556, 567, 620, 624, 1121, 1151, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 575, 568, 593, 626, 635, 1146, 1160, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 621, 620, 625, 643, 644, 1152, 1161, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 649, 648, 714, 858, 879, 1263, 1381, 1, 279, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 722, 715, 810, 881, 997, 1362, 1416, 1, 383, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 859, 858, 880, 1035, 1047, 1382, 1422, 1, 457, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 891, 881, 998, 1049, 1102, 1417, 1437, 1, 508, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1036, 1035, 1048, 1116, 1122, 1423, 1438, 1, 568, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1117, 1116, 1123, 1150, 1153, 1439, 1443, 1, 626, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1164, 1163, 1184, 1252, 1264, 1445, 1478, 2, 715, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1253, 1252, 1265, 1377, 1383, 1479, 1493, 2, 881, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_g(pbuffer, 1378, 1377, 1384, 1421, 1424, 1494, 1498, 2, 1049, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 549, 540, 543, 615, 617, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 603, 590, 594, 633, 636, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 851, 842, 845, 1030, 1032, 1, 430, 432, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1008, 995, 999, 1100, 1103, 1, 540, 543, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1062, 1045, 1050, 1120, 1124, 1, 590, 594, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1109, 1100, 1103, 1145, 1147, 1, 615, 617, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1133, 1120, 1124, 1151, 1154, 1, 633, 636, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1225, 1213, 1216, 1361, 1364, 2, 842, 845, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1278, 1262, 1266, 1381, 1385, 2, 995, 999, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1370, 1361, 1364, 1416, 1418, 2, 1030, 1032, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1394, 1381, 1385, 1422, 1425, 2, 1100, 1103, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1431, 1422, 1425, 1438, 1440, 2, 1145, 1147, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1452, 1444, 1446, 1478, 1480, 3, 1361, 1364, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 1486, 1478, 1480, 1493, 1495, 3, 1416, 1418, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 424, 415, 418, 432, 541, 546, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 518, 505, 509, 543, 591, 597, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 609, 591, 597, 617, 634, 639, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 820, 808, 811, 845, 996, 1002, 1, 415, 418, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 894, 878, 882, 999, 1046, 1053, 1, 505, 509, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1014, 996, 1002, 1032, 1101, 1106, 1, 541, 546, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1068, 1046, 1053, 1103, 1121, 1127, 1, 591, 597, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1139, 1121, 1127, 1147, 1152, 1157, 1, 634, 639, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1231, 1214, 1219, 1249, 1362, 1367, 2, 843, 848, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1284, 1263, 1269, 1364, 1382, 1388, 2, 996, 1002, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1400, 1382, 1388, 1418, 1423, 1428, 2, 1101, 1106, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 1458, 1445, 1449, 1475, 1479, 1483, 3, 1362, 1367, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 393, 381, 384, 418, 506, 512, 540, 543, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 467, 455, 458, 509, 566, 569, 590, 594, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 524, 506, 512, 546, 592, 600, 615, 617, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 578, 566, 569, 597, 624, 627, 633, 636, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 725, 713, 716, 811, 879, 885, 995, 999, 1, 381, 384, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 826, 809, 814, 848, 997, 1005, 1030, 1032, 1, 416, 421, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 900, 879, 885, 1002, 1047, 1056, 1100, 1103, 1, 506, 512, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1074, 1047, 1056, 1106, 1122, 1130, 1145, 1147, 1, 592, 600, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1191, 1183, 1185, 1219, 1264, 1272, 1361, 1364, 2, 809, 814, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 1290, 1264, 1272, 1367, 1383, 1391, 1416, 1418, 2, 997, 1005, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 286, 278, 280, 384, 456, 461, 505, 509, 999, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 399, 382, 387, 421, 507, 515, 541, 546, 1032, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 473, 456, 461, 512, 567, 572, 591, 597, 1103, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 584, 567, 572, 600, 625, 630, 634, 639, 1147, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 731, 714, 719, 814, 880, 888, 996, 1002, 1364, 1, 382, 387, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 906, 880, 888, 1005, 1048, 1059, 1101, 1106, 1418, 1, 507, 515, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 1197, 1184, 1188, 1222, 1265, 1275, 1362, 1367, 1475, 2, 810, 817, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 4, 0, 1, 280, 435, 436, 455, 458, 882, 1045, 1050, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 292, 279, 283, 387, 457, 464, 506, 512, 1002, 1100, 1103, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 439, 435, 436, 461, 556, 557, 566, 569, 1053, 1120, 1124, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 479, 457, 464, 515, 568, 575, 592, 600, 1106, 1145, 1147, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 560, 556, 557, 572, 620, 621, 624, 627, 1127, 1151, 1154, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 652, 648, 649, 719, 858, 859, 879, 885, 1269, 1381, 1385, 1, 279, 283, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 737, 715, 722, 817, 881, 891, 997, 1005, 1367, 1416, 1418, 1, 383, 390, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 862, 858, 859, 888, 1035, 1036, 1047, 1056, 1388, 1422, 1425, 1, 457, 464, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 1039, 1035, 1036, 1059, 1116, 1117, 1122, 1130, 1428, 1438, 1440, 1, 568, 575, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 1167, 1163, 1164, 1188, 1252, 1253, 1264, 1272, 1449, 1478, 1480, 2, 715, 722, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_g(pbuffer, 1256, 1252, 1253, 1275, 1377, 1378, 1383, 1391, 1483, 1493, 1495, 2, 881, 891, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1020, 999, 1008, 1103, 1109, 1, 543, 549, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1080, 1050, 1062, 1124, 1133, 1, 594, 603, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1237, 1216, 1225, 1364, 1370, 2, 845, 851, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1296, 1266, 1278, 1385, 1394, 2, 999, 1008, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1406, 1385, 1394, 1425, 1431, 2, 1103, 1109, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 1464, 1446, 1452, 1480, 1486, 3, 1364, 1370, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 530, 509, 518, 549, 597, 609, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 832, 811, 820, 851, 1002, 1014, 1, 418, 424, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 912, 882, 894, 1008, 1053, 1068, 1, 509, 518, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1090, 1053, 1068, 1109, 1127, 1139, 1, 597, 609, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 1306, 1269, 1284, 1370, 1388, 1400, 2, 1002, 1014, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 405, 384, 393, 424, 512, 524, 543, 549, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 485, 458, 467, 518, 569, 578, 594, 603, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 743, 716, 725, 820, 885, 900, 999, 1008, 1, 384, 393, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 922, 885, 900, 1014, 1056, 1074, 1103, 1109, 1, 512, 524, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 1203, 1185, 1191, 1231, 1272, 1290, 1364, 1370, 2, 814, 826, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 298, 280, 286, 393, 461, 473, 509, 518, 1008, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 495, 461, 473, 524, 572, 584, 597, 609, 1109, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 753, 719, 731, 826, 888, 906, 1002, 1014, 1370, 1, 387, 399, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 10, 1, 4, 286, 436, 439, 458, 467, 894, 1050, 1062, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 308, 283, 292, 399, 464, 479, 512, 524, 1014, 1103, 1109, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 445, 436, 439, 473, 557, 560, 569, 578, 1068, 1124, 1133, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 658, 649, 652, 731, 859, 862, 885, 900, 1284, 1385, 1394, 1, 283, 292, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 868, 859, 862, 906, 1036, 1039, 1056, 1074, 1400, 1425, 1431, 1, 464, 479, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_g(pbuffer, 1173, 1164, 1167, 1197, 1253, 1256, 1272, 1290, 1458, 1480, 1486, 2, 722, 737, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 1316, 1278, 1296, 1394, 1406, 2, 1008, 1020, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 932, 894, 912, 1020, 1068, 1090, 1, 518, 530, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 763, 725, 743, 832, 900, 922, 1008, 1020, 1, 393, 405, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 318, 286, 298, 405, 473, 495, 518, 530, 1020, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 20, 4, 10, 298, 439, 445, 467, 485, 912, 1062, 1080, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_g(pbuffer, 668, 652, 658, 753, 862, 868, 900, 922, 1306, 1394, 1406, 1, 292, 308, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 947, 894, 999, 1008, 1284, 1, 820, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_f(pbuffer, 333, 286, 384, 393, 731, 820, 999, 1008, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_g(pbuffer, 35, 4, 280, 286, 652, 725, 882, 894, 1278, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 1331, 1296, 1464, 1, 1237, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 965, 912, 1008, 1020, 1306, 1, 832, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 778, 743, 820, 832, 1203, 1237, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 351, 298, 393, 405, 753, 832, 1008, 1020, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 53, 10, 286, 298, 658, 743, 894, 912, 1296, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_g(pbuffer, 683, 658, 731, 753, 1173, 1203, 1284, 1306, 1464, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_g(pbuffer, 173, 20, 298, 318, 668, 763, 912, 932, 1316, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_g(pbuffer, 218, 10, 53, 333, 351, 658, 683, 743, 778, 947, 965, 1296, 1331, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2cgeom::comp_prim_op_geom_010_pf(pbuffer, 83, 10, 35, 173, 218, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, 0, pbuffer, 83, 90, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<1, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 1, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2pecp namespace

#endif /* ProjectedCorePotentialGeom010PFForG_hpp */
