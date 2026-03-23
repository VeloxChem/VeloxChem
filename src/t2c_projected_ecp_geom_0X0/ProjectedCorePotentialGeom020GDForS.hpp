#ifndef ProjectedCorePotentialGeom020GDForS_hpp
#define ProjectedCorePotentialGeom020GDForS_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "SimdArray.hpp"
#include "ProjectedCorePotentialPrimRecDDForS.hpp"
#include "ProjectedCorePotentialPrimRecDFForS.hpp"
#include "ProjectedCorePotentialPrimRecDGForS.hpp"
#include "ProjectedCorePotentialPrimRecDPForS.hpp"
#include "ProjectedCorePotentialPrimRecDSForS.hpp"
#include "ProjectedCorePotentialPrimRecFDForS.hpp"
#include "ProjectedCorePotentialPrimRecFFForS.hpp"
#include "ProjectedCorePotentialPrimRecFGForS.hpp"
#include "ProjectedCorePotentialPrimRecFPForS.hpp"
#include "ProjectedCorePotentialPrimRecFSForS.hpp"
#include "ProjectedCorePotentialPrimRecGDForS.hpp"
#include "ProjectedCorePotentialPrimRecGFForS.hpp"
#include "ProjectedCorePotentialPrimRecGGForS.hpp"
#include "ProjectedCorePotentialPrimRecGPForS.hpp"
#include "ProjectedCorePotentialPrimRecGSForS.hpp"
#include "ProjectedCorePotentialPrimRecHDForS.hpp"
#include "ProjectedCorePotentialPrimRecHFForS.hpp"
#include "ProjectedCorePotentialPrimRecHPForS.hpp"
#include "ProjectedCorePotentialPrimRecIDForS.hpp"
#include "ProjectedCorePotentialPrimRecPDForS.hpp"
#include "ProjectedCorePotentialPrimRecPFForS.hpp"
#include "ProjectedCorePotentialPrimRecPGForS.hpp"
#include "ProjectedCorePotentialPrimRecPPForS.hpp"
#include "ProjectedCorePotentialPrimRecPSForS.hpp"
#include "ProjectedCorePotentialPrimRecSDForS.hpp"
#include "ProjectedCorePotentialPrimRecSFForS.hpp"
#include "ProjectedCorePotentialPrimRecSGForS.hpp"
#include "ProjectedCorePotentialPrimRecSPForS.hpp"
#include "ProjectedCorePotentialPrimRecSS.hpp"
#include "GeometricalDerivatives020ForGD.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes (G|U_l|D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_projected_core_potential_geom_020_gd_for_s(T& distributor,
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

    CSimdArray<double> i_values(9, ket_npgtos);

    CSimdArray<double> l_values(1, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(4981, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(540, 1);

    CSimdArray<double> sbuffer(270, 1);

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

                    t2cfunc::comp_i_vals(i_values, 8, pfactors, 8);

                    t2cfunc::comp_l_vals(l_values, 0, pfactors, 8, 6);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 0, pbuffer, 0, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 1, pbuffer, 2332, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 2, pbuffer, 2352, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 3, pbuffer, 2362, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 4, pbuffer, 2366, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 0, pbuffer, 2367, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 1, pbuffer, 2387, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 2, pbuffer, 2397, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 3, pbuffer, 2401, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 0, pbuffer, 2402, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 1, pbuffer, 2412, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 2, pbuffer, 2416, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 0, pbuffer, 2417, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 1, pbuffer, 2421, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 0, pbuffer, 2422, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 0, pbuffer, 2423, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 1, pbuffer, 3534, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 2, pbuffer, 3554, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 3, pbuffer, 3564, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 0, pbuffer, 3568, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 1, pbuffer, 3588, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 2, pbuffer, 3598, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 3, pbuffer, 3602, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 0, pbuffer, 3603, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 1, pbuffer, 3613, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 2, pbuffer, 3617, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 0, pbuffer, 3618, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 1, pbuffer, 3622, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 0, pbuffer, 3623, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 0, pbuffer, 3624, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 1, pbuffer, 4254, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 2, pbuffer, 4274, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 3, pbuffer, 4284, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 0, pbuffer, 4285, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 1, pbuffer, 4305, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 2, pbuffer, 4315, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 0, pbuffer, 4319, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 1, pbuffer, 4329, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 2, pbuffer, 4333, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 0, pbuffer, 4334, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 1, pbuffer, 4338, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 0, pbuffer, 4339, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 0, pbuffer, 4340, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 1, pbuffer, 4654, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 2, pbuffer, 4674, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 0, pbuffer, 4678, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 1, pbuffer, 4698, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 2, pbuffer, 4708, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 0, pbuffer, 4709, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 1, pbuffer, 4719, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 0, pbuffer, 4723, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 1, pbuffer, 4727, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 0, pbuffer, 4728, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 0, pbuffer, 4729, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 1, pbuffer, 4857, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 2, pbuffer, 4867, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 0, pbuffer, 4868, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 1, pbuffer, 4888, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 0, pbuffer, 4892, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 1, pbuffer, 4902, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 0, pbuffer, 4903, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 4, 0, pbuffer, 4907, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 0, 0, pbuffer, 4908, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 0, 1, pbuffer, 4946, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 0, pbuffer, 4950, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 1, pbuffer, 4960, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 2, 0, pbuffer, 4961, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 3, 0, pbuffer, 4965, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 0, 0, pbuffer, 4966, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 1, 0, pbuffer, 4976, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 6, 2, 0, pbuffer, 4980, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 1, 0, 2367, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2333, 2332, 2387, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2353, 2352, 2397, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2363, 2362, 2401, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2368, 2367, 2402, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2388, 2387, 2412, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2398, 2397, 2416, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2403, 2402, 2417, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2413, 2412, 2421, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2418, 2417, 2422, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 2424, 2423, 3568, 1, 2332, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3535, 3534, 3588, 1, 2352, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3555, 3554, 3598, 1, 2362, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3565, 3564, 3602, 1, 2366, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3569, 3568, 3603, 1, 2387, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3589, 3588, 3613, 1, 2397, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3599, 3598, 3617, 1, 2401, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3604, 3603, 3618, 1, 2412, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3614, 3613, 3622, 1, 2416, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3619, 3618, 3623, 1, 2421, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 3625, 3624, 4285, 2, 3534, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4255, 4254, 4305, 2, 3554, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4275, 4274, 4315, 2, 3564, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4286, 4285, 4319, 2, 3588, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4306, 4305, 4329, 2, 3598, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4316, 4315, 4333, 2, 3602, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4320, 4319, 4334, 2, 3613, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4330, 4329, 4338, 2, 3617, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4335, 4334, 4339, 2, 3622, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4341, 4340, 4678, 3, 4254, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4655, 4654, 4698, 3, 4274, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4675, 4674, 4708, 3, 4284, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4679, 4678, 4709, 3, 4305, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4699, 4698, 4719, 3, 4315, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4710, 4709, 4723, 3, 4329, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4720, 4719, 4727, 3, 4333, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4724, 4723, 4728, 3, 4338, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4730, 4729, 4868, 4, 4654, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4858, 4857, 4888, 4, 4674, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4869, 4868, 4892, 4, 4698, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4889, 4888, 4902, 4, 4708, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4893, 4892, 4903, 4, 4719, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4904, 4903, 4907, 4, 4727, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4909, 4908, 4950, 5, 4857, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4947, 4946, 4960, 5, 4867, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4951, 4950, 4961, 5, 4888, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4962, 4961, 4965, 5, 4902, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4967, 4966, 4976, 6, 4946, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4977, 4976, 4980, 6, 4960, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4, 0, 1, 2367, 2368, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2336, 2332, 2333, 2387, 2388, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2356, 2352, 2353, 2397, 2398, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2371, 2367, 2368, 2402, 2403, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2391, 2387, 2388, 2412, 2413, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2406, 2402, 2403, 2417, 2418, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 2427, 2423, 2424, 3568, 3569, 1, 2332, 2333, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3538, 3534, 3535, 3588, 3589, 1, 2352, 2353, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3558, 3554, 3555, 3598, 3599, 1, 2362, 2363, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3572, 3568, 3569, 3603, 3604, 1, 2387, 2388, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3592, 3588, 3589, 3613, 3614, 1, 2397, 2398, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3607, 3603, 3604, 3618, 3619, 1, 2412, 2413, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 3628, 3624, 3625, 4285, 4286, 2, 3534, 3535, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4258, 4254, 4255, 4305, 4306, 2, 3554, 3555, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4278, 4274, 4275, 4315, 4316, 2, 3564, 3565, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4289, 4285, 4286, 4319, 4320, 2, 3588, 3589, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4309, 4305, 4306, 4329, 4330, 2, 3598, 3599, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4323, 4319, 4320, 4334, 4335, 2, 3613, 3614, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4344, 4340, 4341, 4678, 4679, 3, 4254, 4255, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4658, 4654, 4655, 4698, 4699, 3, 4274, 4275, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4682, 4678, 4679, 4709, 4710, 3, 4305, 4306, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4702, 4698, 4699, 4719, 4720, 3, 4315, 4316, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4713, 4709, 4710, 4723, 4724, 3, 4329, 4330, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4733, 4729, 4730, 4868, 4869, 4, 4654, 4655, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4861, 4857, 4858, 4888, 4889, 4, 4674, 4675, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4872, 4868, 4869, 4892, 4893, 4, 4698, 4699, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4896, 4892, 4893, 4903, 4904, 4, 4719, 4720, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4912, 4908, 4909, 4950, 4951, 5, 4857, 4858, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4954, 4950, 4951, 4961, 4962, 5, 4888, 4889, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4970, 4966, 4967, 4976, 4977, 6, 4946, 4947, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 10, 1, 4, 2368, 2371, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2342, 2333, 2336, 2388, 2391, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2377, 2368, 2371, 2403, 2406, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 2433, 2424, 2427, 3569, 3572, 1, 2333, 2336, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3544, 3535, 3538, 3589, 3592, 1, 2353, 2356, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3578, 3569, 3572, 3604, 3607, 1, 2388, 2391, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 3634, 3625, 3628, 4286, 4289, 2, 3535, 3538, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4264, 4255, 4258, 4306, 4309, 2, 3555, 3558, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4295, 4286, 4289, 4320, 4323, 2, 3589, 3592, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4350, 4341, 4344, 4679, 4682, 3, 4255, 4258, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4664, 4655, 4658, 4699, 4702, 3, 4275, 4278, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4688, 4679, 4682, 4710, 4713, 3, 4306, 4309, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4739, 4730, 4733, 4869, 4872, 4, 4655, 4658, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4878, 4869, 4872, 4893, 4896, 4, 4699, 4702, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4918, 4909, 4912, 4951, 4954, 5, 4858, 4861, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 20, 4, 10, 2371, 2377, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 2443, 2427, 2433, 3572, 3578, 1, 2336, 2342, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 3644, 3628, 3634, 4289, 4295, 2, 3538, 3544, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4360, 4344, 4350, 4682, 4688, 3, 4258, 4264, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4749, 4733, 4739, 4872, 4878, 4, 4658, 4664, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_s(pbuffer, 35, 0, 2423, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_s(pbuffer, 2458, 2423, 3624, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_s(pbuffer, 3659, 3624, 4340, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_s(pbuffer, 4375, 4340, 4729, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 38, 1, 2424, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 2461, 2424, 3625, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 3662, 3625, 4341, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 4378, 4341, 4730, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 4764, 4730, 4909, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 47, 4, 2427, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 2470, 2427, 3628, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 3671, 3628, 4344, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 4387, 4344, 4733, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 4773, 4733, 4912, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 4928, 4912, 4970, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 65, 10, 2433, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 2488, 2433, 3634, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 3689, 3634, 4350, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 4405, 4350, 4739, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 4791, 4739, 4918, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 95, 20, 2443, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 2518, 2443, 3644, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 3719, 3644, 4360, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 4435, 4360, 4749, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ds_s(pbuffer, 140, 0, 35, 2423, 2458, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ds_s(pbuffer, 2563, 2423, 2458, 3624, 3659, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ds_s(pbuffer, 3764, 3624, 3659, 4340, 4375, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_s(pbuffer, 146, 1, 38, 2424, 2461, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_s(pbuffer, 2569, 2424, 2461, 3625, 3662, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_s(pbuffer, 3770, 3625, 3662, 4341, 4378, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_s(pbuffer, 4480, 4341, 4378, 4730, 4764, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 164, 4, 47, 2427, 2470, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 2587, 2427, 2470, 3628, 3671, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 3788, 3628, 3671, 4344, 4387, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 4498, 4344, 4387, 4733, 4773, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 4821, 4733, 4773, 4912, 4928, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 200, 10, 65, 2433, 2488, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 2623, 2433, 2488, 3634, 3689, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 3824, 3634, 3689, 4350, 4405, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 4534, 4350, 4405, 4739, 4791, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 260, 20, 95, 2443, 2518, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 2683, 2443, 2518, 3644, 3719, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 3884, 3644, 3719, 4360, 4435, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fs_s(pbuffer, 350, 35, 140, 2458, 2563, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fs_s(pbuffer, 2773, 2458, 2563, 3659, 3764, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fp_s(pbuffer, 360, 38, 146, 2461, 2569, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fp_s(pbuffer, 2783, 2461, 2569, 3662, 3770, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fp_s(pbuffer, 3974, 3662, 3770, 4378, 4480, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_s(pbuffer, 390, 47, 164, 2470, 2587, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_s(pbuffer, 2813, 2470, 2587, 3671, 3788, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_s(pbuffer, 4004, 3671, 3788, 4387, 4498, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_s(pbuffer, 4594, 4387, 4498, 4773, 4821, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_s(pbuffer, 450, 65, 200, 2488, 2623, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_s(pbuffer, 2873, 2488, 2623, 3689, 3824, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_s(pbuffer, 4064, 3689, 3824, 4405, 4534, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_s(pbuffer, 550, 95, 260, 2518, 2683, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_s(pbuffer, 2973, 2518, 2683, 3719, 3884, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gs_s(pbuffer, 700, 140, 350, 2563, 2773, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gp_s(pbuffer, 715, 146, 360, 2569, 2783, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gp_s(pbuffer, 3123, 2569, 2783, 3770, 3974, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gd_s(pbuffer, 760, 164, 390, 2587, 2813, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gd_s(pbuffer, 3168, 2587, 2813, 3788, 4004, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gd_s(pbuffer, 4164, 3788, 4004, 4498, 4594, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gf_s(pbuffer, 1390, 200, 450, 2623, 2873, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gf_s(pbuffer, 3258, 2623, 2873, 3824, 4064, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gg_s(pbuffer, 1540, 260, 550, 2683, 2973, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_hp_s(pbuffer, 1765, 360, 715, 2783, 3123, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_hd_s(pbuffer, 1828, 390, 760, 2813, 3168, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_hd_s(pbuffer, 3408, 2813, 3168, 4004, 4164, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_hf_s(pbuffer, 1954, 450, 1390, 2873, 3258, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_id_s(pbuffer, 2164, 760, 1828, 3168, 3408, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2cgeom::comp_prim_op_geom_020_gd(pbuffer, 850, 164, 360, 450, 700, 760, 1540, 1765, 1954, 2164, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, 0, pbuffer, 850, 540, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2pecp namespace

#endif /* ProjectedCorePotentialGeom020GDForS_hpp */
