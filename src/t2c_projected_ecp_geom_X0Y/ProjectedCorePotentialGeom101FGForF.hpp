#ifndef ProjectedCorePotentialGeom101FGForF_hpp
#define ProjectedCorePotentialGeom101FGForF_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "SimdArray.hpp"
#include "ProjectedCorePotentialPrimRecDDForD.hpp"
#include "ProjectedCorePotentialPrimRecDDForP.hpp"
#include "ProjectedCorePotentialPrimRecDDForS.hpp"
#include "ProjectedCorePotentialPrimRecDFForD.hpp"
#include "ProjectedCorePotentialPrimRecDFForF.hpp"
#include "ProjectedCorePotentialPrimRecDFForP.hpp"
#include "ProjectedCorePotentialPrimRecDFForS.hpp"
#include "ProjectedCorePotentialPrimRecDGForD.hpp"
#include "ProjectedCorePotentialPrimRecDGForP.hpp"
#include "ProjectedCorePotentialPrimRecDGForS.hpp"
#include "ProjectedCorePotentialPrimRecDHForD.hpp"
#include "ProjectedCorePotentialPrimRecDHForF.hpp"
#include "ProjectedCorePotentialPrimRecDHForP.hpp"
#include "ProjectedCorePotentialPrimRecDHForS.hpp"
#include "ProjectedCorePotentialPrimRecDPForP.hpp"
#include "ProjectedCorePotentialPrimRecFDForD.hpp"
#include "ProjectedCorePotentialPrimRecFDForS.hpp"
#include "ProjectedCorePotentialPrimRecFFForD.hpp"
#include "ProjectedCorePotentialPrimRecFFForF.hpp"
#include "ProjectedCorePotentialPrimRecFFForP.hpp"
#include "ProjectedCorePotentialPrimRecFFForS.hpp"
#include "ProjectedCorePotentialPrimRecFGForD.hpp"
#include "ProjectedCorePotentialPrimRecFGForS.hpp"
#include "ProjectedCorePotentialPrimRecFHForD.hpp"
#include "ProjectedCorePotentialPrimRecFHForF.hpp"
#include "ProjectedCorePotentialPrimRecFHForP.hpp"
#include "ProjectedCorePotentialPrimRecFHForS.hpp"
#include "ProjectedCorePotentialPrimRecGFForF.hpp"
#include "ProjectedCorePotentialPrimRecGHForF.hpp"
#include "ProjectedCorePotentialPrimRecPDForD.hpp"
#include "ProjectedCorePotentialPrimRecPDForP.hpp"
#include "ProjectedCorePotentialPrimRecPDForS.hpp"
#include "ProjectedCorePotentialPrimRecPFForD.hpp"
#include "ProjectedCorePotentialPrimRecPFForF.hpp"
#include "ProjectedCorePotentialPrimRecPFForP.hpp"
#include "ProjectedCorePotentialPrimRecPFForS.hpp"
#include "ProjectedCorePotentialPrimRecPGForD.hpp"
#include "ProjectedCorePotentialPrimRecPGForP.hpp"
#include "ProjectedCorePotentialPrimRecPGForS.hpp"
#include "ProjectedCorePotentialPrimRecPHForD.hpp"
#include "ProjectedCorePotentialPrimRecPHForF.hpp"
#include "ProjectedCorePotentialPrimRecPHForP.hpp"
#include "ProjectedCorePotentialPrimRecPHForS.hpp"
#include "ProjectedCorePotentialPrimRecPPForP.hpp"
#include "ProjectedCorePotentialPrimRecPPForS.hpp"
#include "ProjectedCorePotentialPrimRecPSForS.hpp"
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
#include "ProjectedCorePotentialPrimRecSPForD.hpp"
#include "ProjectedCorePotentialPrimRecSPForF.hpp"
#include "ProjectedCorePotentialPrimRecSPForP.hpp"
#include "ProjectedCorePotentialPrimRecSPForS.hpp"
#include "ProjectedCorePotentialPrimRecSS.hpp"
#include "GeometricalDerivatives1X1ForFG.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes (d^(1)/dA^(1)F|U_l|d^(1)/dB^(1)G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_projected_core_potential_geom_101_fg_for_f(T& distributor,
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

    CSimdArray<double> i_values(13, ket_npgtos);

    CSimdArray<double> l_values(4, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(13896, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(1350, 1);

    CSimdArray<double> sbuffer(567, 1);

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

                    t2cfunc::comp_i_vals(i_values, 12, pfactors, 8);

                    t2cfunc::comp_l_vals(l_values, 3, pfactors, 8, 6);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 0, pbuffer, 0, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 1, pbuffer, 2460, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 1, pbuffer, 2461, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 2, pbuffer, 3539, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 2, pbuffer, 3540, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 2, pbuffer, 3541, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 3, pbuffer, 4145, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 3, pbuffer, 4146, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 3, pbuffer, 4147, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 3, pbuffer, 4148, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 4, pbuffer, 4434, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 4, pbuffer, 4435, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 4, pbuffer, 4436, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 0, 4, pbuffer, 4437, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 5, pbuffer, 4503, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 0, 5, pbuffer, 4504, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 0, 5, pbuffer, 4505, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 0, 6, pbuffer, 4518, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 0, pbuffer, 4519, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 1, pbuffer, 4554, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 1, pbuffer, 4555, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 1, pbuffer, 4556, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 2, pbuffer, 4644, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 2, pbuffer, 4645, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 2, pbuffer, 4646, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 2, pbuffer, 4647, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 3, pbuffer, 4744, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 3, pbuffer, 4745, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 3, pbuffer, 4746, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 3, pbuffer, 4747, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 4, pbuffer, 4813, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 4, pbuffer, 4814, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 1, 4, pbuffer, 4815, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 1, 4, pbuffer, 4816, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 1, 5, pbuffer, 4848, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 1, 5, pbuffer, 4849, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 0, pbuffer, 4853, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 1, pbuffer, 4873, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 1, pbuffer, 4874, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 1, pbuffer, 4875, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 2, pbuffer, 4923, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 2, pbuffer, 4924, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 2, pbuffer, 4925, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 2, pbuffer, 4926, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 3, pbuffer, 4977, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 3, pbuffer, 4978, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 3, pbuffer, 4979, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 2, 3, pbuffer, 4980, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 4, pbuffer, 5012, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 2, 4, pbuffer, 5013, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 2, 4, pbuffer, 5014, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 2, 5, pbuffer, 5027, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 0, pbuffer, 5028, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 1, pbuffer, 5038, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 1, pbuffer, 5039, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 1, pbuffer, 5040, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 2, pbuffer, 5062, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 2, pbuffer, 5063, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 2, pbuffer, 5064, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 3, 2, pbuffer, 5065, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 3, pbuffer, 5087, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 3, pbuffer, 5088, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 3, 3, pbuffer, 5089, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 3, 4, pbuffer, 5102, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 3, 4, pbuffer, 5103, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 0, pbuffer, 5107, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 1, pbuffer, 5111, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 1, pbuffer, 5112, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 4, 1, pbuffer, 5113, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 2, pbuffer, 5120, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 2, pbuffer, 5121, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 4, 2, pbuffer, 5122, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 3, pbuffer, 5129, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 4, 3, pbuffer, 5130, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 4, 4, pbuffer, 5134, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 0, 5, 0, pbuffer, 5135, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 5, 1, pbuffer, 5136, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 0, 5, 1, pbuffer, 5137, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 5, 2, pbuffer, 5138, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 0, 5, 2, pbuffer, 5139, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 0, 5, 3, pbuffer, 5140, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 0, pbuffer, 5141, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 1, pbuffer, 5786, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 1, pbuffer, 5787, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 1, pbuffer, 5788, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 2, pbuffer, 6990, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 2, pbuffer, 6991, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 2, pbuffer, 6992, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 2, pbuffer, 6993, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 3, pbuffer, 7790, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 3, pbuffer, 7791, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 3, pbuffer, 7792, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 0, 3, pbuffer, 7793, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 4, pbuffer, 7911, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 4, pbuffer, 7912, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 0, 4, pbuffer, 7913, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 0, 5, pbuffer, 7945, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 0, 5, pbuffer, 7946, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 0, pbuffer, 7950, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 1, pbuffer, 7985, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 1, pbuffer, 7986, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 1, pbuffer, 7987, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 1, pbuffer, 7988, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 2, pbuffer, 9119, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 2, pbuffer, 9120, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 2, pbuffer, 9121, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 2, pbuffer, 9122, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 3, pbuffer, 9219, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 3, pbuffer, 9220, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 3, pbuffer, 9221, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 1, 3, pbuffer, 9222, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 4, pbuffer, 9288, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 1, 4, pbuffer, 9289, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 1, 4, pbuffer, 9290, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 1, 5, pbuffer, 9303, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 0, pbuffer, 9304, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 1, pbuffer, 9324, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 1, pbuffer, 9325, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 1, pbuffer, 9326, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 1, pbuffer, 9327, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 2, pbuffer, 9409, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 2, pbuffer, 9410, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 2, pbuffer, 9411, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 2, pbuffer, 9412, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 3, pbuffer, 9463, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 3, pbuffer, 9464, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 2, 3, pbuffer, 9465, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 2, 3, pbuffer, 9466, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 2, 4, pbuffer, 9498, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 2, 4, pbuffer, 9499, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 0, pbuffer, 9503, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 1, pbuffer, 9513, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 1, pbuffer, 9514, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 1, pbuffer, 9515, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 1, pbuffer, 9516, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 2, pbuffer, 9557, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 2, pbuffer, 9558, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 2, pbuffer, 9559, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 3, 2, pbuffer, 9560, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 3, pbuffer, 9582, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 3, 3, pbuffer, 9583, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 3, 3, pbuffer, 9584, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 3, 4, pbuffer, 9597, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 4, 0, pbuffer, 9598, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 1, pbuffer, 9602, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 1, pbuffer, 9603, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 1, pbuffer, 9604, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 4, 1, pbuffer, 9605, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 2, pbuffer, 9621, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 2, pbuffer, 9622, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 4, 2, pbuffer, 9623, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 4, 3, pbuffer, 9630, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 4, 3, pbuffer, 9631, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 1, 5, 0, pbuffer, 9635, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 1, pbuffer, 9636, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 1, pbuffer, 9637, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 1, 5, 1, pbuffer, 9638, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 2, pbuffer, 9642, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 1, 5, 2, pbuffer, 9643, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 5, 3, pbuffer, 9644, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 1, 6, 1, pbuffer, 9645, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 0, pbuffer, 9646, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 1, pbuffer, 9981, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 1, pbuffer, 9982, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 1, pbuffer, 9983, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 2, pbuffer, 10563, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 2, pbuffer, 10564, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 2, pbuffer, 10565, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 0, 2, pbuffer, 10566, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 3, pbuffer, 10886, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 3, pbuffer, 10887, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 0, 3, pbuffer, 10888, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 4, pbuffer, 10951, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 0, 4, pbuffer, 10952, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 0, 5, pbuffer, 10965, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 0, pbuffer, 10966, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 1, pbuffer, 11001, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 1, pbuffer, 11002, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 1, pbuffer, 11003, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 1, pbuffer, 11004, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 2, pbuffer, 11615, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 2, pbuffer, 11616, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 2, pbuffer, 11617, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 1, 2, pbuffer, 11618, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 3, pbuffer, 11715, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 3, pbuffer, 11716, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 1, 3, pbuffer, 11717, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 1, 4, pbuffer, 11749, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 1, 4, pbuffer, 11750, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 0, pbuffer, 11754, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 1, pbuffer, 11774, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 1, pbuffer, 11775, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 1, pbuffer, 11776, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 1, pbuffer, 11777, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 2, pbuffer, 11859, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 2, pbuffer, 11860, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 2, pbuffer, 11861, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 2, 2, pbuffer, 11862, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 3, pbuffer, 11913, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 2, 3, pbuffer, 11914, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 2, 3, pbuffer, 11915, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 2, 4, pbuffer, 11928, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 0, pbuffer, 11929, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 1, pbuffer, 11939, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 1, pbuffer, 11940, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 1, pbuffer, 11941, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 1, pbuffer, 11942, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 2, pbuffer, 11983, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 2, pbuffer, 11984, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 3, 2, pbuffer, 11985, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 3, 2, pbuffer, 11986, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 3, 3, pbuffer, 12008, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 3, 3, pbuffer, 12009, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 4, 0, pbuffer, 12013, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 1, pbuffer, 12017, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 1, pbuffer, 12018, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 4, 1, pbuffer, 12019, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 4, 1, pbuffer, 12020, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 2, pbuffer, 12036, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 4, 2, pbuffer, 12037, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 4, 2, pbuffer, 12038, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 4, 3, pbuffer, 12045, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 2, 5, 0, pbuffer, 12046, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 1, pbuffer, 12047, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 5, 1, pbuffer, 12048, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 2, 5, 1, pbuffer, 12049, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 5, 2, pbuffer, 12053, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 2, 5, 2, pbuffer, 12054, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 2, 6, 1, pbuffer, 12055, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 0, pbuffer, 12056, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 1, pbuffer, 12205, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 1, pbuffer, 12206, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 0, 1, pbuffer, 12207, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 2, pbuffer, 12445, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 2, pbuffer, 12446, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 0, 2, pbuffer, 12447, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 3, pbuffer, 12556, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 0, 3, pbuffer, 12557, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 0, 4, pbuffer, 12586, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 0, pbuffer, 12590, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 1, pbuffer, 12625, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 1, pbuffer, 12626, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 1, pbuffer, 12627, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 1, 1, pbuffer, 12628, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 2, pbuffer, 12927, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 2, pbuffer, 12928, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 1, 2, pbuffer, 12929, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 3, pbuffer, 12992, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 1, 3, pbuffer, 12993, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 1, 4, pbuffer, 13006, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 2, 0, pbuffer, 13007, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 1, pbuffer, 13027, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 1, pbuffer, 13028, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 1, pbuffer, 13029, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 2, 1, pbuffer, 13030, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 2, pbuffer, 13112, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 2, pbuffer, 13113, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 2, 2, pbuffer, 13114, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 2, 3, pbuffer, 13146, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 2, 3, pbuffer, 13147, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 3, 0, pbuffer, 13151, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 1, pbuffer, 13161, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 1, pbuffer, 13162, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 3, 1, pbuffer, 13163, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 3, 1, pbuffer, 13164, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 2, pbuffer, 13205, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 3, 2, pbuffer, 13206, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 3, 2, pbuffer, 13207, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 3, 3, pbuffer, 13220, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 4, 0, pbuffer, 13221, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 1, pbuffer, 13225, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 4, 1, pbuffer, 13226, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 4, 1, pbuffer, 13227, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 4, 1, pbuffer, 13228, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 4, 2, pbuffer, 13244, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 4, 2, pbuffer, 13245, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 3, 5, 0, pbuffer, 13249, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 5, 1, pbuffer, 13250, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 3, 5, 1, pbuffer, 13251, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 3, 5, 1, pbuffer, 13252, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 5, 2, pbuffer, 13256, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 3, 6, 1, pbuffer, 13257, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 0, 0, pbuffer, 13258, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 0, 1, pbuffer, 13314, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 0, 1, pbuffer, 13315, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 2, pbuffer, 13405, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 0, 2, pbuffer, 13406, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 0, 3, pbuffer, 13460, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 1, 0, pbuffer, 13470, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 1, pbuffer, 13505, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 1, 1, pbuffer, 13506, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 1, 1, pbuffer, 13507, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 2, pbuffer, 13616, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 1, 2, pbuffer, 13617, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 1, 3, pbuffer, 13646, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 2, 0, pbuffer, 13650, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 1, pbuffer, 13670, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 2, 1, pbuffer, 13671, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 2, 1, pbuffer, 13672, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 2, pbuffer, 13735, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 2, 2, pbuffer, 13736, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 2, 3, pbuffer, 13749, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 3, 0, pbuffer, 13750, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 1, pbuffer, 13760, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 3, 1, pbuffer, 13761, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 3, 1, pbuffer, 13762, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 3, 2, pbuffer, 13794, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 3, 2, pbuffer, 13795, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 4, 0, pbuffer, 13799, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 4, 1, pbuffer, 13803, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 4, 1, pbuffer, 13804, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(2, 4, 4, 1, pbuffer, 13805, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 4, 2, pbuffer, 13818, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(3, 4, 5, 0, pbuffer, 13819, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 5, 1, pbuffer, 13820, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(1, 4, 5, 1, pbuffer, 13821, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 4, 6, 1, pbuffer, 13825, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 1, 1, pbuffer, 13826, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 2, 1, pbuffer, 13861, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 3, 1, pbuffer, 13881, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 4, 1, pbuffer, 13891, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_ss(0, 5, 5, 1, pbuffer, 13895, i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4149, 4145, 4744, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4438, 4434, 4813, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4506, 4503, 4848, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4648, 4644, 4923, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4748, 4744, 4977, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4817, 4813, 5012, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4850, 4848, 5027, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4927, 4923, 5062, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 4981, 4977, 5087, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5015, 5012, 5102, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5066, 5062, 5120, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5090, 5087, 5129, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5104, 5102, 5134, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5123, 5120, 5138, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 5131, 5129, 5140, 0, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 6994, 6990, 9119, 1, 4145, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7794, 7790, 9219, 1, 4434, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7914, 7911, 9288, 1, 4503, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7947, 7945, 9303, 1, 4518, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 7989, 7985, 9324, 1, 4644, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9123, 9119, 9409, 1, 4744, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9223, 9219, 9463, 1, 4813, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9291, 9288, 9498, 1, 4848, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9328, 9324, 9513, 1, 4923, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9413, 9409, 9557, 1, 4977, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9467, 9463, 9582, 1, 5012, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9500, 9498, 9597, 1, 5027, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9517, 9513, 9602, 1, 5062, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9561, 9557, 9621, 1, 5087, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9585, 9582, 9630, 1, 5102, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9606, 9602, 9636, 1, 5120, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9624, 9621, 9642, 1, 5129, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9632, 9630, 9644, 1, 5134, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 9639, 9636, 9645, 1, 5138, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 10567, 10563, 11615, 2, 7790, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 10889, 10886, 11715, 2, 7911, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 10953, 10951, 11749, 2, 7945, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11005, 11001, 11774, 2, 9119, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11619, 11615, 11859, 2, 9219, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11718, 11715, 11913, 2, 9288, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11751, 11749, 11928, 2, 9303, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11778, 11774, 11939, 2, 9409, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11863, 11859, 11983, 2, 9463, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11916, 11913, 12008, 2, 9498, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11943, 11939, 12017, 2, 9557, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 11987, 11983, 12036, 2, 9582, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12010, 12008, 12045, 2, 9597, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12021, 12017, 12047, 2, 9621, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12039, 12036, 12053, 2, 9630, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12050, 12047, 12055, 2, 9642, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12448, 12445, 12927, 3, 10886, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12558, 12556, 12992, 3, 10951, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12587, 12586, 13006, 3, 10965, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12629, 12625, 13027, 3, 11615, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12930, 12927, 13112, 3, 11715, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 12994, 12992, 13146, 3, 11749, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13031, 13027, 13161, 3, 11859, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13115, 13112, 13205, 3, 11913, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13148, 13146, 13220, 3, 11928, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13165, 13161, 13225, 3, 11983, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13208, 13205, 13244, 3, 12008, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13229, 13225, 13250, 3, 12036, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13246, 13244, 13256, 3, 12045, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13253, 13250, 13257, 3, 12053, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13407, 13405, 13616, 4, 12556, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13461, 13460, 13646, 4, 12586, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13508, 13505, 13670, 4, 12927, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13618, 13616, 13735, 4, 12992, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13647, 13646, 13749, 4, 13006, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13673, 13670, 13760, 4, 13112, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13737, 13735, 13794, 4, 13146, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13763, 13760, 13803, 4, 13205, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13796, 13794, 13818, 4, 13220, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13806, 13803, 13820, 4, 13244, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13822, 13820, 13825, 4, 13256, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13827, 13826, 13861, 5, 13616, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13862, 13861, 13881, 5, 13735, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13882, 13881, 13891, 5, 13794, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_s(pbuffer, 13892, 13891, 13895, 5, 13818, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 3542, 3539, 4145, 4645, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4152, 4146, 4434, 4745, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4441, 4435, 4503, 4814, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4509, 4504, 4518, 4849, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4557, 4554, 4644, 4873, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4651, 4645, 4744, 4924, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4751, 4745, 4813, 4978, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4820, 4814, 4848, 5013, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4876, 4873, 4923, 5038, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4930, 4924, 4977, 5063, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 4984, 4978, 5012, 5088, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5018, 5013, 5027, 5103, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5041, 5038, 5062, 5111, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5069, 5063, 5087, 5121, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5093, 5088, 5102, 5130, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5114, 5111, 5120, 5136, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5126, 5121, 5129, 5139, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 5789, 5786, 6990, 7986, 1, 3539, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 6997, 6991, 7790, 9120, 1, 4146, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7797, 7791, 7911, 9220, 1, 4435, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7917, 7912, 7945, 9289, 1, 4504, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 7992, 7986, 9119, 9325, 1, 4645, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9126, 9120, 9219, 9410, 1, 4745, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9226, 9220, 9288, 9464, 1, 4814, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9294, 9289, 9303, 9499, 1, 4849, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9331, 9325, 9409, 9514, 1, 4924, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9416, 9410, 9463, 9558, 1, 4978, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9470, 9464, 9498, 9583, 1, 5013, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9520, 9514, 9557, 9603, 1, 5063, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9564, 9558, 9582, 9622, 1, 5088, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9588, 9583, 9597, 9631, 1, 5103, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9609, 9603, 9621, 9637, 1, 5121, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9627, 9622, 9630, 9643, 1, 5130, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 9984, 9981, 10563, 11002, 2, 6991, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 10570, 10564, 10886, 11616, 2, 7791, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 10892, 10887, 10951, 11716, 2, 7912, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 10956, 10952, 10965, 11750, 2, 7946, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11008, 11002, 11615, 11775, 2, 9120, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11622, 11616, 11715, 11860, 2, 9220, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11721, 11716, 11749, 11914, 2, 9289, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11781, 11775, 11859, 11940, 2, 9410, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11866, 11860, 11913, 11984, 2, 9464, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11919, 11914, 11928, 12009, 2, 9499, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11946, 11940, 11983, 12018, 2, 9558, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 11990, 11984, 12008, 12037, 2, 9583, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12024, 12018, 12036, 12048, 2, 9622, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12042, 12037, 12045, 12054, 2, 9631, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12208, 12205, 12445, 12626, 3, 10564, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12451, 12446, 12556, 12928, 3, 10887, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12561, 12557, 12586, 12993, 3, 10952, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12632, 12626, 12927, 13028, 3, 11616, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12933, 12928, 12992, 13113, 3, 11716, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 12997, 12993, 13006, 13147, 3, 11750, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13034, 13028, 13112, 13162, 3, 11860, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13118, 13113, 13146, 13206, 3, 11914, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13168, 13162, 13205, 13226, 3, 11984, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13211, 13206, 13220, 13245, 3, 12009, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13232, 13226, 13244, 13251, 3, 12037, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13316, 13314, 13405, 13506, 4, 12446, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13410, 13406, 13460, 13617, 4, 12557, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13511, 13506, 13616, 13671, 4, 12928, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13621, 13617, 13646, 13736, 4, 12993, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13676, 13671, 13735, 13761, 4, 13113, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13740, 13736, 13749, 13795, 4, 13147, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13766, 13761, 13794, 13804, 4, 13206, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_p(pbuffer, 13809, 13804, 13818, 13821, 4, 13245, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 2462, 2460, 3539, 4555, 4644, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 3545, 3540, 4146, 4646, 4744, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4155, 4147, 4435, 4746, 4813, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4444, 4436, 4504, 4815, 4848, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4560, 4555, 4645, 4874, 4923, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4654, 4646, 4745, 4925, 4977, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4754, 4746, 4814, 4979, 5012, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4823, 4815, 4849, 5014, 5027, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4879, 4874, 4924, 5039, 5062, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4933, 4925, 4978, 5064, 5087, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 4987, 4979, 5013, 5089, 5102, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5044, 5039, 5063, 5112, 5120, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5072, 5064, 5088, 5122, 5129, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5117, 5112, 5121, 5137, 5138, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 5792, 5787, 6991, 7987, 9119, 1, 3540, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7000, 6992, 7791, 9121, 9219, 1, 4147, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7800, 7792, 7912, 9221, 9288, 1, 4436, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7920, 7913, 7946, 9290, 9303, 1, 4505, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 7995, 7987, 9120, 9326, 9409, 1, 4646, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9129, 9121, 9220, 9411, 9463, 1, 4746, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9229, 9221, 9289, 9465, 9498, 1, 4815, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9334, 9326, 9410, 9515, 9557, 1, 4925, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9419, 9411, 9464, 9559, 9582, 1, 4979, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9473, 9465, 9499, 9584, 9597, 1, 5014, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9523, 9515, 9558, 9604, 9621, 1, 5064, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9567, 9559, 9583, 9623, 9630, 1, 5089, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9612, 9604, 9622, 9638, 9642, 1, 5122, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 9987, 9982, 10564, 11003, 11615, 2, 6992, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 10573, 10565, 10887, 11617, 11715, 2, 7792, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 10895, 10888, 10952, 11717, 11749, 2, 7913, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11011, 11003, 11616, 11776, 11859, 2, 9121, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11625, 11617, 11716, 11861, 11913, 2, 9221, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11724, 11717, 11750, 11915, 11928, 2, 9290, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11784, 11776, 11860, 11941, 11983, 2, 9411, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11869, 11861, 11914, 11985, 12008, 2, 9465, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11949, 11941, 11984, 12019, 12036, 2, 9559, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 11993, 11985, 12009, 12038, 12045, 2, 9584, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 12027, 12019, 12037, 12049, 12053, 2, 9623, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 12211, 12206, 12446, 12627, 12927, 3, 10565, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 12454, 12447, 12557, 12929, 12992, 3, 10888, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 12635, 12627, 12928, 13029, 13112, 3, 11617, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 12936, 12929, 12993, 13114, 13146, 3, 11717, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13037, 13029, 13113, 13163, 13205, 3, 11861, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13121, 13114, 13147, 13207, 13220, 3, 11915, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13171, 13163, 13206, 13227, 13244, 3, 11985, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13235, 13227, 13245, 13252, 13256, 3, 12038, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13319, 13315, 13406, 13507, 13616, 4, 12447, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13514, 13507, 13617, 13672, 13735, 4, 12929, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13679, 13672, 13736, 13762, 13794, 4, 13114, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_d(pbuffer, 13769, 13762, 13795, 13805, 13818, 4, 13207, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 1, 0, 2460, 4519, 4554, 7985, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 2465, 2461, 3540, 4556, 4645, 9119, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 3548, 3541, 4147, 4647, 4745, 9219, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4158, 4148, 4436, 4747, 4814, 9288, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4447, 4437, 4505, 4816, 4849, 9303, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4520, 4519, 4555, 4853, 4873, 9324, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4563, 4556, 4646, 4875, 4924, 9409, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4657, 4647, 4746, 4926, 4978, 9463, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4757, 4747, 4815, 4980, 5013, 9498, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4854, 4853, 4874, 5028, 5038, 9513, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4882, 4875, 4925, 5040, 5063, 9557, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 4936, 4926, 4979, 5065, 5088, 9582, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5029, 5028, 5039, 5107, 5111, 9602, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5047, 5040, 5064, 5113, 5121, 9621, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5108, 5107, 5112, 5135, 5136, 9636, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5142, 5141, 5787, 7950, 7986, 11001, 1, 2461, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 5795, 5788, 6992, 7988, 9120, 11615, 1, 3541, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7003, 6993, 7792, 9122, 9220, 11715, 1, 4148, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7803, 7793, 7913, 9222, 9289, 11749, 1, 4437, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7951, 7950, 7987, 9304, 9325, 11774, 1, 4556, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 7998, 7988, 9121, 9327, 9410, 11859, 1, 4647, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9132, 9122, 9221, 9412, 9464, 11913, 1, 4747, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9232, 9222, 9290, 9466, 9499, 11928, 1, 4816, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9305, 9304, 9326, 9503, 9514, 11939, 1, 4875, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9337, 9327, 9411, 9516, 9558, 11983, 1, 4926, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9422, 9412, 9465, 9560, 9583, 12008, 1, 4980, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9504, 9503, 9515, 9598, 9603, 12017, 1, 5040, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9526, 9516, 9559, 9605, 9622, 12036, 1, 5065, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9599, 9598, 9604, 9635, 9637, 12047, 1, 5113, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9647, 9646, 9982, 10966, 11002, 12625, 2, 5788, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 9990, 9983, 10565, 11004, 11616, 12927, 2, 6993, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 10576, 10566, 10888, 11618, 11716, 12992, 2, 7793, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 10967, 10966, 11003, 11754, 11775, 13027, 2, 7988, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11014, 11004, 11617, 11777, 11860, 13112, 2, 9122, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11628, 11618, 11717, 11862, 11914, 13146, 2, 9222, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11755, 11754, 11776, 11929, 11940, 13161, 2, 9327, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11787, 11777, 11861, 11942, 11984, 13205, 2, 9412, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11872, 11862, 11915, 11986, 12009, 13220, 2, 9466, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11930, 11929, 11941, 12013, 12018, 13225, 2, 9516, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 11952, 11942, 11985, 12020, 12037, 13244, 2, 9560, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 12014, 12013, 12019, 12046, 12048, 13250, 2, 9605, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 12057, 12056, 12206, 12590, 12626, 13505, 3, 9983, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 12214, 12207, 12447, 12628, 12928, 13616, 3, 10566, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 12591, 12590, 12627, 13007, 13028, 13670, 3, 11004, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 12638, 12628, 12929, 13030, 13113, 13735, 3, 11618, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13008, 13007, 13029, 13151, 13162, 13760, 3, 11777, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13040, 13030, 13114, 13164, 13206, 13794, 3, 11862, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13152, 13151, 13163, 13221, 13226, 13803, 3, 11942, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13174, 13164, 13207, 13228, 13245, 13818, 3, 11986, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13222, 13221, 13227, 13249, 13251, 13820, 3, 12020, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13259, 13258, 13315, 13470, 13506, 13826, 4, 12207, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13471, 13470, 13507, 13650, 13671, 13861, 4, 12628, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13651, 13650, 13672, 13750, 13761, 13881, 4, 13030, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13751, 13750, 13762, 13799, 13804, 13891, 4, 13164, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sp_f(pbuffer, 13800, 13799, 13805, 13819, 13821, 13895, 4, 13228, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4161, 4145, 4149, 4744, 4748, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4450, 4434, 4438, 4813, 4817, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4512, 4503, 4506, 4848, 4850, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4660, 4644, 4648, 4923, 4927, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4760, 4744, 4748, 4977, 4981, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4826, 4813, 4817, 5012, 5015, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4939, 4923, 4927, 5062, 5066, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 4990, 4977, 4981, 5087, 5090, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5021, 5012, 5015, 5102, 5104, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5075, 5062, 5066, 5120, 5123, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 5096, 5087, 5090, 5129, 5131, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7006, 6990, 6994, 9119, 9123, 1, 4145, 4149, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7806, 7790, 7794, 9219, 9223, 1, 4434, 4438, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 7923, 7911, 7914, 9288, 9291, 1, 4503, 4506, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 8001, 7985, 7989, 9324, 9328, 1, 4644, 4648, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9135, 9119, 9123, 9409, 9413, 1, 4744, 4748, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9235, 9219, 9223, 9463, 9467, 1, 4813, 4817, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9297, 9288, 9291, 9498, 9500, 1, 4848, 4850, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9340, 9324, 9328, 9513, 9517, 1, 4923, 4927, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9425, 9409, 9413, 9557, 9561, 1, 4977, 4981, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9476, 9463, 9467, 9582, 9585, 1, 5012, 5015, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9529, 9513, 9517, 9602, 9606, 1, 5062, 5066, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9570, 9557, 9561, 9621, 9624, 1, 5087, 5090, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9591, 9582, 9585, 9630, 9632, 1, 5102, 5104, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 9615, 9602, 9606, 9636, 9639, 1, 5120, 5123, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 10579, 10563, 10567, 11615, 11619, 2, 7790, 7794, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 10898, 10886, 10889, 11715, 11718, 2, 7911, 7914, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 10959, 10951, 10953, 11749, 11751, 2, 7945, 7947, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11017, 11001, 11005, 11774, 11778, 2, 9119, 9123, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11631, 11615, 11619, 11859, 11863, 2, 9219, 9223, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11727, 11715, 11718, 11913, 11916, 2, 9288, 9291, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11790, 11774, 11778, 11939, 11943, 2, 9409, 9413, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11875, 11859, 11863, 11983, 11987, 2, 9463, 9467, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11922, 11913, 11916, 12008, 12010, 2, 9498, 9500, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11955, 11939, 11943, 12017, 12021, 2, 9557, 9561, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 11996, 11983, 11987, 12036, 12039, 2, 9582, 9585, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 12030, 12017, 12021, 12047, 12050, 2, 9621, 9624, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 12457, 12445, 12448, 12927, 12930, 3, 10886, 10889, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 12564, 12556, 12558, 12992, 12994, 3, 10951, 10953, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 12641, 12625, 12629, 13027, 13031, 3, 11615, 11619, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 12939, 12927, 12930, 13112, 13115, 3, 11715, 11718, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13000, 12992, 12994, 13146, 13148, 3, 11749, 11751, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13043, 13027, 13031, 13161, 13165, 3, 11859, 11863, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13124, 13112, 13115, 13205, 13208, 3, 11913, 11916, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13177, 13161, 13165, 13225, 13229, 3, 11983, 11987, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13214, 13205, 13208, 13244, 13246, 3, 12008, 12010, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13238, 13225, 13229, 13250, 13253, 3, 12036, 12039, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13413, 13405, 13407, 13616, 13618, 4, 12556, 12558, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13464, 13460, 13461, 13646, 13647, 4, 12586, 12587, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13517, 13505, 13508, 13670, 13673, 4, 12927, 12930, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13624, 13616, 13618, 13735, 13737, 4, 12992, 12994, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13682, 13670, 13673, 13760, 13763, 4, 13112, 13115, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13743, 13735, 13737, 13794, 13796, 4, 13146, 13148, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13772, 13760, 13763, 13803, 13806, 4, 13205, 13208, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13812, 13803, 13806, 13820, 13822, 4, 13244, 13246, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13830, 13826, 13827, 13861, 13862, 5, 13616, 13618, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13865, 13861, 13862, 13881, 13882, 5, 13735, 13737, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_s(pbuffer, 13885, 13881, 13882, 13891, 13892, 5, 13794, 13796, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 3551, 3539, 3542, 4149, 4645, 4651, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4167, 4146, 4152, 4438, 4745, 4751, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4456, 4435, 4441, 4506, 4814, 4820, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4566, 4554, 4557, 4648, 4873, 4876, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4666, 4645, 4651, 4748, 4924, 4930, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4766, 4745, 4751, 4817, 4978, 4984, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4832, 4814, 4820, 4850, 5013, 5018, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4885, 4873, 4876, 4927, 5038, 5041, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4945, 4924, 4930, 4981, 5063, 5069, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 4996, 4978, 4984, 5015, 5088, 5093, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5050, 5038, 5041, 5066, 5111, 5114, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5081, 5063, 5069, 5090, 5121, 5126, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 5798, 5786, 5789, 6994, 7986, 7992, 1, 3539, 3542, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7012, 6991, 6997, 7794, 9120, 9126, 1, 4146, 4152, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7812, 7791, 7797, 7914, 9220, 9226, 1, 4435, 4441, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 7929, 7912, 7917, 7947, 9289, 9294, 1, 4504, 4509, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 8007, 7986, 7992, 9123, 9325, 9331, 1, 4645, 4651, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9141, 9120, 9126, 9223, 9410, 9416, 1, 4745, 4751, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9241, 9220, 9226, 9291, 9464, 9470, 1, 4814, 4820, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9346, 9325, 9331, 9413, 9514, 9520, 1, 4924, 4930, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9431, 9410, 9416, 9467, 9558, 9564, 1, 4978, 4984, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9482, 9464, 9470, 9500, 9583, 9588, 1, 5013, 5018, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9535, 9514, 9520, 9561, 9603, 9609, 1, 5063, 5069, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9576, 9558, 9564, 9585, 9622, 9627, 1, 5088, 5093, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 9993, 9981, 9984, 10567, 11002, 11008, 2, 6991, 6997, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 10585, 10564, 10570, 10889, 11616, 11622, 2, 7791, 7797, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 10904, 10887, 10892, 10953, 11716, 11721, 2, 7912, 7917, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11023, 11002, 11008, 11619, 11775, 11781, 2, 9120, 9126, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11637, 11616, 11622, 11718, 11860, 11866, 2, 9220, 9226, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11733, 11716, 11721, 11751, 11914, 11919, 2, 9289, 9294, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11796, 11775, 11781, 11863, 11940, 11946, 2, 9410, 9416, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11881, 11860, 11866, 11916, 11984, 11990, 2, 9464, 9470, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 11961, 11940, 11946, 11987, 12018, 12024, 2, 9558, 9564, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 12002, 11984, 11990, 12010, 12037, 12042, 2, 9583, 9588, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 12217, 12205, 12208, 12448, 12626, 12632, 3, 10564, 10570, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 12463, 12446, 12451, 12558, 12928, 12933, 3, 10887, 10892, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 12570, 12557, 12561, 12587, 12993, 12997, 3, 10952, 10956, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 12647, 12626, 12632, 12930, 13028, 13034, 3, 11616, 11622, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 12945, 12928, 12933, 12994, 13113, 13118, 3, 11716, 11721, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13049, 13028, 13034, 13115, 13162, 13168, 3, 11860, 11866, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13130, 13113, 13118, 13148, 13206, 13211, 3, 11914, 11919, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13183, 13162, 13168, 13208, 13226, 13232, 3, 11984, 11990, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13322, 13314, 13316, 13407, 13506, 13511, 4, 12446, 12451, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13419, 13406, 13410, 13461, 13617, 13621, 4, 12557, 12561, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13523, 13506, 13511, 13618, 13671, 13676, 4, 12928, 12933, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13630, 13617, 13621, 13647, 13736, 13740, 4, 12993, 12997, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13688, 13671, 13676, 13737, 13761, 13766, 4, 13113, 13118, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_p(pbuffer, 13778, 13761, 13766, 13796, 13804, 13809, 4, 13206, 13211, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 2468, 2460, 2462, 3542, 4555, 4560, 4644, 4648, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 3557, 3540, 3545, 4152, 4646, 4654, 4744, 4748, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4173, 4147, 4155, 4441, 4746, 4754, 4813, 4817, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4462, 4436, 4444, 4509, 4815, 4823, 4848, 4850, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4572, 4555, 4560, 4651, 4874, 4879, 4923, 4927, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4672, 4646, 4654, 4751, 4925, 4933, 4977, 4981, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4772, 4746, 4754, 4820, 4979, 4987, 5012, 5015, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4891, 4874, 4879, 4930, 5039, 5044, 5062, 5066, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 4951, 4925, 4933, 4984, 5064, 5072, 5087, 5090, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5056, 5039, 5044, 5069, 5112, 5117, 5120, 5123, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 5804, 5787, 5792, 6997, 7987, 7995, 9119, 9123, 1, 3540, 3545, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7018, 6992, 7000, 7797, 9121, 9129, 9219, 9223, 1, 4147, 4155, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 7818, 7792, 7800, 7917, 9221, 9229, 9288, 9291, 1, 4436, 4444, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 8013, 7987, 7995, 9126, 9326, 9334, 9409, 9413, 1, 4646, 4654, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 9147, 9121, 9129, 9226, 9411, 9419, 9463, 9467, 1, 4746, 4754, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 9247, 9221, 9229, 9294, 9465, 9473, 9498, 9500, 1, 4815, 4823, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 9352, 9326, 9334, 9416, 9515, 9523, 9557, 9561, 1, 4925, 4933, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 9437, 9411, 9419, 9470, 9559, 9567, 9582, 9585, 1, 4979, 4987, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 9541, 9515, 9523, 9564, 9604, 9612, 9621, 9624, 1, 5064, 5072, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 9999, 9982, 9987, 10570, 11003, 11011, 11615, 11619, 2, 6992, 7000, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 10591, 10565, 10573, 10892, 11617, 11625, 11715, 11718, 2, 7792, 7800, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 10910, 10888, 10895, 10956, 11717, 11724, 11749, 11751, 2, 7913, 7920, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 11029, 11003, 11011, 11622, 11776, 11784, 11859, 11863, 2, 9121, 9129, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 11643, 11617, 11625, 11721, 11861, 11869, 11913, 11916, 2, 9221, 9229, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 11802, 11776, 11784, 11866, 11941, 11949, 11983, 11987, 2, 9411, 9419, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 11887, 11861, 11869, 11919, 11985, 11993, 12008, 12010, 2, 9465, 9473, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 11967, 11941, 11949, 11990, 12019, 12027, 12036, 12039, 2, 9559, 9567, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 12223, 12206, 12211, 12451, 12627, 12635, 12927, 12930, 3, 10565, 10573, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 12469, 12447, 12454, 12561, 12929, 12936, 12992, 12994, 3, 10888, 10895, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 12653, 12627, 12635, 12933, 13029, 13037, 13112, 13115, 3, 11617, 11625, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 12951, 12929, 12936, 12997, 13114, 13121, 13146, 13148, 3, 11717, 11724, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 13055, 13029, 13037, 13118, 13163, 13171, 13205, 13208, 3, 11861, 11869, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 13189, 13163, 13171, 13211, 13227, 13235, 13244, 13246, 3, 11985, 11993, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 13328, 13315, 13319, 13410, 13507, 13514, 13616, 13618, 4, 12447, 12454, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 13529, 13507, 13514, 13621, 13672, 13679, 13735, 13737, 4, 12929, 12936, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_d(pbuffer, 13694, 13672, 13679, 13740, 13762, 13769, 13794, 13796, 4, 13114, 13121, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4, 0, 1, 2462, 4519, 4520, 4554, 4557, 7989, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 2474, 2461, 2465, 3545, 4556, 4563, 4645, 4651, 9123, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 3563, 3541, 3548, 4155, 4647, 4657, 4745, 4751, 9223, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4179, 4148, 4158, 4444, 4747, 4757, 4814, 4820, 9291, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4523, 4519, 4520, 4560, 4853, 4854, 4873, 4876, 9328, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4578, 4556, 4563, 4654, 4875, 4882, 4924, 4930, 9413, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4678, 4647, 4657, 4754, 4926, 4936, 4978, 4984, 9467, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4857, 4853, 4854, 4879, 5028, 5029, 5038, 5041, 9517, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 4897, 4875, 4882, 4933, 5040, 5047, 5063, 5069, 9561, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5032, 5028, 5029, 5044, 5107, 5108, 5111, 5114, 9606, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5145, 5141, 5142, 5792, 7950, 7951, 7986, 7992, 11005, 1, 2461, 2465, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 5810, 5788, 5795, 7000, 7988, 7998, 9120, 9126, 11619, 1, 3541, 3548, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 7024, 6993, 7003, 7800, 9122, 9132, 9220, 9226, 11718, 1, 4148, 4158, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 7824, 7793, 7803, 7920, 9222, 9232, 9289, 9294, 11751, 1, 4437, 4447, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 7954, 7950, 7951, 7995, 9304, 9305, 9325, 9331, 11778, 1, 4556, 4563, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 8019, 7988, 7998, 9129, 9327, 9337, 9410, 9416, 11863, 1, 4647, 4657, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 9153, 9122, 9132, 9229, 9412, 9422, 9464, 9470, 11916, 1, 4747, 4757, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 9308, 9304, 9305, 9334, 9503, 9504, 9514, 9520, 11943, 1, 4875, 4882, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 9358, 9327, 9337, 9419, 9516, 9526, 9558, 9564, 11987, 1, 4926, 4936, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 9507, 9503, 9504, 9523, 9598, 9599, 9603, 9609, 12021, 1, 5040, 5047, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 9650, 9646, 9647, 9987, 10966, 10967, 11002, 11008, 12629, 2, 5788, 5795, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 10005, 9983, 9990, 10573, 11004, 11014, 11616, 11622, 12930, 2, 6993, 7003, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 10597, 10566, 10576, 10895, 11618, 11628, 11716, 11721, 12994, 2, 7793, 7803, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 10970, 10966, 10967, 11011, 11754, 11755, 11775, 11781, 13031, 2, 7988, 7998, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 11035, 11004, 11014, 11625, 11777, 11787, 11860, 11866, 13115, 2, 9122, 9132, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 11649, 11618, 11628, 11724, 11862, 11872, 11914, 11919, 13148, 2, 9222, 9232, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 11758, 11754, 11755, 11784, 11929, 11930, 11940, 11946, 13165, 2, 9327, 9337, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 11808, 11777, 11787, 11869, 11942, 11952, 11984, 11990, 13208, 2, 9412, 9422, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 11933, 11929, 11930, 11949, 12013, 12014, 12018, 12024, 13229, 2, 9516, 9526, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 12060, 12056, 12057, 12211, 12590, 12591, 12626, 12632, 13508, 3, 9983, 9990, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 12229, 12207, 12214, 12454, 12628, 12638, 12928, 12933, 13618, 3, 10566, 10576, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 12594, 12590, 12591, 12635, 13007, 13008, 13028, 13034, 13673, 3, 11004, 11014, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 12659, 12628, 12638, 12936, 13030, 13040, 13113, 13118, 13737, 3, 11618, 11628, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13011, 13007, 13008, 13037, 13151, 13152, 13162, 13168, 13763, 3, 11777, 11787, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13061, 13030, 13040, 13121, 13164, 13174, 13206, 13211, 13796, 3, 11862, 11872, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13155, 13151, 13152, 13171, 13221, 13222, 13226, 13232, 13806, 3, 11942, 11952, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13262, 13258, 13259, 13319, 13470, 13471, 13506, 13511, 13827, 4, 12207, 12214, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13474, 13470, 13471, 13514, 13650, 13651, 13671, 13676, 13862, 4, 12628, 12638, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13654, 13650, 13651, 13679, 13750, 13751, 13761, 13766, 13882, 4, 13030, 13040, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sd_f(pbuffer, 13754, 13750, 13751, 13769, 13799, 13800, 13804, 13809, 13892, 4, 13164, 13174, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4185, 4149, 4161, 4748, 4760, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4468, 4438, 4450, 4817, 4826, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4684, 4648, 4660, 4927, 4939, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4778, 4748, 4760, 4981, 4990, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4838, 4817, 4826, 5015, 5021, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 4957, 4927, 4939, 5066, 5075, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 5002, 4981, 4990, 5090, 5096, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7030, 6994, 7006, 9123, 9135, 1, 4149, 4161, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7830, 7794, 7806, 9223, 9235, 1, 4438, 4450, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 7935, 7914, 7923, 9291, 9297, 1, 4506, 4512, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 8025, 7989, 8001, 9328, 9340, 1, 4648, 4660, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 9159, 9123, 9135, 9413, 9425, 1, 4748, 4760, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 9253, 9223, 9235, 9467, 9476, 1, 4817, 4826, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 9364, 9328, 9340, 9517, 9529, 1, 4927, 4939, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 9443, 9413, 9425, 9561, 9570, 1, 4981, 4990, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 9488, 9467, 9476, 9585, 9591, 1, 5015, 5021, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 9547, 9517, 9529, 9606, 9615, 1, 5066, 5075, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 10603, 10567, 10579, 11619, 11631, 2, 7794, 7806, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 10916, 10889, 10898, 11718, 11727, 2, 7914, 7923, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 11041, 11005, 11017, 11778, 11790, 2, 9123, 9135, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 11655, 11619, 11631, 11863, 11875, 2, 9223, 9235, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 11739, 11718, 11727, 11916, 11922, 2, 9291, 9297, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 11814, 11778, 11790, 11943, 11955, 2, 9413, 9425, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 11893, 11863, 11875, 11987, 11996, 2, 9467, 9476, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 11973, 11943, 11955, 12021, 12030, 2, 9561, 9570, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 12475, 12448, 12457, 12930, 12939, 3, 10889, 10898, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 12576, 12558, 12564, 12994, 13000, 3, 10953, 10959, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 12665, 12629, 12641, 13031, 13043, 3, 11619, 11631, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 12957, 12930, 12939, 13115, 13124, 3, 11718, 11727, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13067, 13031, 13043, 13165, 13177, 3, 11863, 11875, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13136, 13115, 13124, 13208, 13214, 3, 11916, 11922, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13195, 13165, 13177, 13229, 13238, 3, 11987, 11996, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13425, 13407, 13413, 13618, 13624, 4, 12558, 12564, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13535, 13508, 13517, 13673, 13682, 4, 12930, 12939, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13636, 13618, 13624, 13737, 13743, 4, 12994, 13000, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13700, 13673, 13682, 13763, 13772, 4, 13115, 13124, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13784, 13763, 13772, 13806, 13812, 4, 13208, 13214, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13836, 13827, 13830, 13862, 13865, 5, 13618, 13624, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_s(pbuffer, 13871, 13862, 13865, 13882, 13885, 5, 13737, 13743, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 3569, 3542, 3551, 4161, 4651, 4666, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4195, 4152, 4167, 4450, 4751, 4766, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4478, 4441, 4456, 4512, 4820, 4832, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4584, 4557, 4566, 4660, 4876, 4885, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4694, 4651, 4666, 4760, 4930, 4945, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4788, 4751, 4766, 4826, 4984, 4996, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4903, 4876, 4885, 4939, 5041, 5050, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 4967, 4930, 4945, 4990, 5069, 5081, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 5816, 5789, 5798, 7006, 7992, 8007, 1, 3542, 3551, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 7040, 6997, 7012, 7806, 9126, 9141, 1, 4152, 4167, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 7840, 7797, 7812, 7923, 9226, 9241, 1, 4441, 4456, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 8035, 7992, 8007, 9135, 9331, 9346, 1, 4651, 4666, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 9169, 9126, 9141, 9235, 9416, 9431, 1, 4751, 4766, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 9263, 9226, 9241, 9297, 9470, 9482, 1, 4820, 4832, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 9374, 9331, 9346, 9425, 9520, 9535, 1, 4930, 4945, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 9453, 9416, 9431, 9476, 9564, 9576, 1, 4984, 4996, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 10011, 9984, 9993, 10579, 11008, 11023, 2, 6997, 7012, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 10613, 10570, 10585, 10898, 11622, 11637, 2, 7797, 7812, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 10926, 10892, 10904, 10959, 11721, 11733, 2, 7917, 7929, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 11051, 11008, 11023, 11631, 11781, 11796, 2, 9126, 9141, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 11665, 11622, 11637, 11727, 11866, 11881, 2, 9226, 9241, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 11824, 11781, 11796, 11875, 11946, 11961, 2, 9416, 9431, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 11903, 11866, 11881, 11922, 11990, 12002, 2, 9470, 9482, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 12235, 12208, 12217, 12457, 12632, 12647, 3, 10570, 10585, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 12485, 12451, 12463, 12564, 12933, 12945, 3, 10892, 10904, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 12675, 12632, 12647, 12939, 13034, 13049, 3, 11622, 11637, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 12967, 12933, 12945, 13000, 13118, 13130, 3, 11721, 11733, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 13077, 13034, 13049, 13124, 13168, 13183, 3, 11866, 11881, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 13334, 13316, 13322, 13413, 13511, 13523, 4, 12451, 12463, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 13435, 13410, 13419, 13464, 13621, 13630, 4, 12561, 12570, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 13545, 13511, 13523, 13624, 13676, 13688, 4, 12933, 12945, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_p(pbuffer, 13710, 13676, 13688, 13743, 13766, 13778, 4, 13118, 13130, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 2480, 2462, 2468, 3551, 4560, 4572, 4648, 4660, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 3579, 3545, 3557, 4167, 4654, 4672, 4748, 4760, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4205, 4155, 4173, 4456, 4754, 4772, 4817, 4826, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4594, 4560, 4572, 4666, 4879, 4891, 4927, 4939, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4704, 4654, 4672, 4766, 4933, 4951, 4981, 4990, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 4913, 4879, 4891, 4945, 5044, 5056, 5066, 5075, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 5826, 5792, 5804, 7012, 7995, 8013, 9123, 9135, 1, 3545, 3557, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 7050, 7000, 7018, 7812, 9129, 9147, 9223, 9235, 1, 4155, 4173, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 7850, 7800, 7818, 7929, 9229, 9247, 9291, 9297, 1, 4444, 4462, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 8045, 7995, 8013, 9141, 9334, 9352, 9413, 9425, 1, 4654, 4672, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 9179, 9129, 9147, 9241, 9419, 9437, 9467, 9476, 1, 4754, 4772, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 9384, 9334, 9352, 9431, 9523, 9541, 9561, 9570, 1, 4933, 4951, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 10021, 9987, 9999, 10585, 11011, 11029, 11619, 11631, 2, 7000, 7018, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 10623, 10573, 10591, 10904, 11625, 11643, 11718, 11727, 2, 7800, 7818, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 11061, 11011, 11029, 11637, 11784, 11802, 11863, 11875, 2, 9129, 9147, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 11675, 11625, 11643, 11733, 11869, 11887, 11916, 11922, 2, 9229, 9247, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 11834, 11784, 11802, 11881, 11949, 11967, 11987, 11996, 2, 9419, 9437, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 12245, 12211, 12223, 12463, 12635, 12653, 12930, 12939, 3, 10573, 10591, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 12495, 12454, 12469, 12570, 12936, 12951, 12994, 13000, 3, 10895, 10910, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 12685, 12635, 12653, 12945, 13037, 13055, 13115, 13124, 3, 11625, 11643, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 13087, 13037, 13055, 13130, 13171, 13189, 13208, 13214, 3, 11869, 11887, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 13344, 13319, 13328, 13419, 13514, 13529, 13618, 13624, 4, 12454, 12469, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_d(pbuffer, 13555, 13514, 13529, 13630, 13679, 13694, 13737, 13743, 4, 12936, 12951, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 10, 1, 4, 2468, 4520, 4523, 4557, 4566, 8001, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 2490, 2465, 2474, 3557, 4563, 4578, 4651, 4666, 9135, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 3589, 3548, 3563, 4173, 4657, 4678, 4751, 4766, 9235, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 4529, 4520, 4523, 4572, 4854, 4857, 4876, 4885, 9340, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 4604, 4563, 4578, 4672, 4882, 4897, 4930, 4945, 9425, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 4863, 4854, 4857, 4891, 5029, 5032, 5041, 5050, 9529, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 5151, 5142, 5145, 5804, 7951, 7954, 7992, 8007, 11017, 1, 2465, 2474, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 5836, 5795, 5810, 7018, 7998, 8019, 9126, 9141, 11631, 1, 3548, 3563, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 7060, 7003, 7024, 7818, 9132, 9153, 9226, 9241, 11727, 1, 4158, 4179, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 7960, 7951, 7954, 8013, 9305, 9308, 9331, 9346, 11790, 1, 4563, 4578, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 8055, 7998, 8019, 9147, 9337, 9358, 9416, 9431, 11875, 1, 4657, 4678, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 9314, 9305, 9308, 9352, 9504, 9507, 9520, 9535, 11955, 1, 4882, 4897, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 9656, 9647, 9650, 9999, 10967, 10970, 11008, 11023, 12641, 2, 5795, 5810, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 10031, 9990, 10005, 10591, 11014, 11035, 11622, 11637, 12939, 2, 7003, 7024, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 10633, 10576, 10597, 10910, 11628, 11649, 11721, 11733, 13000, 2, 7803, 7824, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 10976, 10967, 10970, 11029, 11755, 11758, 11781, 11796, 13043, 2, 7998, 8019, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 11071, 11014, 11035, 11643, 11787, 11808, 11866, 11881, 13124, 2, 9132, 9153, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 11764, 11755, 11758, 11802, 11930, 11933, 11946, 11961, 13177, 2, 9337, 9358, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 12066, 12057, 12060, 12223, 12591, 12594, 12632, 12647, 13517, 3, 9990, 10005, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 12255, 12214, 12229, 12469, 12638, 12659, 12933, 12945, 13624, 3, 10576, 10597, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 12600, 12591, 12594, 12653, 13008, 13011, 13034, 13049, 13682, 3, 11014, 11035, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 12695, 12638, 12659, 12951, 13040, 13061, 13118, 13130, 13743, 3, 11628, 11649, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 13017, 13008, 13011, 13055, 13152, 13155, 13168, 13183, 13772, 3, 11787, 11808, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 13268, 13259, 13262, 13328, 13471, 13474, 13511, 13523, 13830, 4, 12214, 12229, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 13480, 13471, 13474, 13529, 13651, 13654, 13676, 13688, 13865, 4, 12638, 12659, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sf_f(pbuffer, 13660, 13651, 13654, 13694, 13751, 13754, 13766, 13778, 13885, 4, 13040, 13061, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4215, 4161, 4185, 4760, 4778, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4488, 4450, 4468, 4826, 4838, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4714, 4660, 4684, 4939, 4957, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 4798, 4760, 4778, 4990, 5002, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 7070, 7006, 7030, 9135, 9159, 1, 4161, 4185, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 7860, 7806, 7830, 9235, 9253, 1, 4450, 4468, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 8065, 8001, 8025, 9340, 9364, 1, 4660, 4684, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 9189, 9135, 9159, 9425, 9443, 1, 4760, 4778, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 9273, 9235, 9253, 9476, 9488, 1, 4826, 4838, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 9394, 9340, 9364, 9529, 9547, 1, 4939, 4957, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 10643, 10579, 10603, 11631, 11655, 2, 7806, 7830, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 10936, 10898, 10916, 11727, 11739, 2, 7923, 7935, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 11081, 11017, 11041, 11790, 11814, 2, 9135, 9159, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 11685, 11631, 11655, 11875, 11893, 2, 9235, 9253, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 11844, 11790, 11814, 11955, 11973, 2, 9425, 9443, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 12505, 12457, 12475, 12939, 12957, 3, 10898, 10916, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 12705, 12641, 12665, 13043, 13067, 3, 11631, 11655, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 12977, 12939, 12957, 13124, 13136, 3, 11727, 11739, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 13097, 13043, 13067, 13177, 13195, 3, 11875, 11893, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 13445, 13413, 13425, 13624, 13636, 4, 12564, 12576, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 13565, 13517, 13535, 13682, 13700, 4, 12939, 12957, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 13720, 13682, 13700, 13772, 13784, 4, 13124, 13136, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_s(pbuffer, 13846, 13830, 13836, 13865, 13871, 5, 13624, 13636, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 3599, 3551, 3569, 4185, 4666, 4694, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 4230, 4167, 4195, 4468, 4766, 4788, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 4614, 4566, 4584, 4684, 4885, 4903, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 4729, 4666, 4694, 4778, 4945, 4967, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 5846, 5798, 5816, 7030, 8007, 8035, 1, 3551, 3569, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 7085, 7012, 7040, 7830, 9141, 9169, 1, 4167, 4195, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 7875, 7812, 7840, 7935, 9241, 9263, 1, 4456, 4478, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 8080, 8007, 8035, 9159, 9346, 9374, 1, 4666, 4694, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 9204, 9141, 9169, 9253, 9431, 9453, 1, 4766, 4788, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 10041, 9993, 10011, 10603, 11023, 11051, 2, 7012, 7040, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 10658, 10585, 10613, 10916, 11637, 11665, 2, 7812, 7840, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 11096, 11023, 11051, 11655, 11796, 11824, 2, 9141, 9169, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 11700, 11637, 11665, 11739, 11881, 11903, 2, 9241, 9263, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 12265, 12217, 12235, 12475, 12647, 12675, 3, 10585, 10613, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 12520, 12463, 12485, 12576, 12945, 12967, 3, 10904, 10926, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 12720, 12647, 12675, 12957, 13049, 13077, 3, 11637, 11665, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 13354, 13322, 13334, 13425, 13523, 13545, 4, 12463, 12485, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_p(pbuffer, 13580, 13523, 13545, 13636, 13688, 13710, 4, 12945, 12967, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 2500, 2468, 2480, 3569, 4572, 4594, 4660, 4684, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 3614, 3557, 3579, 4195, 4672, 4704, 4760, 4778, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 4629, 4572, 4594, 4694, 4891, 4913, 4939, 4957, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 5861, 5804, 5826, 7040, 8013, 8045, 9135, 9159, 1, 3557, 3579, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 7100, 7018, 7050, 7840, 9147, 9179, 9235, 9253, 1, 4173, 4205, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 8095, 8013, 8045, 9169, 9352, 9384, 9425, 9443, 1, 4672, 4704, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 10056, 9999, 10021, 10613, 11029, 11061, 11631, 11655, 2, 7018, 7050, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 10673, 10591, 10623, 10926, 11643, 11675, 11727, 11739, 2, 7818, 7850, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 11111, 11029, 11061, 11665, 11802, 11834, 11875, 11893, 2, 9147, 9179, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 12280, 12223, 12245, 12485, 12653, 12685, 12939, 12957, 3, 10591, 10623, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 12735, 12653, 12685, 12967, 13055, 13087, 13124, 13136, 3, 11643, 11675, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_d(pbuffer, 13369, 13328, 13344, 13435, 13529, 13555, 13624, 13636, 4, 12469, 12495, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 20, 4, 10, 2480, 4523, 4529, 4566, 4584, 8025, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 2515, 2474, 2490, 3579, 4578, 4604, 4666, 4694, 9159, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 4539, 4523, 4529, 4594, 4857, 4863, 4885, 4903, 9364, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 5161, 5145, 5151, 5826, 7954, 7960, 8007, 8035, 11041, 1, 2474, 2490, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 5876, 5810, 5836, 7050, 8019, 8055, 9141, 9169, 11655, 1, 3563, 3589, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 7970, 7954, 7960, 8045, 9308, 9314, 9346, 9374, 11814, 1, 4578, 4604, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 9666, 9650, 9656, 10021, 10970, 10976, 11023, 11051, 12665, 2, 5810, 5836, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 10071, 10005, 10031, 10623, 11035, 11071, 11637, 11665, 12957, 2, 7024, 7060, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 10986, 10970, 10976, 11061, 11758, 11764, 11796, 11824, 13067, 2, 8019, 8055, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 12076, 12060, 12066, 12245, 12594, 12600, 12647, 12675, 13535, 3, 10005, 10031, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 12295, 12229, 12255, 12495, 12659, 12695, 12945, 12967, 13636, 3, 10597, 10633, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 12610, 12594, 12600, 12685, 13011, 13017, 13049, 13077, 13700, 3, 11035, 11071, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 13278, 13262, 13268, 13344, 13474, 13480, 13523, 13545, 13836, 4, 12229, 12255, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sg_f(pbuffer, 13490, 13474, 13480, 13555, 13654, 13660, 13688, 13710, 13871, 4, 12659, 12695, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 4245, 4185, 4215, 4778, 4798, 0, -1, -1, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 7115, 7030, 7070, 9159, 9189, 1, 4185, 4215, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 7890, 7830, 7860, 9253, 9273, 1, 4468, 4488, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 8110, 8025, 8065, 9364, 9394, 1, 4684, 4714, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 10688, 10603, 10643, 11655, 11685, 2, 7830, 7860, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 11126, 11041, 11081, 11814, 11844, 2, 9159, 9189, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 12535, 12475, 12505, 12957, 12977, 3, 10916, 10936, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 12750, 12665, 12705, 13067, 13097, 3, 11655, 11685, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_s(pbuffer, 13595, 13535, 13565, 13700, 13720, 4, 12957, 12977, pfactors, 2, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 3629, 3569, 3599, 4215, 4694, 4729, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 5891, 5816, 5846, 7070, 8035, 8080, 1, 3569, 3599, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 7136, 7040, 7085, 7860, 9169, 9204, 1, 4195, 4230, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 10086, 10011, 10041, 10643, 11051, 11096, 2, 7040, 7085, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 10709, 10613, 10658, 10936, 11665, 11700, 2, 7840, 7875, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 12310, 12235, 12265, 12505, 12675, 12720, 3, 10613, 10658, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_p(pbuffer, 13384, 13334, 13354, 13445, 13545, 13580, 4, 12485, 12520, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 2530, 2480, 2500, 3599, 4594, 4629, 4684, 4714, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 5912, 5826, 5861, 7085, 8045, 8095, 9159, 9189, 1, 3579, 3614, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 10107, 10021, 10056, 10658, 11061, 11111, 11655, 11685, 2, 7050, 7100, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_d(pbuffer, 12331, 12245, 12280, 12520, 12685, 12735, 12957, 12977, 3, 10623, 10673, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 35, 10, 20, 2500, 4529, 4539, 4584, 4614, 8065, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 5176, 5151, 5161, 5861, 7960, 7970, 8035, 8080, 11081, 1, 2490, 2515, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 9681, 9656, 9666, 10056, 10976, 10986, 11051, 11096, 12705, 2, 5836, 5876, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 12091, 12066, 12076, 12280, 12600, 12610, 12675, 12720, 13565, 3, 10031, 10071, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_sh_f(pbuffer, 13293, 13268, 13278, 13369, 13480, 13490, 13545, 13580, 13846, 4, 12255, 12295, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ps_s(pbuffer, 4266, 4145, 7790, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_s(pbuffer, 4269, 4149, 7794, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 3650, 3542, 4145, 4149, 6997, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pp_p(pbuffer, 7157, 6997, 7790, 7794, 10570, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 4278, 4161, 7806, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 7166, 7006, 10579, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 8131, 8001, 11017, 1, 7006, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 10730, 10579, 12457, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 11147, 11017, 12641, 1, 10579, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_s(pbuffer, 12771, 12641, 13517, 1, 12457, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 3659, 3551, 4149, 4161, 7012, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_p(pbuffer, 7184, 7012, 7794, 7806, 10585, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 2551, 2468, 3542, 3551, 5804, 7006, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 5933, 5804, 6997, 7012, 9999, 10579, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pd_d(pbuffer, 10128, 9999, 10570, 10585, 12223, 12457, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 4296, 4185, 7830, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 7202, 7030, 10603, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 8149, 8025, 11041, 1, 7030, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 10748, 10603, 12475, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 11165, 11041, 12665, 1, 10603, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_s(pbuffer, 12789, 12665, 13535, 1, 12475, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 3677, 3569, 4161, 4185, 7040, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 5951, 5816, 7006, 7030, 10011, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 7232, 7040, 7806, 7830, 10613, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 10146, 10011, 10579, 10603, 12235, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_p(pbuffer, 12352, 12235, 12457, 12475, 13334, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 2569, 2480, 3551, 3569, 5826, 7030, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 5981, 5826, 7012, 7040, 10021, 10603, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_d(pbuffer, 10176, 10021, 10585, 10613, 12245, 12475, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 56, 10, 2468, 2480, 5151, 5816, 8001, 8025, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 5197, 5151, 5804, 5826, 9656, 10011, 11017, 11041, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 9702, 9656, 9999, 10021, 12066, 12235, 12641, 12665, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pf_f(pbuffer, 12112, 12066, 12223, 12245, 13268, 13334, 13517, 13535, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 4326, 4215, 7860, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 7262, 7070, 10643, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 8179, 8065, 11081, 1, 7070, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 10778, 10643, 12505, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 11195, 11081, 12705, 1, 10643, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_s(pbuffer, 12819, 12705, 13565, 1, 12505, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 3707, 3599, 4185, 4215, 7085, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_p(pbuffer, 7307, 7085, 7830, 7860, 10658, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 2599, 2500, 3569, 3599, 5861, 7070, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 6011, 5861, 7040, 7085, 10056, 10643, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_pg_d(pbuffer, 10206, 10056, 10613, 10658, 12280, 12505, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 4371, 4245, 7890, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 7352, 7115, 10688, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 8224, 8110, 11126, 1, 7115, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 10823, 10688, 12535, 0, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 11240, 11126, 12750, 1, 10688, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_s(pbuffer, 12864, 12750, 13595, 1, 12535, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 3752, 3629, 4215, 4245, 7136, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 6056, 5891, 7070, 7115, 10086, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 7415, 7136, 7860, 7890, 10709, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 10251, 10086, 10643, 10688, 12310, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_p(pbuffer, 12382, 12310, 12505, 12535, 13384, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 2644, 2530, 3599, 3629, 5912, 7115, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 6119, 5912, 7085, 7136, 10107, 10688, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_d(pbuffer, 10314, 10107, 10658, 10709, 12331, 12535, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 86, 35, 2500, 2530, 5176, 5891, 8065, 8110, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 5227, 5176, 5861, 5912, 9681, 10086, 11081, 11126, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 9732, 9681, 10056, 10107, 12091, 12310, 12705, 12750, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ph_f(pbuffer, 12142, 12091, 12280, 12331, 13293, 13384, 13565, 13595, 0, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dp_p(pbuffer, 3815, 3542, 3650, 4266, 4269, 6997, 7157, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 7478, 7006, 7166, 10579, 10730, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 8287, 8001, 8131, 11017, 11147, 1, 7006, 7166, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_s(pbuffer, 11303, 11017, 11147, 12641, 12771, 1, 10579, 10730, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_p(pbuffer, 3833, 3551, 3659, 4269, 4278, 7012, 7184, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_d(pbuffer, 2707, 2468, 2551, 3650, 3659, 5804, 5933, 7006, 7166, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dd_d(pbuffer, 6182, 5804, 5933, 7157, 7184, 9999, 10128, 10579, 10730, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 7514, 7030, 7202, 10603, 10748, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 8323, 8025, 8149, 11041, 11165, 1, 7030, 7202, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_s(pbuffer, 11339, 11041, 11165, 12665, 12789, 1, 10603, 10748, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 3869, 3569, 3677, 4278, 4296, 7040, 7232, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 6218, 5816, 5951, 7166, 7202, 10011, 10146, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_p(pbuffer, 10377, 10011, 10146, 10730, 10748, 12235, 12352, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 2743, 2480, 2569, 3659, 3677, 5826, 5981, 7030, 7202, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_d(pbuffer, 6278, 5826, 5981, 7184, 7232, 10021, 10176, 10603, 10748, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 149, 10, 56, 2551, 2569, 5151, 5197, 5816, 5951, 8131, 8149, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 5290, 5151, 5197, 5933, 5981, 9656, 9702, 10011, 10146, 11147, 11165, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_df_f(pbuffer, 9795, 9656, 9702, 10128, 10176, 12066, 12112, 12235, 12352, 12771, 12789, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 7574, 7070, 7262, 10643, 10778, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 8383, 8065, 8179, 11081, 11195, 1, 7070, 7262, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_s(pbuffer, 11399, 11081, 11195, 12705, 12819, 1, 10643, 10778, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_p(pbuffer, 3929, 3599, 3707, 4296, 4326, 7085, 7307, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 2803, 2500, 2599, 3677, 3707, 5861, 6011, 7070, 7262, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dg_d(pbuffer, 6338, 5861, 6011, 7232, 7307, 10056, 10206, 10643, 10778, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_s(pbuffer, 7664, 7115, 7352, 10688, 10823, 0, -1, -1, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_s(pbuffer, 8473, 8110, 8224, 11126, 11240, 1, 7115, 7352, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_s(pbuffer, 11489, 11126, 11240, 12750, 12864, 1, 10688, 10823, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_p(pbuffer, 4019, 3629, 3752, 4326, 4371, 7136, 7415, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_p(pbuffer, 6428, 5891, 6056, 7262, 7352, 10086, 10251, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_p(pbuffer, 10437, 10086, 10251, 10778, 10823, 12310, 12382, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_d(pbuffer, 2893, 2530, 2644, 3707, 3752, 5912, 6119, 7115, 7352, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_d(pbuffer, 6554, 5912, 6119, 7307, 7415, 10107, 10314, 10688, 10823, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_f(pbuffer, 209, 35, 86, 2599, 2644, 5176, 5227, 5891, 6056, 8179, 8224, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_f(pbuffer, 5350, 5176, 5227, 6011, 6119, 9681, 9732, 10086, 10251, 11195, 11240, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_dh_f(pbuffer, 9855, 9681, 9732, 10206, 10314, 12091, 12142, 12310, 12382, 12819, 12864, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_s(pbuffer, 8599, 8131, 8287, 11147, 11303, 1, 7166, 7478, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fd_d(pbuffer, 3019, 2551, 2707, 3815, 3833, 5933, 6182, 7166, 7478, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_s(pbuffer, 8659, 8149, 8323, 11165, 11339, 1, 7202, 7514, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_p(pbuffer, 6680, 5951, 6218, 7478, 7514, 10146, 10377, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_d(pbuffer, 3079, 2569, 2743, 3833, 3869, 5981, 6278, 7202, 7514, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_f(pbuffer, 335, 56, 149, 2707, 2743, 5197, 5290, 5951, 6218, 8287, 8323, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_ff_f(pbuffer, 5476, 5197, 5290, 6182, 6278, 9702, 9795, 10146, 10377, 11303, 11339, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_s(pbuffer, 8759, 8179, 8383, 11195, 11399, 1, 7262, 7574, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fg_d(pbuffer, 3179, 2599, 2803, 3869, 3929, 6011, 6338, 7262, 7574, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_s(pbuffer, 8909, 8224, 8473, 11240, 11489, 1, 7352, 7664, pfactors, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_p(pbuffer, 6780, 6056, 6428, 7574, 7664, 10251, 10437, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_d(pbuffer, 3329, 2644, 2893, 3929, 4019, 6119, 6554, 7352, 7664, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_f(pbuffer, 1785, 86, 209, 2803, 2893, 5227, 5350, 6056, 6428, 8383, 8473, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_fh_f(pbuffer, 5576, 5227, 5350, 6338, 6554, 9732, 9855, 10251, 10437, 11399, 11489, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gf_f(pbuffer, 1995, 149, 335, 3019, 3079, 5290, 5476, 6218, 6680, 8599, 8659, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2pecp::comp_prim_projected_core_potential_gh_f(pbuffer, 2145, 209, 1785, 3179, 3329, 5350, 5576, 6428, 6780, 8759, 8909, 0, -1, -1, pfactors, 2, r_a, a_exp, c_exp);

                    t2cgeom::comp_prim_op_geom_11_fg(pbuffer, 435, 149, 209, 1995, 2145, 1, pfactors, a_exp);

                    t2cfunc::reduce(cbuffer, 0, pbuffer, 435, 1350, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<3, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2pecp namespace

#endif /* ProjectedCorePotentialGeom101FGForF_hpp */
