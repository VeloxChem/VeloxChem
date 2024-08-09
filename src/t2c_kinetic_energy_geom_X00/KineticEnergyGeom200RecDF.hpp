#ifndef KineticEnergyGeom200RecDF_hpp
#define KineticEnergyGeom200RecDF_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "KineticEnergyPrimRecSS.hpp"
#include "OverlapPrimRecSP.hpp"
#include "KineticEnergyPrimRecSP.hpp"
#include "OverlapPrimRecSD.hpp"
#include "KineticEnergyPrimRecSD.hpp"
#include "OverlapPrimRecSF.hpp"
#include "KineticEnergyPrimRecSF.hpp"
#include "OverlapPrimRecPS.hpp"
#include "KineticEnergyPrimRecPS.hpp"
#include "OverlapPrimRecPP.hpp"
#include "KineticEnergyPrimRecPP.hpp"
#include "OverlapPrimRecPD.hpp"
#include "KineticEnergyPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "KineticEnergyPrimRecPF.hpp"
#include "OverlapPrimRecDP.hpp"
#include "KineticEnergyPrimRecDP.hpp"
#include "OverlapPrimRecDD.hpp"
#include "KineticEnergyPrimRecDD.hpp"
#include "OverlapPrimRecDF.hpp"
#include "KineticEnergyPrimRecDF.hpp"
#include "OverlapPrimRecFD.hpp"
#include "KineticEnergyPrimRecFD.hpp"
#include "OverlapPrimRecFF.hpp"
#include "KineticEnergyPrimRecFF.hpp"
#include "OverlapPrimRecGF.hpp"
#include "KineticEnergyPrimRecGF.hpp"
#include "GeometricalDerivatives2X0ForDY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes (d^(2)/dA^(2)D|T|F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_kinetic_energy_geom_20_df(T& distributor,
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

    CSimdArray<double> pbuffer(1368, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(360, 1);

    CSimdArray<double> sbuffer(210, 1);

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

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                kinrec::comp_prim_kinetic_energy_ss(pbuffer, 1, 0, factors, a_exp);

                ovlrec::comp_prim_overlap_sp(pbuffer, 2, 0, factors, 11);

                kinrec::comp_prim_kinetic_energy_sp(pbuffer, 5, 1, 2, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sd(pbuffer, 8, 0, 2, factors, 11, a_exp);

                kinrec::comp_prim_kinetic_energy_sd(pbuffer, 14, 0, 1, 5, 8, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sf(pbuffer, 20, 2, 8, factors, 11, a_exp);

                kinrec::comp_prim_kinetic_energy_sf(pbuffer, 30, 2, 5, 14, 20, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_ps(pbuffer, 40, 0, factors, 8);

                kinrec::comp_prim_kinetic_energy_ps(pbuffer, 43, 1, 40, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pp(pbuffer, 46, 0, 2, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pp(pbuffer, 55, 1, 5, 46, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 64, 2, 8, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pd(pbuffer, 82, 5, 14, 64, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pf(pbuffer, 100, 8, 20, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pf(pbuffer, 130, 14, 30, 100, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 160, 2, 40, 46, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dp(pbuffer, 178, 2, 5, 43, 55, 160, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 196, 8, 46, 64, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dd(pbuffer, 232, 8, 14, 55, 82, 196, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_df(pbuffer, 268, 20, 64, 100, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_df(pbuffer, 328, 20, 30, 82, 130, 268, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fd(pbuffer, 388, 64, 160, 196, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fd(pbuffer, 448, 64, 82, 178, 232, 388, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ff(pbuffer, 508, 100, 196, 268, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_ff(pbuffer, 608, 100, 130, 232, 328, 508, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gf(pbuffer, 708, 268, 388, 508, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gf(pbuffer, 858, 268, 328, 448, 608, 708, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_20_dx(pbuffer, 1008, 30, 328, 858, 1, 10, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 0, ket_width, ket_npgtos);
            }

            t2cfunc::transform<2, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // kinrec namespace

#endif /* KineticEnergyGeom200RecDF_hpp */
