#ifndef KineticEnergyRecDG_hpp
#define KineticEnergyRecDG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "GtoBlock.hpp"
#include "KineticEnergyPrimRecDG.hpp"
#include "KineticEnergyPrimRecPF.hpp"
#include "KineticEnergyPrimRecPG.hpp"
#include "KineticEnergyPrimRecSD.hpp"
#include "KineticEnergyPrimRecSF.hpp"
#include "KineticEnergyPrimRecSG.hpp"
#include "KineticEnergyPrimRecSP.hpp"
#include "KineticEnergyPrimRecSS.hpp"
#include "OverlapPrimRecDG.hpp"
#include "OverlapPrimRecPF.hpp"
#include "OverlapPrimRecPG.hpp"
#include "OverlapPrimRecSD.hpp"
#include "OverlapPrimRecSF.hpp"
#include "OverlapPrimRecSG.hpp"
#include "OverlapPrimRecSP.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes (D|T|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_kinetic_energy_dg(T&                               distributor,
                       const CGtoBlock&                 bra_gto_block,
                       const CGtoBlock&                 ket_gto_block,
                       const std::pair<size_t, size_t>& bra_indices,
                       const std::pair<size_t, size_t>& ket_indices,
                       const bool                       bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(400, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(90, 1);

    CSimdArray<double> sbuffer(45, 1);

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

                ovlrec::comp_prim_overlap_sg(pbuffer, 40, 8, 20, factors, 11, a_exp);

                kinrec::comp_prim_kinetic_energy_sg(pbuffer, 55, 8, 14, 30, 40, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_pf(pbuffer, 70, 8, 20, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pf(pbuffer, 100, 14, 30, 70, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pg(pbuffer, 130, 20, 40, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pg(pbuffer, 175, 30, 55, 130, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dg(pbuffer, 220, 40, 70, 130, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dg(pbuffer, 310, 40, 55, 100, 175, 220, factors, 8, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 310, ket_width, ket_npgtos);
            }

            t2cfunc::transform<2, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 4, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace kinrec

#endif /* KineticEnergyRecDG_hpp */
