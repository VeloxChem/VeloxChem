#ifndef ElectricDipoleMomentumRecGD_hpp
#define ElectricDipoleMomentumRecGD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "ElectricDipoleMomentumPrimRecDD.hpp"
#include "ElectricDipoleMomentumPrimRecDP.hpp"
#include "ElectricDipoleMomentumPrimRecDS.hpp"
#include "ElectricDipoleMomentumPrimRecFD.hpp"
#include "ElectricDipoleMomentumPrimRecFP.hpp"
#include "ElectricDipoleMomentumPrimRecGD.hpp"
#include "ElectricDipoleMomentumPrimRecPD.hpp"
#include "ElectricDipoleMomentumPrimRecPP.hpp"
#include "ElectricDipoleMomentumPrimRecPS.hpp"
#include "ElectricDipoleMomentumPrimRecSD.hpp"
#include "ElectricDipoleMomentumPrimRecSP.hpp"
#include "ElectricDipoleMomentumPrimRecSS.hpp"
#include "GtoBlock.hpp"
#include "OverlapPrimRecDD.hpp"
#include "OverlapPrimRecDP.hpp"
#include "OverlapPrimRecFD.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPP.hpp"
#include "OverlapPrimRecPS.hpp"
#include "OverlapPrimRecSD.hpp"
#include "OverlapPrimRecSP.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes (G|r|D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electric_dipole_momentum_gd(T&                               distributor,
                                 const CGtoBlock&                 bra_gto_block,
                                 const CGtoBlock&                 ket_gto_block,
                                 const std::pair<size_t, size_t>& bra_indices,
                                 const std::pair<size_t, size_t>& ket_indices,
                                 const bool                       bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto r_c = distributor.coordinates()[0];

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

    CSimdArray<double> factors(20, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(994, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(270, 1);

    CSimdArray<double> sbuffer(135, 1);

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

                t2cfunc::comp_coordinates_p(factors, 8, 2, r_a, a_exp);

                t2cfunc::comp_distances_pa_from_p(factors, 11, 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14, 8, 2);

                t2cfunc::comp_distances_pc(factors, 17, 8, r_c);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                diprec::comp_prim_electric_dipole_momentum_ss(pbuffer, 1, 0, factors, 17);

                ovlrec::comp_prim_overlap_sp(pbuffer, 4, 0, factors, 14);

                diprec::comp_prim_electric_dipole_momentum_sp(pbuffer, 7, 0, 1, factors, 14, a_exp);

                ovlrec::comp_prim_overlap_sd(pbuffer, 16, 0, 4, factors, 14, a_exp);

                diprec::comp_prim_electric_dipole_momentum_sd(pbuffer, 22, 1, 4, 7, factors, 14, a_exp);

                ovlrec::comp_prim_overlap_ps(pbuffer, 40, 0, factors, 11);

                diprec::comp_prim_electric_dipole_momentum_ps(pbuffer, 43, 0, 1, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_pp(pbuffer, 52, 0, 4, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_pp(pbuffer, 61, 1, 4, 7, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 88, 4, 16, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_pd(pbuffer, 106, 7, 16, 22, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_ds(pbuffer, 160, 1, 40, 43, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 178, 4, 40, 52, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_dp(pbuffer, 196, 7, 43, 52, 61, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 250, 16, 52, 88, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_dd(pbuffer, 286, 22, 61, 88, 106, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_fp(pbuffer, 394, 61, 160, 178, 196, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_fd(pbuffer, 484, 88, 178, 250, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_fd(pbuffer, 544, 106, 196, 250, 286, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_gd(pbuffer, 724, 286, 394, 484, 544, factors, 11, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 724, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace diprec

#endif /* ElectricDipoleMomentumRecGD_hpp */
