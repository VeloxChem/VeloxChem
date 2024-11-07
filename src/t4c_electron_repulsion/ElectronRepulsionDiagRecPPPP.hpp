#ifndef ElectronRepulsionDiagRecPPPP_hpp
#define ElectronRepulsionDiagRecPPPP_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "BoysFunc.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "T4CUtils.hpp"

namespace erirec {  // erirec namespace

/// @brief Computes (PP|1/|r-r'||PP)  integrals for GTOs pair block.
/// @param distributor The pointer to screening data distributor.
/// @param gto_pair_block The GTOs pair block.
/// @param gto_indices The range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_pppp(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
{
    // intialize GTOs pair data

    const auto a_coords = gto_pair_block.bra_coordinates();

    const auto b_coords = gto_pair_block.ket_coordinates();

    const auto a_vec_exps = gto_pair_block.bra_exponents();

    const auto b_vec_exps = gto_pair_block.ket_exponents();

    const auto ab_vec_norms = gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = gto_pair_block.overlap_factors();

    const auto a_indices = gto_pair_block.bra_orbital_indices();

    const auto b_indices = gto_pair_block.ket_orbital_indices();

    const auto ncgtos = gto_pair_block.number_of_contracted_pairs();

    const auto npgtos = gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(29, npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(146, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(81, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(81, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(162, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(81, 1);

    // setup Boys fuction data

    const CBoysFunc<4> bf_table;

    CSimdArray<double> bf_data(6, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dims = gto_indices.second - gto_indices.first;

    std::vector<double> max_values(gto_dims, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices.first; i < gto_indices.second; i++)
    {
        // set up indices on ket side

        auto ket_range = std::pair<size_t, size_t>{i, i + 1};

        pfactors.load(a_vec_exps, ket_range, 0, npgtos);

        pfactors.load(b_vec_exps, ket_range, 1, npgtos);

        pfactors.load(ab_vec_ovls, ket_range, 2, npgtos);

        pfactors.load(ab_vec_norms, ket_range, 3, npgtos);

        pfactors.replicate_points(a_coords, ket_range, 4, npgtos);

        pfactors.replicate_points(b_coords, ket_range, 7, npgtos);

        cfactors.replicate_points(a_coords, ket_range, 0, 1);

        cfactors.replicate_points(b_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // zero integral buffers

        cbuffer.zero();

        ckbuffer.zero();

        skbuffer.zero();

        sbuffer.zero();

        // set up coordinates on bra side

        const auto r_a = a_coords[i];

        const auto r_b = b_coords[i];

        const auto a_xyz = r_a.coordinates();

        const auto b_xyz = r_b.coordinates();

        const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

            const auto p_x = (a_xyz[0] * a_exp + b_xyz[0] * b_exp) / (a_exp + b_exp);

            const auto p_y = (a_xyz[1] * a_exp + b_xyz[1] * b_exp) / (a_exp + b_exp);

            const auto p_z = (a_xyz[2] * a_exp + b_xyz[2] * b_exp) / (a_exp + b_exp);

            const auto r_p = TPoint<double>({p_x, p_y, p_z});

            const auto pb_x = p_x - b_xyz[0];

            const auto pb_y = p_y - b_xyz[1];

            const auto pb_z = p_z - b_xyz[2];

            const auto r_pb = TPoint<double>({pb_x, pb_y, pb_z});

            t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);

            t4cfunc::comp_distances_pq(pfactors, 13, 10, r_p);

            t4cfunc::comp_coordinates_w(pfactors, 17, 10, r_p, a_exp, b_exp);

            t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

            t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

            t4cfunc::comp_distances_wp(pfactors, 26, 17, r_p);

            t4cfunc::comp_boys_args(bf_data, 5, pfactors, 13, a_exp, b_exp);

            bf_table.compute(bf_data, 0, 5);

            t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 5, 0, 1, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 8, 1, 2, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 2, 3, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 3, 4, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 17, 0, 1, 5, 8, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 23, 1, 2, 8, 11, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 29, 2, 3, 11, 14, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 35, 1, 2, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 38, 1, 5, 8, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 47, 2, 8, 11, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 56, 8, 17, 23, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 74, 11, 23, 29, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 92, 5, 8, 35, 38, 47, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 110, 17, 23, 47, 56, 74, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 38, 9, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 9, pbuffer, 56, 18, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 27, pbuffer, 92, 18, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 45, pbuffer, 110, 36, ket_width, npgtos);
        }

        erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 0, cbuffer, 0, 9, cfactors, 6, 0, 1);

        erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 27, cbuffer, 27, 45, cfactors, 6, 0, 2);

        t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 0, 0, 1);

        t4cfunc::ket_transform<1, 1>(skbuffer, 27, ckbuffer, 27, 0, 2);

        erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 81, 0, 27, r_ab, 1, 1);

        t4cfunc::bra_transform<1, 1>(sbuffer, 0, skbuffer, 81, 1, 1);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

}  // namespace erirec

#endif /* ElectronRepulsionDiagRecPPPP_hpp */
