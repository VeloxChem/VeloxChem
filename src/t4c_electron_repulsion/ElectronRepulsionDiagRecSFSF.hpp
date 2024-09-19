#ifndef ElectronRepulsionDiagRecSFSF_hpp
#define ElectronRepulsionDiagRecSFSF_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "BoysFunc.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "T4CUtils.hpp"

namespace erirec {  // erirec namespace

/// @brief Computes (SF|1/|r-r'||SF)  integrals for GTOs pair block.
/// @param distributor The pointer to screening data distributor.
/// @param gto_pair_block The GTOs pair block.
/// @param gto_indices The range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_sfsf(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
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

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(486, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(100, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(70, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(49, 1);

    // setup Boys fuction data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, npgtos);

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

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // zero integral buffers

        cbuffer.zero();

        skbuffer.zero();

        sbuffer.zero();

        // set up coordinates on bra side

        const auto r_a = a_coords[i];

        const auto r_b = b_coords[i];

        const auto a_xyz = r_a.coordinates();

        const auto b_xyz = r_b.coordinates();

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

            t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp);

            bf_table.compute(bf_data, 0, 7);

            t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 7, 0, 1, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 1, 2, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 2, 3, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 3, 4, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 4, 5, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 5, 6, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 55, 7, 10, 25, 31, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 10, 13, 31, 37, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 13, 16, 37, 43, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 16, 19, 43, 49, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 95, 3, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 104, 13, 31, 37, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 122, 16, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 140, 31, 55, 65, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 170, 37, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 200, 43, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 230, 31, 37, 95, 104, 122, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 266, 55, 65, 104, 140, 170, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 326, 65, 75, 122, 170, 200, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 386, 140, 170, 230, 266, 326, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 386, 100, ket_width, npgtos);
        }

        t4cfunc::ket_transform<0, 3>(skbuffer, 0, cbuffer, 0, 0, 3);

        t4cfunc::bra_transform<0, 3>(sbuffer, 0, skbuffer, 0, 0, 3);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

}  // namespace erirec

#endif /* ElectronRepulsionDiagRecSFSF_hpp */
