#ifndef ElectronRepulsionGeom1000RecPPSS_hpp
#define ElectronRepulsionGeom1000RecPPSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "ElectronRepulsionGeomContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "T4CUtils.hpp"

#include "TensorComponents.hpp"
#include <iomanip>
#include <iostream>

namespace erirec {  // erirec namespace

/// @brief Computes (PP|1/|r-r'||SS)  integral derivatives for two basis function pairs blocks.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom_1000_ppss(T&                               distributor,
                                       const CGtoPairBlock&             bra_gto_pair_block,
                                       const CGtoPairBlock&             ket_gto_pair_block,
                                       const std::pair<size_t, size_t>& bra_indices,
                                       const std::pair<size_t, size_t>& ket_indices) -> void
{
    // intialize GTOs pair data on bra side

    const auto a_coords = bra_gto_pair_block.bra_coordinates();

    const auto b_coords = bra_gto_pair_block.ket_coordinates();

    const auto a_vec_exps = bra_gto_pair_block.bra_exponents();

    const auto b_vec_exps = bra_gto_pair_block.ket_exponents();

    const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();

    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();

    const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();

    const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();

    // intialize GTOs data on ket side

    const auto c_coords = ket_gto_pair_block.bra_coordinates();

    const auto d_coords = ket_gto_pair_block.ket_coordinates();

    const auto c_vec_exps = ket_gto_pair_block.bra_exponents();

    const auto d_vec_exps = ket_gto_pair_block.ket_exponents();

    const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();

    const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();

    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();

    const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(23, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(35, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(22, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(76, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(27, 1);

    // setup Boys fuction data

    const CBoysFunc<3> bf_table;

    CSimdArray<double> bf_data(5, ket_npgtos);

    // set up range seperation factor

    const auto use_rs = distributor.need_omega();

    const auto omega = distributor.get_omega();

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

        pfactors.load(c_vec_exps, ket_range, 0, ket_npgtos);

        pfactors.load(d_vec_exps, ket_range, 1, ket_npgtos);

        pfactors.load(cd_vec_ovls, ket_range, 2, ket_npgtos);

        pfactors.load(cd_vec_norms, ket_range, 3, ket_npgtos);

        pfactors.replicate_points(c_coords, ket_range, 4, ket_npgtos);

        pfactors.replicate_points(d_coords, ket_range, 7, ket_npgtos);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            skbuffer.zero();

            sbuffer.zero();

            // set up coordinates on bra side

            const auto r_a = a_coords[j];

            const auto r_b = b_coords[j];

            const auto a_xyz = r_a.coordinates();

            const auto b_xyz = r_b.coordinates();

            const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

            for (int k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = a_vec_exps[k * bra_ncgtos + j];

                const auto b_exp = b_vec_exps[k * bra_ncgtos + j];

                const auto ab_norm = ab_vec_norms[k * bra_ncgtos + j];

                const auto ab_ovl = ab_vec_ovls[k * bra_ncgtos + j];

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

                t4cfunc::comp_distances_wp(pfactors, 20, 17, r_p);

                if (use_rs)
                {
                    t4cfunc::comp_boys_args(bf_data, 4, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 4, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 4, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 4);
                }

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 4, 0, 1, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 7, 1, 2, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 10, 2, 3, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 13, 0, 1, 4, 7, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 19, 1, 2, 7, 10, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 25, 4, 7, 13, 19, pfactors, 20, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 4, 3, ket_width, ket_npgtos);
                
                pbuffer.scale(2.0 * a_exp, {4, 7});
                
                pbuffer.scale(2.0 * a_exp, {13, 19});
                
                pbuffer.scale(2.0 * a_exp, {25, 35});

                t2cfunc::reduce(cbuffer, 3, pbuffer, 4, 3, ket_width, ket_npgtos);
                
                t2cfunc::reduce(cbuffer, 6, pbuffer, 13, 6, ket_width, ket_npgtos);
                
                t2cfunc::reduce(cbuffer, 12, pbuffer, 25, 10, ket_width, ket_npgtos);
            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 1);
            
            t4cfunc::ket_transform<0, 0>(skbuffer, 3, cbuffer, 3, 0, 1);
            
            t4cfunc::ket_transform<0, 0>(skbuffer, 6, cbuffer, 6, 0, 2);
            
            t4cfunc::ket_transform<0, 0>(skbuffer, 12, cbuffer, 12, 0, 3);
            
            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 22, 3, 6, r_ab, 0, 0);
            
            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 31, 6, 12, r_ab, 0, 0);
            
            erirec::comp_bra_geom_hrr_electron_repulsion_ppxx(skbuffer, 49, 22, 31, 0, r_ab, 0, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 0, skbuffer, 49, 0, 0);
            
            t4cfunc::bra_transform<1, 1>(sbuffer, 9, skbuffer, 58, 0, 0);
            
            t4cfunc::bra_transform<1, 1>(sbuffer, 18, skbuffer, 67, 0, 0);
            
            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 1, 0, 0, j, ket_range);
            
//            // *** START DEBUG BLOCK
//
//            const auto [a_angmom, b_angmom]=  bra_gto_pair_block.angular_momentums();
//
//            const auto [c_angmom, d_angmom]=  ket_gto_pair_block.angular_momentums();
//
//            const auto adim = a_indices[0];
//
//            const auto bdim = b_indices[0];
//
//            const auto cdim = c_indices[0];
//
//            const auto ddim = d_indices[0];
//
//            // set up angular components
//
//            const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});
//
//            const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});
//
//            const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});
//
//            const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});
//
//            const auto tcomps = acomps * bcomps * ccomps * dcomps;
//
//            std::cout << std::setprecision(15);
//
//            for (size_t p = 0; p < sbuffer.number_of_active_elements(); p++)
//            {
//                for (int k = 0; k < acomps; k++)
//                {
//                    for (int l = 0; l < bcomps; l++)
//                    {
//                        for (int m = 0; m < ccomps; m++)
//                        {
//                            for (int n = 0; n < dcomps; n++)
//                            {
//                                auto idx = k * bcomps * ccomps * dcomps + l * ccomps * dcomps  + m * dcomps + n;
//                                
//                                auto tint_x = sbuffer.data(idx);
//                                
//                                auto tint_y = sbuffer.data(idx + tcomps);
//                                
//                                auto tint_z = sbuffer.data(idx + 2 * tcomps);
//                                
//                                std::cout << k * adim + a_indices[j + 1] << " " << l * bdim + b_indices[j + 1];
//                                
//                                std::cout << " " << m * cdim + c_indices[ket_indices.first + p + 1] << " " << n * ddim + d_indices[ket_indices.first + p + 1];
//                                
//                                std::cout << " " << tint_x[p] << " " << tint_y[p] << " " << tint_z[p] << std::endl;
//                            }
//                        }
//                    }
//                }
//            }
//
//            // *** END DEBUG BLOCK
        }
    }
}

}  // namespace erirec

#endif /* ElectronRepulsionGeom1000RecPPSS_hpp */