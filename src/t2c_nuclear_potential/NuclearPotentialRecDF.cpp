#include "NuclearPotentialRecDF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "T2CDistributor.hpp"

namespace npotrec { // npotrec namespace

auto
compNuclearPotentialDF(      CSubMatrix* matrix,
                       const CGtoBlock&  bra_gto_block,
                       const CGtoBlock&  ket_gto_block,
                       const bool        ang_order,
                       const int64_t     bra_first,
                       const int64_t     bra_last) -> void
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);

    const double f3_5 = std::sqrt(2.5);

    const double f3_15 = 2.0 * std::sqrt(15.0);

    const double f3_3 = std::sqrt(1.5);

    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.getCoordinates();

    const auto bra_gto_exps = bra_gto_block.getExponents();

    const auto bra_gto_norms = bra_gto_block.getNormalizationFactors();

    const auto bra_gto_indexes = bra_gto_block.getOrbitalIndexes();

    const auto bra_ncgtos = bra_gto_block.getNumberOfBasisFunctions();

    const auto bra_npgtos = bra_gto_block.getNumberOfPrimitives();

    // intialize GTOs data on ket side

    const auto ket_gto_coords = ket_gto_block.getCoordinates();

    const auto ket_gto_exps = ket_gto_block.getExponents();

    const auto ket_gto_norms = ket_gto_block.getNormalizationFactors();

    const auto ket_gto_indexes = ket_gto_block.getOrbitalIndexes();

    const auto ket_ncgtos = ket_gto_block.getNumberOfBasisFunctions();

    const auto ket_npgtos = ket_gto_block.getNumberOfPrimitives();

    // initialize aligned arrays for ket side

    alignas(64) TDoubleArray ket_coords_x;

    alignas(64) TDoubleArray ket_coords_y;

    alignas(64) TDoubleArray ket_coords_z;

    alignas(64) TDoubleArray ket_exps;

    alignas(64) TDoubleArray ket_norms;

    // initialize contracted integrals buffer

    alignas(64) TDoubleArray buffer_xx;

    alignas(64) TDoubleArray buffer_xy;

    alignas(64) TDoubleArray buffer_xz;

    alignas(64) TDoubleArray buffer_yy;

    alignas(64) TDoubleArray buffer_yz;

    alignas(64) TDoubleArray buffer_zz;

    // loop over integral batches

    const auto nbatches = batch::getNumberOfBatches(ket_ncgtos, simd_width);

    for (int64_t i = 0; i < nbatches; i++)
    {
        const auto [ket_first, ket_last] = batch::getBatchRange(i, ket_ncgtos, simd_width);

        const auto ket_dim = ket_last - ket_first;

        simd::loadCoordinates(ket_coords_x,
                              ket_coords_y,
                              ket_coords_z,
                              ket_gto_coords,
                              ket_first,
                              ket_last);

        for (int64_t j = bra_first; j < bra_last; j++) 
        {
            const auto bra_coord = bra_gto_coords[j];

            // compute primitive integrals block (XXX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_XXX(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_XXY(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_XXZ(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_XYY(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_XYZ(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_XZZ(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_YYY(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_YYZ(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_YZZ(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDF_T_ZZZ(buffer_xx,
                                                                   buffer_xy,
                                                                   buffer_xz,
                                                                   buffer_yy,
                                                                   buffer_yz,
                                                                   buffer_zz,
                                                                   bra_exp,
                                                                   bra_norm,
                                                                   bra_coord,
                                                                   ket_exps,
                                                                   ket_norms,
                                                                   ket_coords_x,
                                                                   ket_coords_y,
                                                                   ket_coords_z,
                                                                   ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

        }
    }
}

auto
compPrimitiveNuclearPotentialDF_T_XXX(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * (3.0 * fe_0 * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_x);

        fints_xx[i] += fss * bf_values[0][i] * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x;

        fints_xx[i] += fss * bf_values[1][i] * (-9.0 * fe_0 * rpa_x * rpb_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x);

        fints_xx[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x - 3.0 * fe_0 * fe_0 * rpa_x - (9.0 / 2.0) * fe_0 * fe_0 * rpb_x - (15.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xx[i] += fss * bf_values[1][i] * -3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x;

        fints_xx[i] += fss * bf_values[2][i] * (9.0 * fe_0 * rpa_x * rpb_x * rpc_x + 6.0 * fe_0 * rpa_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x + 9.0 * fe_0 * rpb_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpc_x + 6.0 * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x;

        fints_xx[i] += fss * bf_values[3][i] * (-6.0 * fe_0 * rpa_x * rpc_x * rpc_x - 9.0 * fe_0 * rpb_x * rpc_x * rpc_x - 5.0 * fe_0 * rpc_x * rpc_x * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpc_x - 6.0 * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-rpa_x * rpa_x * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * (5.0 * fe_0 * rpc_x * rpc_x * rpc_x + 2.0 * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y);

        fints_xy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 3.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x - rpa_y * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * -rpa_x * rpb_x * rpb_x * rpb_x * rpc_y;

        fints_xy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x + 3.0 * fe_0 * rpa_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * ((9.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y + 3.0 * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (3.0 * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x + rpb_x * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpc_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpa_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_y * rpc_x * rpc_x + rpa_y * rpc_x * rpc_x * rpc_x * rpc_x + rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[5][i] * -rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z);

        fints_xz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 3.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_x - rpa_z * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * -rpa_x * rpb_x * rpb_x * rpb_x * rpc_z;

        fints_xz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * ((9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 3.0 * rpa_z * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (3.0 * rpa_z * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + rpb_x * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpc_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_x * rpc_x + rpa_z * rpc_x * rpc_x * rpc_x * rpc_x + rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpa_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x + 6.0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y;

        fints_yy[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 6.0 * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x - rpa_y * rpa_y * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x + 2.0 * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_yz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z - 3.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_yz[i] += fss * bf_values[1][i] * (-rpa_z * rpb_x * rpb_x * rpb_x * rpc_y - rpa_y * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y + 3.0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x + 3.0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x + rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x - rpa_z * rpa_y * rpc_x * rpc_x * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-3.0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_zz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpa_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpb_x * rpb_x * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x + 6.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z;

        fints_zz[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 6.0 * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x - rpa_z * rpa_z * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_XXY(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * (2.0 * fe_0 * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xx[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_x * rpb_y * rpb_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y);

        fints_xx[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y);

        fints_xx[i] += fss * bf_values[1][i] * (-2.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x - rpa_x * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_x * rpb_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_x * rpc_y + 3.0 * fe_0 * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y + 3.0 * fe_0 * rpb_y * rpb_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * (4.0 * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x + rpa_x * rpa_x * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x + rpb_y * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_x * rpc_y * rpc_x - 3.0 * fe_0 * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpc_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y);

        fints_xx[i] += fss * bf_values[3][i] * (-2.0 * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x - 4.0 * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - rpa_x * rpa_x * rpc_y * rpc_x * rpc_x - 2.0 * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x - rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y + fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x);

        fints_xy[i] += fss * bf_values[0][i] * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y - fe_0 * rpa_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x - fe_0 * rpa_y * rpb_x * rpc_y);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y - fe_0 * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x - fe_0 * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x - fe_0 * fe_0 * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x - rpa_y * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xy[i] += fss * bf_values[1][i] * (-rpa_y * rpb_y * rpb_x * rpb_x * rpc_x - rpa_x * rpb_y * rpb_x * rpb_x * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x + fe_0 * rpa_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * (fe_0 * rpa_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x + fe_0 * rpb_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (fe_0 * rpb_x * rpc_y * rpc_y + fe_0 * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x);

        fints_xy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x + rpa_y * rpa_x * rpb_y * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x + rpa_y * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (2.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + rpa_x * rpb_x * rpb_x * rpc_y * rpc_y + rpb_y * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x - fe_0 * rpb_x * rpc_y * rpc_y);

        fints_xy[i] += fss * bf_values[3][i] * (-fe_0 * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - rpa_y * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpb_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x - rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * -rpb_x * rpb_x * rpc_y * rpc_y * rpc_x;

        fints_xy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x + rpa_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * 2.0 * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y + fe_0 * rpa_z * rpb_y * rpb_x + rpa_z * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y - fe_0 * rpa_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x - fe_0 * rpa_z * rpb_x * rpc_y);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z - fe_0 * rpb_y * rpb_x * rpc_z - 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_x - rpa_z * rpa_x * rpb_x * rpb_x * rpc_y - rpa_z * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * -rpa_x * rpb_y * rpb_x * rpb_x * rpc_z;

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x + fe_0 * rpa_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y + fe_0 * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x + fe_0 * rpb_x * rpc_z * rpc_y + rpa_z * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (2.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_z * rpb_y * rpb_x * rpc_x * rpc_x + rpa_z * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xz[i] += fss * bf_values[2][i] * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x;

        fints_xz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x - fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_x * rpc_y * rpc_x * rpc_x - rpa_z * rpb_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x - rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-2.0 * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + 2.0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * (fe_0 * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += fss * bf_values[0][i] * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x;

        fints_yy[i] += fss * bf_values[1][i] * (-fe_0 * rpa_y * rpb_y * rpc_y - 2.0 * fe_0 * rpa_y * rpb_x * rpc_x - fe_0 * rpa_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yy[i] += fss * bf_values[1][i] * (-fe_0 * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y - fe_0 * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += fss * bf_values[1][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_x - rpa_y * rpa_y * rpb_x * rpb_x * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (fe_0 * rpa_y * rpb_y * rpc_y + 2.0 * fe_0 * rpa_y * rpb_x * rpc_x + fe_0 * rpa_y * rpc_y * rpc_y + fe_0 * rpa_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (fe_0 * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y + 4.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (rpa_y * rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x + rpb_y * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-fe_0 * rpa_y * rpc_y * rpc_y - fe_0 * rpa_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - 4.0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-rpa_y * rpa_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - rpb_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y + 2.0 * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x + rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y - fe_0 * rpa_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x);

        fints_yz[i] += fss * bf_values[1][i] * (-rpa_z * rpa_y * rpb_x * rpb_x * rpc_y - rpa_z * rpb_y * rpb_x * rpb_x * rpc_y - rpa_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y + fe_0 * rpa_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y + fe_0 * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * fe_0 * rpc_z + rpa_z * rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x + 2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * (rpa_z * rpb_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x + rpa_y * rpb_x * rpb_x * rpc_z * rpc_y + rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y - fe_0 * rpb_x * rpc_z * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - rpa_z * rpa_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-2.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x - rpa_y * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + rpa_z * rpc_y * rpc_y * rpc_x * rpc_x + rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yz[i] += fss * bf_values[4][i] * 2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_zz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpb_y * rpb_x * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-fe_0 * rpa_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y - fe_0 * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_x);

        fints_zz[i] += fss * bf_values[1][i] * -rpa_z * rpa_z * rpb_x * rpb_x * rpc_y;

        fints_zz[i] += fss * bf_values[2][i] * (fe_0 * rpa_z * rpb_y * rpc_z + fe_0 * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y + fe_0 * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x + fe_0 * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * (4.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_z * rpb_x * rpc_y * rpc_x + rpb_y * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[3][i] * (-fe_0 * rpa_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x - fe_0 * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_z * rpb_y * rpc_z * rpc_x * rpc_x - 4.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpc_y * rpc_x * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x - rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x + rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_XXZ(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * (2.0 * fe_0 * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xx[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_x * rpb_z * rpb_x - 3.0 * fe_0 * rpa_x * rpb_z * rpc_x - 2.0 * fe_0 * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z);

        fints_xx[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z);

        fints_xx[i] += fss * bf_values[1][i] * (-2.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x - rpa_x * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xx[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_x * rpb_z * rpc_x + 2.0 * fe_0 * rpa_x * rpb_x * rpc_z + 3.0 * fe_0 * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z + 3.0 * fe_0 * rpb_z * rpb_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z);

        fints_xx[i] += fss * bf_values[2][i] * (4.0 * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_z * rpc_x * rpc_x + 2.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x + rpb_z * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_x * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpc_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z);

        fints_xx[i] += fss * bf_values[3][i] * (-2.0 * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x - 4.0 * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - rpa_x * rpa_x * rpc_z * rpc_x * rpc_x - 2.0 * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x - rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + rpb_z * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_z * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z + fe_0 * rpa_y * rpb_z * rpb_x + rpa_y * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z - fe_0 * rpa_y * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x - fe_0 * rpa_y * rpb_x * rpc_z);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y - fe_0 * rpb_z * rpb_x * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x - rpa_y * rpa_x * rpb_x * rpb_x * rpc_z - rpa_y * rpb_z * rpb_x * rpb_x * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * -rpa_x * rpb_z * rpb_x * rpb_x * rpc_y;

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x + fe_0 * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y + fe_0 * rpb_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x + fe_0 * rpb_x * rpc_z * rpc_y + rpa_y * rpa_x * rpb_z * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x + rpa_y * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + rpa_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x;

        fints_xy[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x - fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpa_x * rpc_z * rpc_x * rpc_x - rpa_y * rpb_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x - rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_x * rpc_x * rpc_x + rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z + fe_0 * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x);

        fints_xz[i] += fss * bf_values[0][i] * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z - fe_0 * rpa_z * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x - fe_0 * rpa_z * rpb_x * rpc_z);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z - fe_0 * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x - fe_0 * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x - fe_0 * fe_0 * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_x - rpa_z * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xz[i] += fss * bf_values[1][i] * (-rpa_z * rpb_z * rpb_x * rpb_x * rpc_x - rpa_x * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x + fe_0 * rpa_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * (fe_0 * rpa_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x + fe_0 * rpb_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (fe_0 * rpb_x * rpc_z * rpc_z + fe_0 * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x);

        fints_xz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x + rpa_z * rpa_x * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_z * rpb_x * rpc_x * rpc_x + rpa_z * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + rpa_x * rpb_x * rpb_x * rpc_z * rpc_z + rpb_z * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x - fe_0 * rpb_x * rpc_z * rpc_z);

        fints_xz[i] += fss * bf_values[3][i] * (-fe_0 * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - rpa_z * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_z * rpb_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x - rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * -rpb_x * rpb_x * rpc_z * rpc_z * rpc_x;

        fints_xz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x + rpa_z * rpc_z * rpc_x * rpc_x * rpc_x + rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * 2.0 * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-fe_0 * rpa_y * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z - fe_0 * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_x);

        fints_yy[i] += fss * bf_values[1][i] * -rpa_y * rpa_y * rpb_x * rpb_x * rpc_z;

        fints_yy[i] += fss * bf_values[2][i] * (fe_0 * rpa_y * rpb_z * rpc_y + fe_0 * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z + fe_0 * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x + fe_0 * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpc_z);

        fints_yy[i] += fss * bf_values[2][i] * (4.0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpb_z * rpc_x * rpc_x + 2.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x + rpb_z * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-fe_0 * rpa_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x - fe_0 * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpc_z * rpc_x * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + 2.0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z - fe_0 * rpa_y * rpb_x * rpc_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x);

        fints_yz[i] += fss * bf_values[1][i] * (-rpa_z * rpa_y * rpb_x * rpb_x * rpc_z - rpa_z * rpb_z * rpb_x * rpb_x * rpc_y - rpa_y * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpb_x * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y + fe_0 * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpc_y + rpa_z * rpa_y * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * (rpa_z * rpb_x * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x + rpa_y * rpb_x * rpb_x * rpc_z * rpc_z + rpb_z * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y - fe_0 * rpb_x * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - rpa_z * rpa_y * rpc_z * rpc_x * rpc_x - rpa_z * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + rpa_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yz[i] += fss * bf_values[4][i] * 2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_zz[i] += fss * bf_values[0][i] * (fe_0 * rpa_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z);

        fints_zz[i] += fss * bf_values[0][i] * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x;

        fints_zz[i] += fss * bf_values[1][i] * (-fe_0 * rpa_z * rpb_z * rpc_z - 2.0 * fe_0 * rpa_z * rpb_x * rpc_x - fe_0 * rpa_z * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z);

        fints_zz[i] += fss * bf_values[1][i] * (-fe_0 * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z - fe_0 * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z);

        fints_zz[i] += fss * bf_values[1][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_x - rpa_z * rpa_z * rpb_x * rpb_x * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (fe_0 * rpa_z * rpb_z * rpc_z + 2.0 * fe_0 * rpa_z * rpb_x * rpc_x + fe_0 * rpa_z * rpc_z * rpc_z + fe_0 * rpa_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (fe_0 * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 4.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (rpa_z * rpa_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_x + rpb_z * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[3][i] * (-fe_0 * rpa_z * rpc_z * rpc_z - fe_0 * rpa_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x - 4.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_z * rpc_z * rpc_x * rpc_x - 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x - rpb_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z + 2.0 * rpa_z * rpc_z * rpc_z * rpc_x * rpc_x + rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_XYY(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * (fe_0 * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x);

        fints_xx[i] += fss * bf_values[0][i] * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x;

        fints_xx[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_x * rpb_y * rpc_y - fe_0 * rpa_x * rpb_y * rpb_y - fe_0 * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x);

        fints_xx[i] += fss * bf_values[1][i] * (-fe_0 * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x - fe_0 * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x);

        fints_xx[i] += fss * bf_values[1][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y - rpa_x * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (2.0 * fe_0 * rpa_x * rpb_y * rpc_y + fe_0 * rpa_x * rpb_x * rpc_x + fe_0 * rpa_x * rpc_y * rpc_y + fe_0 * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (fe_0 * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x + 4.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + 2.0 * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (2.0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x + rpa_x * rpa_x * rpb_x * rpc_y * rpc_y + rpb_y * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-fe_0 * rpa_x * rpc_y * rpc_y - fe_0 * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 4.0 * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-rpa_x * rpa_x * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x + 2.0 * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y + fe_0 * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y);

        fints_xy[i] += fss * bf_values[0][i] * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x;

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x - fe_0 * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * (-fe_0 * rpa_x * rpb_y * rpb_x - fe_0 * rpa_x * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y - fe_0 * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_y - fe_0 * fe_0 * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y - rpa_y * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * (-rpa_y * rpb_y * rpb_y * rpb_x * rpc_x - rpa_x * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x + fe_0 * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (fe_0 * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x + fe_0 * rpb_y * rpb_x * rpc_x + fe_0 * rpb_y * rpc_y * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * (fe_0 * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y);

        fints_xy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_y + 2.0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x + rpa_y * rpa_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x + rpa_y * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (2.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y + rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + rpb_y * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x - fe_0 * rpb_y * rpc_y * rpc_y - fe_0 * rpb_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - rpa_y * rpa_x * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-2.0 * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - rpa_y * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - rpa_x * rpb_x * rpc_y * rpc_y * rpc_y - 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * -rpb_y * rpb_y * rpc_y * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y + rpa_y * rpc_y * rpc_y * rpc_x * rpc_x + rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xy[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x - fe_0 * rpa_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y);

        fints_xz[i] += fss * bf_values[1][i] * (-rpa_z * rpa_x * rpb_y * rpb_y * rpc_x - rpa_z * rpb_y * rpb_y * rpb_x * rpc_x - rpa_x * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x + fe_0 * rpa_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x + fe_0 * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * fe_0 * rpc_z + 2.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x + rpa_z * rpa_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (rpa_z * rpb_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + rpb_y * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x - fe_0 * rpb_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - rpa_z * rpa_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_z * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + rpa_z * rpc_y * rpc_y * rpc_x * rpc_x + rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * (2.0 * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpb_y * rpb_y * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_y * rpb_y * rpb_x - 2.0 * fe_0 * rpa_y * rpb_y * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x);

        fints_yy[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x);

        fints_yy[i] += fss * bf_values[1][i] * (-2.0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y - rpa_y * rpa_y * rpb_y * rpb_y * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * (2.0 * fe_0 * rpa_y * rpb_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x + 3.0 * fe_0 * rpb_y * rpb_x * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x + 3.0 * fe_0 * rpb_x * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * (4.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x + 2.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpb_x * rpc_y * rpc_y + rpb_y * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_x * rpc_y * rpc_y - 3.0 * fe_0 * rpc_y * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-4.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y - rpa_y * rpa_y * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpb_x * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x + fe_0 * rpa_z * rpb_y * rpb_x + rpa_z * rpa_y * rpb_y * rpb_y * rpb_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x - fe_0 * rpa_z * rpb_y * rpb_x - fe_0 * rpa_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z - fe_0 * rpb_y * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y - rpa_z * rpa_y * rpb_y * rpb_y * rpc_x - rpa_z * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * -rpa_y * rpb_y * rpb_y * rpb_x * rpc_z;

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x + fe_0 * rpa_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x + fe_0 * rpb_y * rpb_x * rpc_z + fe_0 * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * (rpa_z * rpa_y * rpb_x * rpc_y * rpc_y + 2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y + rpa_z * rpb_y * rpb_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y;

        fints_yz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x - fe_0 * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x - rpa_z * rpb_x * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_zz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpb_y * rpb_y * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-fe_0 * rpa_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x - fe_0 * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_zz[i] += fss * bf_values[1][i] * -rpa_z * rpa_z * rpb_y * rpb_y * rpc_x;

        fints_zz[i] += fss * bf_values[2][i] * (fe_0 * rpa_z * rpb_x * rpc_z + fe_0 * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x + fe_0 * rpb_y * rpb_x * rpc_y + fe_0 * rpb_y * rpc_y * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * (4.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_z * rpa_z * rpb_y * rpc_y * rpc_x + rpa_z * rpa_z * rpb_x * rpc_y * rpc_y + rpb_y * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[3][i] * (-fe_0 * rpa_z * rpc_z * rpc_x - fe_0 * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - 4.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y - rpa_z * rpa_z * rpc_y * rpc_y * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + 2.0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_XYZ(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * (fe_0 * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x + rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xx[i] += fss * bf_values[1][i] * (-fe_0 * rpa_x * rpb_z * rpb_y - fe_0 * rpa_x * rpb_z * rpc_y - fe_0 * rpa_x * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x);

        fints_xx[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z - 2.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x - rpa_x * rpa_x * rpb_z * rpb_y * rpc_x - rpa_x * rpa_x * rpb_z * rpb_x * rpc_y);

        fints_xx[i] += fss * bf_values[1][i] * -rpa_x * rpa_x * rpb_y * rpb_x * rpc_z;

        fints_xx[i] += fss * bf_values[2][i] * (fe_0 * rpa_x * rpb_z * rpc_y + fe_0 * rpa_x * rpb_y * rpc_z + fe_0 * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y + 2.0 * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (2.0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + rpa_x * rpa_x * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x;

        fints_xx[i] += fss * bf_values[3][i] * (-fe_0 * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-2.0 * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - rpb_z * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_y * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z);

        fints_xy[i] += fss * bf_values[1][i] * (-rpa_y * rpa_x * rpb_z * rpb_y * rpc_x - rpa_y * rpa_x * rpb_z * rpb_x * rpc_y - rpa_y * rpa_x * rpb_y * rpb_x * rpc_z - rpa_y * rpb_z * rpb_y * rpb_x * rpc_x - rpa_x * rpb_z * rpb_y * rpb_x * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpc_z + rpa_y * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (rpa_y * rpa_x * rpb_y * rpc_z * rpc_x + rpa_y * rpa_x * rpb_x * rpc_z * rpc_y + rpa_y * rpb_z * rpb_y * rpc_x * rpc_x + rpa_y * rpb_z * rpb_x * rpc_y * rpc_x + rpa_y * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (rpa_x * rpb_z * rpb_y * rpc_y * rpc_x + rpa_x * rpb_z * rpb_x * rpc_y * rpc_y + rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + rpb_z * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y);

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - rpa_y * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - rpa_y * rpb_y * rpc_z * rpc_x * rpc_x - rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - rpa_x * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * (rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_z * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y);

        fints_xz[i] += fss * bf_values[1][i] * (-rpa_z * rpa_x * rpb_z * rpb_y * rpc_x - rpa_z * rpa_x * rpb_z * rpb_x * rpc_y - rpa_z * rpa_x * rpb_y * rpb_x * rpc_z - rpa_z * rpb_z * rpb_y * rpb_x * rpc_x - rpa_x * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpc_y + rpa_z * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (rpa_z * rpa_x * rpb_y * rpc_z * rpc_x + rpa_z * rpa_x * rpb_x * rpc_z * rpc_y + rpa_z * rpb_z * rpb_y * rpc_x * rpc_x + rpa_z * rpb_z * rpb_x * rpc_y * rpc_x + rpa_z * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + rpa_x * rpb_y * rpb_x * rpc_z * rpc_z + rpb_z * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - rpa_z * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_z * rpb_z * rpc_y * rpc_x * rpc_x - rpa_z * rpb_y * rpc_z * rpc_x * rpc_x - rpa_z * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - rpa_x * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + rpa_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * (rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * (fe_0 * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x + rpa_y * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-fe_0 * rpa_y * rpb_z * rpb_x - fe_0 * rpa_y * rpb_z * rpc_x - fe_0 * rpa_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x);

        fints_yy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z - 2.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y - rpa_y * rpa_y * rpb_z * rpb_y * rpc_x - rpa_y * rpa_y * rpb_z * rpb_x * rpc_y);

        fints_yy[i] += fss * bf_values[1][i] * -rpa_y * rpa_y * rpb_y * rpb_x * rpc_z;

        fints_yy[i] += fss * bf_values[2][i] * (fe_0 * rpa_y * rpb_z * rpc_x + fe_0 * rpa_y * rpb_x * rpc_z + fe_0 * rpa_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * (2.0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpb_z * rpc_y * rpc_x + rpa_y * rpa_y * rpb_y * rpc_z * rpc_x + rpa_y * rpa_y * rpb_x * rpc_z * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y;

        fints_yy[i] += fss * bf_values[3][i] * (-fe_0 * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-2.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y - rpa_y * rpa_y * rpc_z * rpc_y * rpc_x - rpb_z * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-rpb_z * rpb_x * rpc_y * rpc_y * rpc_y - rpb_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + 2.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_z * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x);

        fints_yz[i] += fss * bf_values[1][i] * (-rpa_z * rpa_y * rpb_z * rpb_y * rpc_x - rpa_z * rpa_y * rpb_z * rpb_x * rpc_y - rpa_z * rpa_y * rpb_y * rpb_x * rpc_z - rpa_z * rpb_z * rpb_y * rpb_x * rpc_y - rpa_y * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpc_x + rpa_z * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * (rpa_z * rpa_y * rpb_y * rpc_z * rpc_x + rpa_z * rpa_y * rpb_x * rpc_z * rpc_y + rpa_z * rpb_z * rpb_y * rpc_y * rpc_x + rpa_z * rpb_z * rpb_x * rpc_y * rpc_y + rpa_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * (rpa_y * rpb_z * rpb_y * rpc_z * rpc_x + rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + rpa_y * rpb_y * rpb_x * rpc_z * rpc_z + rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z);

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - rpa_z * rpa_y * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-rpa_z * rpb_z * rpc_y * rpc_y * rpc_x - rpa_z * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpb_x * rpc_z * rpc_y * rpc_y - rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - rpa_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-rpa_y * rpb_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - rpb_z * rpb_x * rpc_z * rpc_y * rpc_y - rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + rpa_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[4][i] * (rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_zz[i] += fss * bf_values[0][i] * (fe_0 * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x + rpa_z * rpa_z * rpb_z * rpb_y * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-fe_0 * rpa_z * rpb_y * rpb_x - fe_0 * rpa_z * rpb_y * rpc_x - fe_0 * rpa_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x);

        fints_zz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z - 2.0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z - rpa_z * rpa_z * rpb_z * rpb_y * rpc_x - rpa_z * rpa_z * rpb_z * rpb_x * rpc_y);

        fints_zz[i] += fss * bf_values[1][i] * -rpa_z * rpa_z * rpb_y * rpb_x * rpc_z;

        fints_zz[i] += fss * bf_values[2][i] * (fe_0 * rpa_z * rpb_y * rpc_x + fe_0 * rpa_z * rpb_x * rpc_y + fe_0 * rpa_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * (2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z + rpa_z * rpa_z * rpb_z * rpc_y * rpc_x + rpa_z * rpa_z * rpb_y * rpc_z * rpc_x + rpa_z * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z;

        fints_zz[i] += fss * bf_values[3][i] * (-fe_0 * rpa_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y - rpa_z * rpa_z * rpc_z * rpc_y * rpc_x - rpb_z * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - rpb_y * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_XZZ(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * (fe_0 * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x);

        fints_xx[i] += fss * bf_values[0][i] * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x;

        fints_xx[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_x * rpb_z * rpc_z - fe_0 * rpa_x * rpb_z * rpb_z - fe_0 * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x);

        fints_xx[i] += fss * bf_values[1][i] * (-fe_0 * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x - fe_0 * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x);

        fints_xx[i] += fss * bf_values[1][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z - rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (2.0 * fe_0 * rpa_x * rpb_z * rpc_z + fe_0 * rpa_x * rpb_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z + fe_0 * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (fe_0 * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x + 4.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + 2.0 * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * (2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x + rpa_x * rpa_x * rpb_x * rpc_z * rpc_z + rpb_z * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-fe_0 * rpa_x * rpc_z * rpc_z - fe_0 * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 4.0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - 2.0 * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - rpb_z * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_y * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x - fe_0 * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z);

        fints_xy[i] += fss * bf_values[1][i] * (-rpa_y * rpa_x * rpb_z * rpb_z * rpc_x - rpa_y * rpb_z * rpb_z * rpb_x * rpc_x - rpa_x * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x + fe_0 * rpa_y * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x + fe_0 * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpc_y + 2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x + rpa_y * rpa_x * rpb_x * rpc_z * rpc_z + 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (rpa_y * rpb_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + rpb_z * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x - fe_0 * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - rpa_y * rpa_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - rpb_z * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_xy[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z + fe_0 * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z);

        fints_xz[i] += fss * bf_values[0][i] * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x;

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * (-fe_0 * rpa_x * rpb_z * rpb_x - fe_0 * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z - fe_0 * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z - fe_0 * fe_0 * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z - rpa_z * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * (-rpa_z * rpb_z * rpb_z * rpb_x * rpc_x - rpa_x * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x + fe_0 * rpa_z * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (fe_0 * rpa_x * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x + fe_0 * rpb_z * rpb_x * rpc_x + fe_0 * rpb_z * rpc_z * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * (fe_0 * rpb_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z);

        fints_xz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 2.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_x + rpa_z * rpa_x * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x + rpa_z * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z + rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + rpb_z * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x - fe_0 * rpb_z * rpc_z * rpc_z - fe_0 * rpb_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - rpa_z * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-2.0 * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x - rpa_z * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - rpa_x * rpb_x * rpc_z * rpc_z * rpc_z - 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * -rpb_z * rpb_z * rpc_z * rpc_x * rpc_x;

        fints_xz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z + rpa_z * rpc_z * rpc_z * rpc_x * rpc_x + rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpb_z * rpb_z * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-fe_0 * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x);

        fints_yy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z);

        fints_yy[i] += fss * bf_values[1][i] * -rpa_y * rpa_y * rpb_z * rpb_z * rpc_x;

        fints_yy[i] += fss * bf_values[2][i] * (fe_0 * rpa_y * rpb_x * rpc_y + fe_0 * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x + fe_0 * rpb_z * rpb_x * rpc_z + fe_0 * rpb_z * rpc_z * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpc_x);

        fints_yy[i] += fss * bf_values[2][i] * (4.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x + rpa_y * rpa_y * rpb_x * rpc_z * rpc_z + rpb_z * rpb_z * rpb_x * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-fe_0 * rpa_y * rpc_y * rpc_x - fe_0 * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - 4.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpc_z * rpc_z * rpc_x);

        fints_yy[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_yy[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_yz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x + fe_0 * rpa_y * rpb_z * rpb_x + rpa_z * rpa_y * rpb_z * rpb_z * rpb_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y - fe_0 * rpa_y * rpb_z * rpb_x - fe_0 * rpa_y * rpb_z * rpc_x);

        fints_yz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z - fe_0 * rpb_z * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z - rpa_z * rpa_y * rpb_z * rpb_z * rpc_x - rpa_z * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * -rpa_y * rpb_z * rpb_z * rpb_x * rpc_z;

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x + fe_0 * rpa_y * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x + fe_0 * rpb_z * rpb_x * rpc_y + fe_0 * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * (rpa_z * rpa_y * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y + rpa_z * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z + rpa_y * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_yz[i] += fss * bf_values[2][i] * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y;

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x - fe_0 * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - rpa_y * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_yz[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_zz[i] += fss * bf_values[0][i] * (2.0 * fe_0 * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpb_z * rpb_z * rpb_x);

        fints_zz[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_z * rpb_z * rpb_x - 2.0 * fe_0 * rpa_z * rpb_z * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x);

        fints_zz[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x);

        fints_zz[i] += fss * bf_values[1][i] * (-2.0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_z - rpa_z * rpa_z * rpb_z * rpb_z * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * (2.0 * fe_0 * rpa_z * rpb_z * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x + 3.0 * fe_0 * rpb_z * rpb_x * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x);

        fints_zz[i] += fss * bf_values[2][i] * (4.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x + 2.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_x + rpa_z * rpa_z * rpb_x * rpc_z * rpc_z + rpb_z * rpb_z * rpb_x * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_z - 3.0 * fe_0 * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x);

        fints_zz[i] += fss * bf_values[3][i] * (-4.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z - rpa_z * rpa_z * rpc_z * rpc_z * rpc_x - 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_zz[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_x + 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpb_x * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_z * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_YYY(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_xx[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpa_x * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y);

        fints_xx[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_x * rpb_y * rpc_x + 3.0 * fe_0 * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y + 6.0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x;

        fints_xx[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y);

        fints_xx[i] += fss * bf_values[3][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 6.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y + 2.0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_y * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_xy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y);

        fints_xy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 3.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y - rpa_y * rpb_y * rpb_y * rpb_y * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * -rpa_x * rpb_y * rpb_y * rpb_y * rpc_y;

        fints_xy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y + 3.0 * fe_0 * rpa_x * rpc_y * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((9.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x + 3.0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * (3.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y + rpb_y * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpc_y * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpa_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_y * rpc_y * rpc_x + rpa_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y + rpa_z * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_xz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z - 3.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_xz[i] += fss * bf_values[1][i] * (-rpa_z * rpb_y * rpb_y * rpb_y * rpc_x - rpa_x * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y);

        fints_xz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x + 3.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y + 3.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y + rpb_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_xz[i] += fss * bf_values[3][i] * (-3.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * (3.0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += fss * bf_values[0][i] * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y;

        fints_yy[i] += fss * bf_values[1][i] * (-9.0 * fe_0 * rpa_y * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_yy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y - 3.0 * fe_0 * fe_0 * rpa_y - (9.0 / 2.0) * fe_0 * fe_0 * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yy[i] += fss * bf_values[1][i] * -3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y;

        fints_yy[i] += fss * bf_values[2][i] * (9.0 * fe_0 * rpa_y * rpb_y * rpc_y + 6.0 * fe_0 * rpa_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y + 9.0 * fe_0 * rpb_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpc_y + 6.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y;

        fints_yy[i] += fss * bf_values[3][i] * (-6.0 * fe_0 * rpa_y * rpc_y * rpc_y - 9.0 * fe_0 * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpc_y * rpc_y * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[4][i] * (5.0 * fe_0 * rpc_y * rpc_y * rpc_y + 2.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_y * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z);

        fints_yz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y - rpa_z * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * -rpa_y * rpb_y * rpb_y * rpb_y * rpc_z;

        fints_yz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y + 3.0 * fe_0 * rpa_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * ((9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 3.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * (3.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y + rpb_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z);

        fints_yz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_y * rpc_y + rpa_z * rpc_y * rpc_y * rpc_y * rpc_y + rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_zz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpb_y * rpb_y * rpb_y);

        fints_zz[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpa_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y);

        fints_zz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_z * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y + 6.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpb_y * rpc_y * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z;

        fints_zz[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y);

        fints_zz[i] += fss * bf_values[3][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 6.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y - rpa_z * rpa_z * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y + 2.0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_YYZ(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_xx[i] += fss * bf_values[1][i] * (-fe_0 * rpa_x * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z - fe_0 * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y);

        fints_xx[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y);

        fints_xx[i] += fss * bf_values[1][i] * -rpa_x * rpa_x * rpb_y * rpb_y * rpc_z;

        fints_xx[i] += fss * bf_values[2][i] * (fe_0 * rpa_x * rpb_z * rpc_x + fe_0 * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z + fe_0 * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x + fe_0 * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpc_z);

        fints_xx[i] += fss * bf_values[2][i] * (4.0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x + 2.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpb_z * rpc_y * rpc_y + 2.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y + rpb_z * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-fe_0 * rpa_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x - fe_0 * rpb_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y);

        fints_xx[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - 4.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_xx[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z + fe_0 * rpa_x * rpb_z * rpb_y + rpa_y * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x - fe_0 * rpa_x * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y);

        fints_xy[i] += fss * bf_values[1][i] * (-fe_0 * rpa_x * rpb_y * rpc_z - fe_0 * rpb_z * rpb_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y - rpa_y * rpa_x * rpb_y * rpb_y * rpc_z - rpa_y * rpb_z * rpb_y * rpb_y * rpc_x);

        fints_xy[i] += fss * bf_values[1][i] * -rpa_x * rpb_z * rpb_y * rpb_y * rpc_y;

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y + fe_0 * rpa_x * rpb_y * rpc_z);

        fints_xy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x + fe_0 * rpb_y * rpc_z * rpc_x + rpa_y * rpa_x * rpb_z * rpc_y * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * (2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x + rpa_y * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y + rpa_x * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x;

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x - fe_0 * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpa_x * rpc_z * rpc_y * rpc_y - rpa_y * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpb_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_xy[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x - rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_z * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z - fe_0 * rpa_x * rpb_y * rpc_y);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y);

        fints_xz[i] += fss * bf_values[1][i] * (-rpa_z * rpa_x * rpb_y * rpb_y * rpc_z - rpa_z * rpb_z * rpb_y * rpb_y * rpc_x - rpa_x * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z + fe_0 * rpa_x * rpb_y * rpc_y);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x + fe_0 * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpc_x + rpa_z * rpa_x * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[2][i] * (rpa_z * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + rpa_x * rpb_y * rpb_y * rpc_z * rpc_z + rpb_z * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x - fe_0 * rpb_y * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - rpa_z * rpa_x * rpc_z * rpc_y * rpc_y - rpa_z * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + rpa_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * (2.0 * fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yy[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_y * rpb_z * rpb_y - 3.0 * fe_0 * rpa_y * rpb_z * rpc_y - 2.0 * fe_0 * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z);

        fints_yy[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z);

        fints_yy[i] += fss * bf_values[1][i] * (-2.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - rpa_y * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_yy[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_y * rpb_z * rpc_y + 2.0 * fe_0 * rpa_y * rpb_y * rpc_z + 3.0 * fe_0 * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z + 3.0 * fe_0 * rpb_z * rpb_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpb_z * rpc_y * rpc_y + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z);

        fints_yy[i] += fss * bf_values[2][i] * (4.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y + rpb_z * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_y * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z);

        fints_yy[i] += fss * bf_values[3][i] * (-2.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y - 4.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - 2.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_y * rpc_y + 2.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z + fe_0 * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y);

        fints_yz[i] += fss * bf_values[0][i] * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y;

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z - fe_0 * rpa_z * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y - fe_0 * rpa_z * rpb_y * rpc_z);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z - fe_0 * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y - fe_0 * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_y - fe_0 * fe_0 * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y - rpa_z * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_yz[i] += fss * bf_values[1][i] * (-rpa_z * rpb_z * rpb_y * rpb_y * rpc_y - rpa_y * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y + fe_0 * rpa_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * (fe_0 * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y + fe_0 * rpb_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * (fe_0 * rpb_y * rpc_z * rpc_z + fe_0 * rpb_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y);

        fints_yz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_y + rpa_z * rpa_y * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y + rpa_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * (2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpb_y * rpb_y * rpc_z * rpc_z + rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y - fe_0 * rpb_y * rpc_z * rpc_z);

        fints_yz[i] += fss * bf_values[3][i] * (-fe_0 * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - rpa_z * rpa_y * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-rpa_z * rpb_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y - rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * -rpb_y * rpb_y * rpc_z * rpc_z * rpc_y;

        fints_yz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y + rpa_z * rpc_z * rpc_y * rpc_y * rpc_y + rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[4][i] * 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_zz[i] += fss * bf_values[0][i] * (fe_0 * rpa_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z);

        fints_zz[i] += fss * bf_values[0][i] * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y;

        fints_zz[i] += fss * bf_values[1][i] * (-fe_0 * rpa_z * rpb_z * rpc_z - 2.0 * fe_0 * rpa_z * rpb_y * rpc_y - fe_0 * rpa_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z);

        fints_zz[i] += fss * bf_values[1][i] * (-fe_0 * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z - fe_0 * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z);

        fints_zz[i] += fss * bf_values[1][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_y - rpa_z * rpa_z * rpb_y * rpb_y * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (fe_0 * rpa_z * rpb_z * rpc_z + 2.0 * fe_0 * rpa_z * rpb_y * rpc_y + fe_0 * rpa_z * rpc_z * rpc_z + fe_0 * rpa_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (fe_0 * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 4.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (rpa_z * rpa_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_y + rpb_z * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[3][i] * (-fe_0 * rpa_z * rpc_z * rpc_z - fe_0 * rpa_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y - 4.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_z * rpc_z * rpc_y * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - rpb_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z + 2.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_YZZ(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_xx[i] += fss * bf_values[1][i] * (-fe_0 * rpa_x * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y - fe_0 * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y);

        fints_xx[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_xx[i] += fss * bf_values[1][i] * -rpa_x * rpa_x * rpb_z * rpb_z * rpc_y;

        fints_xx[i] += fss * bf_values[2][i] * (fe_0 * rpa_x * rpb_y * rpc_x + fe_0 * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y + fe_0 * rpb_z * rpb_y * rpc_z + fe_0 * rpb_z * rpc_z * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpc_y);

        fints_xx[i] += fss * bf_values[2][i] * (4.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + 2.0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + rpa_x * rpa_x * rpb_y * rpc_z * rpc_z + rpb_z * rpb_z * rpb_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[3][i] * (-fe_0 * rpa_x * rpc_y * rpc_x - fe_0 * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y);

        fints_xx[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_y - 4.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_xx[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - rpb_z * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_y * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x - fe_0 * rpa_x * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z);

        fints_xy[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_xy[i] += fss * bf_values[1][i] * (-rpa_y * rpa_x * rpb_z * rpb_z * rpc_y - rpa_y * rpb_z * rpb_z * rpb_y * rpc_x - rpa_x * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x + fe_0 * rpa_x * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y + fe_0 * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpc_x + 2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y + rpa_y * rpa_x * rpb_y * rpc_z * rpc_z + 2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xy[i] += fss * bf_values[2][i] * (rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + rpa_x * rpb_z * rpb_z * rpc_y * rpc_y + rpb_z * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y - fe_0 * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpc_x - rpa_y * rpa_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-rpa_y * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y - rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - rpb_z * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * ((1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x + rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_xy[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y + fe_0 * rpa_x * rpb_z * rpb_y + rpa_z * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_xz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x - fe_0 * rpa_x * rpb_z * rpb_y - fe_0 * rpa_x * rpb_z * rpc_y);

        fints_xz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z - fe_0 * rpb_z * rpb_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z - rpa_z * rpa_x * rpb_z * rpb_z * rpc_y - rpa_z * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * -rpa_x * rpb_z * rpb_z * rpb_y * rpc_z;

        fints_xz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x + fe_0 * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_x + fe_0 * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_xz[i] += fss * bf_values[2][i] * (rpa_z * rpa_x * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x + rpa_z * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z + rpa_x * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_xz[i] += fss * bf_values[2][i] * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x;

        fints_xz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y - fe_0 * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - rpa_x * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_xz[i] += fss * bf_values[3][i] * (-2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x - rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpb_y * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * (fe_0 * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += fss * bf_values[0][i] * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y;

        fints_yy[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_y * rpb_z * rpc_z - fe_0 * rpa_y * rpb_z * rpb_z - fe_0 * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yy[i] += fss * bf_values[1][i] * (-fe_0 * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y - fe_0 * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y);

        fints_yy[i] += fss * bf_values[1][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 2.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - rpa_y * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (2.0 * fe_0 * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpb_y * rpc_y + fe_0 * rpa_y * rpc_z * rpc_z + fe_0 * rpa_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (fe_0 * rpb_z * rpb_y * rpc_z + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y + 4.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * (2.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + rpa_y * rpa_y * rpb_y * rpc_z * rpc_z + rpb_z * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-fe_0 * rpa_y * rpc_z * rpc_z - fe_0 * rpa_y * rpc_y * rpc_y - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 4.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yy[i] += fss * bf_values[3][i] * (-rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y + 2.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_yz[i] += fss * bf_values[0][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z + fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z);

        fints_yz[i] += fss * bf_values[0][i] * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y;

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y - fe_0 * rpa_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * (-fe_0 * rpa_y * rpb_z * rpb_y - fe_0 * rpa_y * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z - fe_0 * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z);

        fints_yz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z - fe_0 * fe_0 * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z - rpa_z * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * (-rpa_z * rpb_z * rpb_z * rpb_y * rpc_y - rpa_y * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y + fe_0 * rpa_z * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * (fe_0 * rpa_y * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_y + fe_0 * rpb_z * rpc_z * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * (fe_0 * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z);

        fints_yz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y + rpa_z * rpa_y * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y + rpa_z * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[2][i] * (2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z + rpa_y * rpb_z * rpb_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y - fe_0 * rpb_z * rpc_z * rpc_z - fe_0 * rpb_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - rpa_z * rpa_y * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y - rpa_z * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - rpa_y * rpb_y * rpc_z * rpc_z * rpc_z - 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * -rpb_z * rpb_z * rpc_z * rpc_y * rpc_y;

        fints_yz[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z + rpa_z * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yz[i] += fss * bf_values[4][i] * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_zz[i] += fss * bf_values[0][i] * (2.0 * fe_0 * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpb_z * rpb_z * rpb_y);

        fints_zz[i] += fss * bf_values[1][i] * (-2.0 * fe_0 * rpa_z * rpb_z * rpb_y - 2.0 * fe_0 * rpa_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y);

        fints_zz[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y);

        fints_zz[i] += fss * bf_values[1][i] * (-2.0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z - 2.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_z - rpa_z * rpa_z * rpb_z * rpb_z * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * (2.0 * fe_0 * rpa_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_z * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpb_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y);

        fints_zz[i] += fss * bf_values[2][i] * (4.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y + 2.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_y + rpa_z * rpa_z * rpb_y * rpc_z * rpc_z + rpb_z * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_z * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y);

        fints_zz[i] += fss * bf_values[3][i] * (-4.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z - rpa_z * rpa_z * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_zz[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_y + 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpb_y * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_z * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialDF_T_ZZZ(      TDoubleArray& buffer_xx,
                                            TDoubleArray& buffer_xy,
                                            TDoubleArray& buffer_xz,
                                            TDoubleArray& buffer_yy,
                                            TDoubleArray& buffer_yz,
                                            TDoubleArray& buffer_zz,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xx = buffer_xx.data();

    auto fints_xy = buffer_xy.data();

    auto fints_xz = buffer_xz.data();

    auto fints_yy = buffer_yy.data();

    auto fints_yz = buffer_yz.data();

    auto fints_zz = buffer_zz.data();

    #pragma omp simd aligned(fints_xx,\
                             fints_xy,\
                             fints_xz,\
                             fints_yy,\
                             fints_yz,\
                             fints_zz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xx[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_xx[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z);

        fints_xx[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_xx[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_x * rpb_z * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 6.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_xx[i] += fss * bf_values[2][i] * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x;

        fints_xx[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z);

        fints_xx[i] += fss * bf_values[3][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 6.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z + 2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xx[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_xy[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z + rpa_y * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_xy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y - 3.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_xy[i] += fss * bf_values[1][i] * (-rpa_y * rpb_z * rpb_z * rpb_z * rpc_x - rpa_x * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y);

        fints_xy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x + 3.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z + 3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + rpb_z * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_xy[i] += fss * bf_values[3][i] * (-3.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - 3.0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xy[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_xz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_z * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_xz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z);

        fints_xz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x - 3.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z - rpa_z * rpb_z * rpb_z * rpb_z * rpc_x);

        fints_xz[i] += fss * bf_values[1][i] * -rpa_x * rpb_z * rpb_z * rpb_z * rpc_z;

        fints_xz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_x * rpc_z * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * ((9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x + 3.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_xz[i] += fss * bf_values[2][i] * (3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + rpb_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x - 3.0 * fe_0 * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpc_x);

        fints_xz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_x * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_x + rpa_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_yy[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);

        fints_yy[i] += fss * bf_values[1][i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z);

        fints_yy[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z);

        fints_yy[i] += fss * bf_values[2][i] * (3.0 * fe_0 * rpa_y * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z + 6.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_yy[i] += fss * bf_values[2][i] * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y;

        fints_yy[i] += fss * bf_values[3][i] * (-3.0 * fe_0 * rpa_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z);

        fints_yy[i] += fss * bf_values[3][i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpc_z - 6.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[4][i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z + 2.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yy[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_yz[i] += fss * bf_values[0][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_y * rpb_z * rpb_z * rpb_z);

        fints_yz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_yz[i] += fss * bf_values[1][i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y - 3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z - rpa_z * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_yz[i] += fss * bf_values[1][i] * -rpa_y * rpb_z * rpb_z * rpb_z * rpc_z;

        fints_yz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z + 3.0 * fe_0 * rpa_y * rpc_z * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * ((9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y + 3.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_yz[i] += fss * bf_values[2][i] * (3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z + rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpc_y);

        fints_yz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_y * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[4][i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_y + rpa_z * rpc_z * rpc_z * rpc_z * rpc_y + rpa_y * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_zz[i] += fss * bf_values[0][i] * (3.0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z + (9.0 / 4.0) * fe_0 * fe_0 * rpb_z);

        fints_zz[i] += fss * bf_values[0][i] * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z;

        fints_zz[i] += fss * bf_values[1][i] * (-9.0 * fe_0 * rpa_z * rpb_z * rpc_z - 3.0 * fe_0 * rpa_z * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z - (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z);

        fints_zz[i] += fss * bf_values[1][i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z - 3.0 * fe_0 * fe_0 * rpa_z - (9.0 / 2.0) * fe_0 * fe_0 * rpb_z - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z - 2.0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zz[i] += fss * bf_values[1][i] * -3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpc_z;

        fints_zz[i] += fss * bf_values[2][i] * (9.0 * fe_0 * rpa_z * rpb_z * rpc_z + 6.0 * fe_0 * rpa_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z + 9.0 * fe_0 * rpb_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z + (9.0 / 4.0) * fe_0 * fe_0 * rpb_z + (15.0 / 2.0) * fe_0 * fe_0 * rpc_z + 6.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[2][i] * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z;

        fints_zz[i] += fss * bf_values[3][i] * (-6.0 * fe_0 * rpa_z * rpc_z * rpc_z - 9.0 * fe_0 * rpb_z * rpc_z * rpc_z - 5.0 * fe_0 * rpc_z * rpc_z * rpc_z - (15.0 / 4.0) * fe_0 * fe_0 * rpc_z - 6.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[3][i] * (-rpa_z * rpa_z * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[4][i] * (5.0 * fe_0 * rpc_z * rpc_z * rpc_z + 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zz[i] += fss * bf_values[5][i] * -rpc_z * rpc_z * rpc_z * rpc_z * rpc_z;

    }
}

} // npotrec namespace

