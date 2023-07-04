#include "QuadrupoleRecGD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveQuadrupoleGD_XXXX_XX.hpp"
#include "PrimitiveQuadrupoleGD_XXXX_XY.hpp"
#include "PrimitiveQuadrupoleGD_XXXX_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXX_YY.hpp"
#include "PrimitiveQuadrupoleGD_XXXX_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXX_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXY_XX.hpp"
#include "PrimitiveQuadrupoleGD_XXXY_XY.hpp"
#include "PrimitiveQuadrupoleGD_XXXY_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXY_YY.hpp"
#include "PrimitiveQuadrupoleGD_XXXY_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXY_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_XXXZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_XXXZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_XXXZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XXXZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XXYY_XX.hpp"
#include "PrimitiveQuadrupoleGD_XXYY_XY.hpp"
#include "PrimitiveQuadrupoleGD_XXYY_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XXYY_YY.hpp"
#include "PrimitiveQuadrupoleGD_XXYY_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XXYY_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XXYZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_XXYZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_XXYZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XXYZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_XXYZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XXYZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XXZZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_XXZZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_XXZZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XXZZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_XXZZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XXZZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XYYY_XX.hpp"
#include "PrimitiveQuadrupoleGD_XYYY_XY.hpp"
#include "PrimitiveQuadrupoleGD_XYYY_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XYYY_YY.hpp"
#include "PrimitiveQuadrupoleGD_XYYY_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XYYY_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XYYZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_XYYZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_XYYZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XYYZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_XYYZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XYYZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XYZZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_XYZZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_XYZZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XYZZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_XYZZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XYZZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_XZZZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_XZZZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_XZZZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_XZZZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_XZZZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_XZZZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_YYYY_XX.hpp"
#include "PrimitiveQuadrupoleGD_YYYY_XY.hpp"
#include "PrimitiveQuadrupoleGD_YYYY_XZ.hpp"
#include "PrimitiveQuadrupoleGD_YYYY_YY.hpp"
#include "PrimitiveQuadrupoleGD_YYYY_YZ.hpp"
#include "PrimitiveQuadrupoleGD_YYYY_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_YYYZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_YYYZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_YYYZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_YYYZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_YYYZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_YYYZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_YYZZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_YYZZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_YYZZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_YYZZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_YYZZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_YYZZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_YZZZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_YZZZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_YZZZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_YZZZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_YZZZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_YZZZ_ZZ.hpp"
#include "PrimitiveQuadrupoleGD_ZZZZ_XX.hpp"
#include "PrimitiveQuadrupoleGD_ZZZZ_XY.hpp"
#include "PrimitiveQuadrupoleGD_ZZZZ_XZ.hpp"
#include "PrimitiveQuadrupoleGD_ZZZZ_YY.hpp"
#include "PrimitiveQuadrupoleGD_ZZZZ_YZ.hpp"
#include "PrimitiveQuadrupoleGD_ZZZZ_ZZ.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compQuadrupoleGD(CSubMatrix*      matrix_xx,
                 CSubMatrix*      matrix_xy,
                 CSubMatrix*      matrix_xz,
                 CSubMatrix*      matrix_yy,
                 CSubMatrix*      matrix_yz,
                 CSubMatrix*      matrix_zz,
                 const TPoint3D&  point,
                 const CGtoBlock& bra_gto_block,
                 const CGtoBlock& ket_gto_block,
                 const bool       ang_order,
                 const int64_t    bra_first,
                 const int64_t    bra_last) -> void
{
    // spherical transformation factors

    const double f4_35 = 4.0 * std::sqrt(35);

    const double f4_17 = 4.0 * std::sqrt(17.5);

    const double f4_5 = 4.0 * std::sqrt(5.0);

    const double f4_2 = 4.0 * std::sqrt(2.5);

    const double f2_3 = 2.0 * std::sqrt(3.0);

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

        simd::loadCoordinates(ket_coords_x, ket_coords_y, ket_coords_z, ket_gto_coords, ket_first, ket_last);

        for (int64_t j = bra_first; j < bra_last; j++)
        {
            const auto bra_coord = bra_gto_coords[j];

            // compute primitive integrals block (XXXX_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXX_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXX_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXX_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXX_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXX_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXX_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXY_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXY_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXY_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXY_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXY_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXY_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXXZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYY_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYY_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYY_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYY_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYY_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYY_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXYZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XXZZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXZZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXZZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XXZZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXZZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XXZZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYY_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYY_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYY_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYY_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYY_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYY_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYYZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XYZZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XYZZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYZZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XYZZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYZZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XYZZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_XZZZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_XZZZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XZZZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_XZZZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XZZZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_XZZZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYY_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYY_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYY_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYY_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYY_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYY_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYYZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_YYZZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_YYZZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYZZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_YYZZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYZZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YYZZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_YZZZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_YZZZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YZZZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_YZZZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YZZZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_YZZZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XX)

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

                    mpol::compPrimitiveQuadrupoleGD_ZZZZ_XX(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XY)

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

                    mpol::compPrimitiveQuadrupoleGD_ZZZZ_XY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XZ)

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

                    mpol::compPrimitiveQuadrupoleGD_ZZZZ_XZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_YY)

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

                    mpol::compPrimitiveQuadrupoleGD_ZZZZ_YY(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_YZ)

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

                    mpol::compPrimitiveQuadrupoleGD_ZZZZ_YZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_ZZ)

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

                    mpol::compPrimitiveQuadrupoleGD_ZZZZ_ZZ(buffer_xx,
                                                            buffer_xy,
                                                            buffer_xz,
                                                            buffer_yy,
                                                            buffer_yz,
                                                            buffer_zz,
                                                            point,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace mpol
