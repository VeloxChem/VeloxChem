#include "QuadrupoleRecDG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveQuadrupoleDG_XX_XXXX.hpp"
#include "PrimitiveQuadrupoleDG_XX_XXXY.hpp"
#include "PrimitiveQuadrupoleDG_XX_XXXZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_XXYY.hpp"
#include "PrimitiveQuadrupoleDG_XX_XXYZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_XXZZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_XYYY.hpp"
#include "PrimitiveQuadrupoleDG_XX_XYYZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_XYZZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_XZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_YYYY.hpp"
#include "PrimitiveQuadrupoleDG_XX_YYYZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_YYZZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_YZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XX_ZZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_XXXX.hpp"
#include "PrimitiveQuadrupoleDG_XY_XXXY.hpp"
#include "PrimitiveQuadrupoleDG_XY_XXXZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_XXYY.hpp"
#include "PrimitiveQuadrupoleDG_XY_XXYZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_XXZZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_XYYY.hpp"
#include "PrimitiveQuadrupoleDG_XY_XYYZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_XYZZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_XZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_YYYY.hpp"
#include "PrimitiveQuadrupoleDG_XY_YYYZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_YYZZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_YZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XY_ZZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XXXX.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XXXY.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XXXZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XXYY.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XXYZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XXZZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XYYY.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XYYZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XYZZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_XZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_YYYY.hpp"
#include "PrimitiveQuadrupoleDG_XZ_YYYZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_YYZZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_YZZZ.hpp"
#include "PrimitiveQuadrupoleDG_XZ_ZZZZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_XXXX.hpp"
#include "PrimitiveQuadrupoleDG_YY_XXXY.hpp"
#include "PrimitiveQuadrupoleDG_YY_XXXZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_XXYY.hpp"
#include "PrimitiveQuadrupoleDG_YY_XXYZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_XXZZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_XYYY.hpp"
#include "PrimitiveQuadrupoleDG_YY_XYYZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_XYZZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_XZZZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_YYYY.hpp"
#include "PrimitiveQuadrupoleDG_YY_YYYZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_YYZZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_YZZZ.hpp"
#include "PrimitiveQuadrupoleDG_YY_ZZZZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XXXX.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XXXY.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XXXZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XXYY.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XXYZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XXZZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XYYY.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XYYZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XYZZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_XZZZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_YYYY.hpp"
#include "PrimitiveQuadrupoleDG_YZ_YYYZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_YYZZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_YZZZ.hpp"
#include "PrimitiveQuadrupoleDG_YZ_ZZZZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XXXX.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XXXY.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XXXZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XXYY.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XXYZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XXZZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XYYY.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XYYZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XYZZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_XZZZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_YYYY.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_YYYZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_YYZZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_YZZZ.hpp"
#include "PrimitiveQuadrupoleDG_ZZ_ZZZZ.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compQuadrupoleDG(CSubMatrix*      matrix_xx,
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

    const double f2_3 = 2.0 * std::sqrt(3.0);

    const double f4_35 = 4.0 * std::sqrt(35);

    const double f4_17 = 4.0 * std::sqrt(17.5);

    const double f4_5 = 4.0 * std::sqrt(5.0);

    const double f4_2 = 4.0 * std::sqrt(2.5);

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

            // compute primitive integrals block (XX_XXXX)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XXXY)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XXXZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XXYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XXYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XXZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XYZZ(buffer_xx,
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
                matrix_xx, buffer_xx, 0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_YYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_YYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_YYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_YZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_ZZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XX_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXXX)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXXY)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXXZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_YYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_YYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_YYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_YZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_ZZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XY_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXXX)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXXY)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXXZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_YYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_YYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_YYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_YZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_ZZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_XZ_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXXX)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXXY)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXXZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXYY)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XYZZ(buffer_xx,
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
                matrix_xx, buffer_xx, -0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_YYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_YYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_YYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_YZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_ZZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YY_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXXX)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXXY)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXXZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXYY)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -f2_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_YYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_YYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_YYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f2_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_YZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_ZZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_YZ_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f2_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXXX)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXXY)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXXZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXYY)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_YYYY)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_YYYZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_YYZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_YZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_ZZZZ)

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

                    mpol::compPrimitiveQuadrupoleDG_ZZ_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace mpol
