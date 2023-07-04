#include "DipoleRecFG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveDipoleFG_XXX_XXXX.hpp"
#include "PrimitiveDipoleFG_XXX_XXXY.hpp"
#include "PrimitiveDipoleFG_XXX_XXXZ.hpp"
#include "PrimitiveDipoleFG_XXX_XXYY.hpp"
#include "PrimitiveDipoleFG_XXX_XXYZ.hpp"
#include "PrimitiveDipoleFG_XXX_XXZZ.hpp"
#include "PrimitiveDipoleFG_XXX_XYYY.hpp"
#include "PrimitiveDipoleFG_XXX_XYYZ.hpp"
#include "PrimitiveDipoleFG_XXX_XYZZ.hpp"
#include "PrimitiveDipoleFG_XXX_XZZZ.hpp"
#include "PrimitiveDipoleFG_XXX_YYYY.hpp"
#include "PrimitiveDipoleFG_XXX_YYYZ.hpp"
#include "PrimitiveDipoleFG_XXX_YYZZ.hpp"
#include "PrimitiveDipoleFG_XXX_YZZZ.hpp"
#include "PrimitiveDipoleFG_XXX_ZZZZ.hpp"
#include "PrimitiveDipoleFG_XXY_XXXX.hpp"
#include "PrimitiveDipoleFG_XXY_XXXY.hpp"
#include "PrimitiveDipoleFG_XXY_XXXZ.hpp"
#include "PrimitiveDipoleFG_XXY_XXYY.hpp"
#include "PrimitiveDipoleFG_XXY_XXYZ.hpp"
#include "PrimitiveDipoleFG_XXY_XXZZ.hpp"
#include "PrimitiveDipoleFG_XXY_XYYY.hpp"
#include "PrimitiveDipoleFG_XXY_XYYZ.hpp"
#include "PrimitiveDipoleFG_XXY_XYZZ.hpp"
#include "PrimitiveDipoleFG_XXY_XZZZ.hpp"
#include "PrimitiveDipoleFG_XXY_YYYY.hpp"
#include "PrimitiveDipoleFG_XXY_YYYZ.hpp"
#include "PrimitiveDipoleFG_XXY_YYZZ.hpp"
#include "PrimitiveDipoleFG_XXY_YZZZ.hpp"
#include "PrimitiveDipoleFG_XXY_ZZZZ.hpp"
#include "PrimitiveDipoleFG_XXZ_XXXX.hpp"
#include "PrimitiveDipoleFG_XXZ_XXXY.hpp"
#include "PrimitiveDipoleFG_XXZ_XXXZ.hpp"
#include "PrimitiveDipoleFG_XXZ_XXYY.hpp"
#include "PrimitiveDipoleFG_XXZ_XXYZ.hpp"
#include "PrimitiveDipoleFG_XXZ_XXZZ.hpp"
#include "PrimitiveDipoleFG_XXZ_XYYY.hpp"
#include "PrimitiveDipoleFG_XXZ_XYYZ.hpp"
#include "PrimitiveDipoleFG_XXZ_XYZZ.hpp"
#include "PrimitiveDipoleFG_XXZ_XZZZ.hpp"
#include "PrimitiveDipoleFG_XXZ_YYYY.hpp"
#include "PrimitiveDipoleFG_XXZ_YYYZ.hpp"
#include "PrimitiveDipoleFG_XXZ_YYZZ.hpp"
#include "PrimitiveDipoleFG_XXZ_YZZZ.hpp"
#include "PrimitiveDipoleFG_XXZ_ZZZZ.hpp"
#include "PrimitiveDipoleFG_XYY_XXXX.hpp"
#include "PrimitiveDipoleFG_XYY_XXXY.hpp"
#include "PrimitiveDipoleFG_XYY_XXXZ.hpp"
#include "PrimitiveDipoleFG_XYY_XXYY.hpp"
#include "PrimitiveDipoleFG_XYY_XXYZ.hpp"
#include "PrimitiveDipoleFG_XYY_XXZZ.hpp"
#include "PrimitiveDipoleFG_XYY_XYYY.hpp"
#include "PrimitiveDipoleFG_XYY_XYYZ.hpp"
#include "PrimitiveDipoleFG_XYY_XYZZ.hpp"
#include "PrimitiveDipoleFG_XYY_XZZZ.hpp"
#include "PrimitiveDipoleFG_XYY_YYYY.hpp"
#include "PrimitiveDipoleFG_XYY_YYYZ.hpp"
#include "PrimitiveDipoleFG_XYY_YYZZ.hpp"
#include "PrimitiveDipoleFG_XYY_YZZZ.hpp"
#include "PrimitiveDipoleFG_XYY_ZZZZ.hpp"
#include "PrimitiveDipoleFG_XYZ_XXXX.hpp"
#include "PrimitiveDipoleFG_XYZ_XXXY.hpp"
#include "PrimitiveDipoleFG_XYZ_XXXZ.hpp"
#include "PrimitiveDipoleFG_XYZ_XXYY.hpp"
#include "PrimitiveDipoleFG_XYZ_XXYZ.hpp"
#include "PrimitiveDipoleFG_XYZ_XXZZ.hpp"
#include "PrimitiveDipoleFG_XYZ_XYYY.hpp"
#include "PrimitiveDipoleFG_XYZ_XYYZ.hpp"
#include "PrimitiveDipoleFG_XYZ_XYZZ.hpp"
#include "PrimitiveDipoleFG_XYZ_XZZZ.hpp"
#include "PrimitiveDipoleFG_XYZ_YYYY.hpp"
#include "PrimitiveDipoleFG_XYZ_YYYZ.hpp"
#include "PrimitiveDipoleFG_XYZ_YYZZ.hpp"
#include "PrimitiveDipoleFG_XYZ_YZZZ.hpp"
#include "PrimitiveDipoleFG_XYZ_ZZZZ.hpp"
#include "PrimitiveDipoleFG_XZZ_XXXX.hpp"
#include "PrimitiveDipoleFG_XZZ_XXXY.hpp"
#include "PrimitiveDipoleFG_XZZ_XXXZ.hpp"
#include "PrimitiveDipoleFG_XZZ_XXYY.hpp"
#include "PrimitiveDipoleFG_XZZ_XXYZ.hpp"
#include "PrimitiveDipoleFG_XZZ_XXZZ.hpp"
#include "PrimitiveDipoleFG_XZZ_XYYY.hpp"
#include "PrimitiveDipoleFG_XZZ_XYYZ.hpp"
#include "PrimitiveDipoleFG_XZZ_XYZZ.hpp"
#include "PrimitiveDipoleFG_XZZ_XZZZ.hpp"
#include "PrimitiveDipoleFG_XZZ_YYYY.hpp"
#include "PrimitiveDipoleFG_XZZ_YYYZ.hpp"
#include "PrimitiveDipoleFG_XZZ_YYZZ.hpp"
#include "PrimitiveDipoleFG_XZZ_YZZZ.hpp"
#include "PrimitiveDipoleFG_XZZ_ZZZZ.hpp"
#include "PrimitiveDipoleFG_YYY_XXXX.hpp"
#include "PrimitiveDipoleFG_YYY_XXXY.hpp"
#include "PrimitiveDipoleFG_YYY_XXXZ.hpp"
#include "PrimitiveDipoleFG_YYY_XXYY.hpp"
#include "PrimitiveDipoleFG_YYY_XXYZ.hpp"
#include "PrimitiveDipoleFG_YYY_XXZZ.hpp"
#include "PrimitiveDipoleFG_YYY_XYYY.hpp"
#include "PrimitiveDipoleFG_YYY_XYYZ.hpp"
#include "PrimitiveDipoleFG_YYY_XYZZ.hpp"
#include "PrimitiveDipoleFG_YYY_XZZZ.hpp"
#include "PrimitiveDipoleFG_YYY_YYYY.hpp"
#include "PrimitiveDipoleFG_YYY_YYYZ.hpp"
#include "PrimitiveDipoleFG_YYY_YYZZ.hpp"
#include "PrimitiveDipoleFG_YYY_YZZZ.hpp"
#include "PrimitiveDipoleFG_YYY_ZZZZ.hpp"
#include "PrimitiveDipoleFG_YYZ_XXXX.hpp"
#include "PrimitiveDipoleFG_YYZ_XXXY.hpp"
#include "PrimitiveDipoleFG_YYZ_XXXZ.hpp"
#include "PrimitiveDipoleFG_YYZ_XXYY.hpp"
#include "PrimitiveDipoleFG_YYZ_XXYZ.hpp"
#include "PrimitiveDipoleFG_YYZ_XXZZ.hpp"
#include "PrimitiveDipoleFG_YYZ_XYYY.hpp"
#include "PrimitiveDipoleFG_YYZ_XYYZ.hpp"
#include "PrimitiveDipoleFG_YYZ_XYZZ.hpp"
#include "PrimitiveDipoleFG_YYZ_XZZZ.hpp"
#include "PrimitiveDipoleFG_YYZ_YYYY.hpp"
#include "PrimitiveDipoleFG_YYZ_YYYZ.hpp"
#include "PrimitiveDipoleFG_YYZ_YYZZ.hpp"
#include "PrimitiveDipoleFG_YYZ_YZZZ.hpp"
#include "PrimitiveDipoleFG_YYZ_ZZZZ.hpp"
#include "PrimitiveDipoleFG_YZZ_XXXX.hpp"
#include "PrimitiveDipoleFG_YZZ_XXXY.hpp"
#include "PrimitiveDipoleFG_YZZ_XXXZ.hpp"
#include "PrimitiveDipoleFG_YZZ_XXYY.hpp"
#include "PrimitiveDipoleFG_YZZ_XXYZ.hpp"
#include "PrimitiveDipoleFG_YZZ_XXZZ.hpp"
#include "PrimitiveDipoleFG_YZZ_XYYY.hpp"
#include "PrimitiveDipoleFG_YZZ_XYYZ.hpp"
#include "PrimitiveDipoleFG_YZZ_XYZZ.hpp"
#include "PrimitiveDipoleFG_YZZ_XZZZ.hpp"
#include "PrimitiveDipoleFG_YZZ_YYYY.hpp"
#include "PrimitiveDipoleFG_YZZ_YYYZ.hpp"
#include "PrimitiveDipoleFG_YZZ_YYZZ.hpp"
#include "PrimitiveDipoleFG_YZZ_YZZZ.hpp"
#include "PrimitiveDipoleFG_YZZ_ZZZZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_XXXX.hpp"
#include "PrimitiveDipoleFG_ZZZ_XXXY.hpp"
#include "PrimitiveDipoleFG_ZZZ_XXXZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_XXYY.hpp"
#include "PrimitiveDipoleFG_ZZZ_XXYZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_XXZZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_XYYY.hpp"
#include "PrimitiveDipoleFG_ZZZ_XYYZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_XYZZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_XZZZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_YYYY.hpp"
#include "PrimitiveDipoleFG_ZZZ_YYYZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_YYZZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_YZZZ.hpp"
#include "PrimitiveDipoleFG_ZZZ_ZZZZ.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compDipoleFG(CSubMatrix*      matrix_x,
             CSubMatrix*      matrix_y,
             CSubMatrix*      matrix_z,
             const TPoint3D&  point,
             const CGtoBlock& bra_gto_block,
             const CGtoBlock& ket_gto_block,
             const bool       ang_order,
             const int64_t    bra_first,
             const int64_t    bra_last) -> void
{
    // spherical transformation factors

    const double f3_5 = std::sqrt(2.5);

    const double f3_15 = 2.0 * std::sqrt(15.0);

    const double f3_3 = std::sqrt(1.5);

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

    alignas(64) TDoubleArray buffer_x;

    alignas(64) TDoubleArray buffer_y;

    alignas(64) TDoubleArray buffer_z;

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

            // compute primitive integrals block (XXX_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXX_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXY_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XXZ_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYY_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XYZ_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_XZZ_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYY_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YYZ_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_YZZ_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XXXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XXXX(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XXXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XXXY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XXXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XXXZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XXYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XXYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XXYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XXYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XXZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XXZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_XZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_YYYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_YYYY(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_YYYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_YYYZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_YYZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_YYZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_YZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_YZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_ZZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFG_ZZZ_ZZZZ(buffer_x,
                                                         buffer_y,
                                                         buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace mpol
