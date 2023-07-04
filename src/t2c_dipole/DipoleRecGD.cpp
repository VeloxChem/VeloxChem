#include "DipoleRecGD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveDipoleGD_XXXX_XX.hpp"
#include "PrimitiveDipoleGD_XXXX_XY.hpp"
#include "PrimitiveDipoleGD_XXXX_XZ.hpp"
#include "PrimitiveDipoleGD_XXXX_YY.hpp"
#include "PrimitiveDipoleGD_XXXX_YZ.hpp"
#include "PrimitiveDipoleGD_XXXX_ZZ.hpp"
#include "PrimitiveDipoleGD_XXXY_XX.hpp"
#include "PrimitiveDipoleGD_XXXY_XY.hpp"
#include "PrimitiveDipoleGD_XXXY_XZ.hpp"
#include "PrimitiveDipoleGD_XXXY_YY.hpp"
#include "PrimitiveDipoleGD_XXXY_YZ.hpp"
#include "PrimitiveDipoleGD_XXXY_ZZ.hpp"
#include "PrimitiveDipoleGD_XXXZ_XX.hpp"
#include "PrimitiveDipoleGD_XXXZ_XY.hpp"
#include "PrimitiveDipoleGD_XXXZ_XZ.hpp"
#include "PrimitiveDipoleGD_XXXZ_YY.hpp"
#include "PrimitiveDipoleGD_XXXZ_YZ.hpp"
#include "PrimitiveDipoleGD_XXXZ_ZZ.hpp"
#include "PrimitiveDipoleGD_XXYY_XX.hpp"
#include "PrimitiveDipoleGD_XXYY_XY.hpp"
#include "PrimitiveDipoleGD_XXYY_XZ.hpp"
#include "PrimitiveDipoleGD_XXYY_YY.hpp"
#include "PrimitiveDipoleGD_XXYY_YZ.hpp"
#include "PrimitiveDipoleGD_XXYY_ZZ.hpp"
#include "PrimitiveDipoleGD_XXYZ_XX.hpp"
#include "PrimitiveDipoleGD_XXYZ_XY.hpp"
#include "PrimitiveDipoleGD_XXYZ_XZ.hpp"
#include "PrimitiveDipoleGD_XXYZ_YY.hpp"
#include "PrimitiveDipoleGD_XXYZ_YZ.hpp"
#include "PrimitiveDipoleGD_XXYZ_ZZ.hpp"
#include "PrimitiveDipoleGD_XXZZ_XX.hpp"
#include "PrimitiveDipoleGD_XXZZ_XY.hpp"
#include "PrimitiveDipoleGD_XXZZ_XZ.hpp"
#include "PrimitiveDipoleGD_XXZZ_YY.hpp"
#include "PrimitiveDipoleGD_XXZZ_YZ.hpp"
#include "PrimitiveDipoleGD_XXZZ_ZZ.hpp"
#include "PrimitiveDipoleGD_XYYY_XX.hpp"
#include "PrimitiveDipoleGD_XYYY_XY.hpp"
#include "PrimitiveDipoleGD_XYYY_XZ.hpp"
#include "PrimitiveDipoleGD_XYYY_YY.hpp"
#include "PrimitiveDipoleGD_XYYY_YZ.hpp"
#include "PrimitiveDipoleGD_XYYY_ZZ.hpp"
#include "PrimitiveDipoleGD_XYYZ_XX.hpp"
#include "PrimitiveDipoleGD_XYYZ_XY.hpp"
#include "PrimitiveDipoleGD_XYYZ_XZ.hpp"
#include "PrimitiveDipoleGD_XYYZ_YY.hpp"
#include "PrimitiveDipoleGD_XYYZ_YZ.hpp"
#include "PrimitiveDipoleGD_XYYZ_ZZ.hpp"
#include "PrimitiveDipoleGD_XYZZ_XX.hpp"
#include "PrimitiveDipoleGD_XYZZ_XY.hpp"
#include "PrimitiveDipoleGD_XYZZ_XZ.hpp"
#include "PrimitiveDipoleGD_XYZZ_YY.hpp"
#include "PrimitiveDipoleGD_XYZZ_YZ.hpp"
#include "PrimitiveDipoleGD_XYZZ_ZZ.hpp"
#include "PrimitiveDipoleGD_XZZZ_XX.hpp"
#include "PrimitiveDipoleGD_XZZZ_XY.hpp"
#include "PrimitiveDipoleGD_XZZZ_XZ.hpp"
#include "PrimitiveDipoleGD_XZZZ_YY.hpp"
#include "PrimitiveDipoleGD_XZZZ_YZ.hpp"
#include "PrimitiveDipoleGD_XZZZ_ZZ.hpp"
#include "PrimitiveDipoleGD_YYYY_XX.hpp"
#include "PrimitiveDipoleGD_YYYY_XY.hpp"
#include "PrimitiveDipoleGD_YYYY_XZ.hpp"
#include "PrimitiveDipoleGD_YYYY_YY.hpp"
#include "PrimitiveDipoleGD_YYYY_YZ.hpp"
#include "PrimitiveDipoleGD_YYYY_ZZ.hpp"
#include "PrimitiveDipoleGD_YYYZ_XX.hpp"
#include "PrimitiveDipoleGD_YYYZ_XY.hpp"
#include "PrimitiveDipoleGD_YYYZ_XZ.hpp"
#include "PrimitiveDipoleGD_YYYZ_YY.hpp"
#include "PrimitiveDipoleGD_YYYZ_YZ.hpp"
#include "PrimitiveDipoleGD_YYYZ_ZZ.hpp"
#include "PrimitiveDipoleGD_YYZZ_XX.hpp"
#include "PrimitiveDipoleGD_YYZZ_XY.hpp"
#include "PrimitiveDipoleGD_YYZZ_XZ.hpp"
#include "PrimitiveDipoleGD_YYZZ_YY.hpp"
#include "PrimitiveDipoleGD_YYZZ_YZ.hpp"
#include "PrimitiveDipoleGD_YYZZ_ZZ.hpp"
#include "PrimitiveDipoleGD_YZZZ_XX.hpp"
#include "PrimitiveDipoleGD_YZZZ_XY.hpp"
#include "PrimitiveDipoleGD_YZZZ_XZ.hpp"
#include "PrimitiveDipoleGD_YZZZ_YY.hpp"
#include "PrimitiveDipoleGD_YZZZ_YZ.hpp"
#include "PrimitiveDipoleGD_YZZZ_ZZ.hpp"
#include "PrimitiveDipoleGD_ZZZZ_XX.hpp"
#include "PrimitiveDipoleGD_ZZZZ_XY.hpp"
#include "PrimitiveDipoleGD_ZZZZ_XZ.hpp"
#include "PrimitiveDipoleGD_ZZZZ_YY.hpp"
#include "PrimitiveDipoleGD_ZZZZ_YZ.hpp"
#include "PrimitiveDipoleGD_ZZZZ_ZZ.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compDipoleGD(CSubMatrix*      matrix_x,
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

            // compute primitive integrals block (XXXX_XX)

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

                    mpol::compPrimitiveDipoleGD_XXXX_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XY)

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

                    mpol::compPrimitiveDipoleGD_XXXX_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XZ)

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

                    mpol::compPrimitiveDipoleGD_XXXX_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_YY)

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

                    mpol::compPrimitiveDipoleGD_XXXX_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_YZ)

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

                    mpol::compPrimitiveDipoleGD_XXXX_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XXXX_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XX)

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

                    mpol::compPrimitiveDipoleGD_XXXY_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XY)

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

                    mpol::compPrimitiveDipoleGD_XXXY_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XZ)

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

                    mpol::compPrimitiveDipoleGD_XXXY_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_YY)

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

                    mpol::compPrimitiveDipoleGD_XXXY_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_YZ)

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

                    mpol::compPrimitiveDipoleGD_XXXY_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XXXY_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XX)

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

                    mpol::compPrimitiveDipoleGD_XXXZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XY)

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

                    mpol::compPrimitiveDipoleGD_XXXZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_XXXZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_YY)

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

                    mpol::compPrimitiveDipoleGD_XXXZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_XXXZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XXXZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XX)

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

                    mpol::compPrimitiveDipoleGD_XXYY_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XY)

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

                    mpol::compPrimitiveDipoleGD_XXYY_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XZ)

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

                    mpol::compPrimitiveDipoleGD_XXYY_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_YY)

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

                    mpol::compPrimitiveDipoleGD_XXYY_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_YZ)

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

                    mpol::compPrimitiveDipoleGD_XXYY_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XXYY_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XX)

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

                    mpol::compPrimitiveDipoleGD_XXYZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XY)

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

                    mpol::compPrimitiveDipoleGD_XXYZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_XXYZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_YY)

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

                    mpol::compPrimitiveDipoleGD_XXYZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_XXYZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XXYZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XX)

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

                    mpol::compPrimitiveDipoleGD_XXZZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XY)

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

                    mpol::compPrimitiveDipoleGD_XXZZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_XXZZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_YY)

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

                    mpol::compPrimitiveDipoleGD_XXZZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_XXZZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XXZZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XX)

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

                    mpol::compPrimitiveDipoleGD_XYYY_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XY)

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

                    mpol::compPrimitiveDipoleGD_XYYY_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XZ)

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

                    mpol::compPrimitiveDipoleGD_XYYY_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_YY)

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

                    mpol::compPrimitiveDipoleGD_XYYY_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_YZ)

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

                    mpol::compPrimitiveDipoleGD_XYYY_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XYYY_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XX)

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

                    mpol::compPrimitiveDipoleGD_XYYZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XY)

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

                    mpol::compPrimitiveDipoleGD_XYYZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_XYYZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_YY)

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

                    mpol::compPrimitiveDipoleGD_XYYZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_XYYZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XYYZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XX)

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

                    mpol::compPrimitiveDipoleGD_XYZZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XY)

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

                    mpol::compPrimitiveDipoleGD_XYZZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_XYZZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_YY)

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

                    mpol::compPrimitiveDipoleGD_XYZZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_XYZZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XYZZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XX)

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

                    mpol::compPrimitiveDipoleGD_XZZZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XY)

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

                    mpol::compPrimitiveDipoleGD_XZZZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_XZZZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_YY)

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

                    mpol::compPrimitiveDipoleGD_XZZZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_XZZZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_XZZZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XX)

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

                    mpol::compPrimitiveDipoleGD_YYYY_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XY)

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

                    mpol::compPrimitiveDipoleGD_YYYY_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XZ)

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

                    mpol::compPrimitiveDipoleGD_YYYY_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_YY)

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

                    mpol::compPrimitiveDipoleGD_YYYY_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_YZ)

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

                    mpol::compPrimitiveDipoleGD_YYYY_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_ZZ)

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

                    mpol::compPrimitiveDipoleGD_YYYY_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XX)

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

                    mpol::compPrimitiveDipoleGD_YYYZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XY)

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

                    mpol::compPrimitiveDipoleGD_YYYZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_YYYZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_YY)

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

                    mpol::compPrimitiveDipoleGD_YYYZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_YYYZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_YYYZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XX)

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

                    mpol::compPrimitiveDipoleGD_YYZZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XY)

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

                    mpol::compPrimitiveDipoleGD_YYZZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_YYZZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_YY)

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

                    mpol::compPrimitiveDipoleGD_YYZZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_YYZZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_YYZZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XX)

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

                    mpol::compPrimitiveDipoleGD_YZZZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XY)

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

                    mpol::compPrimitiveDipoleGD_YZZZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_YZZZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_YY)

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

                    mpol::compPrimitiveDipoleGD_YZZZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_YZZZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_YZZZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XX)

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

                    mpol::compPrimitiveDipoleGD_ZZZZ_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XY)

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

                    mpol::compPrimitiveDipoleGD_ZZZZ_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XZ)

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

                    mpol::compPrimitiveDipoleGD_ZZZZ_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_YY)

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

                    mpol::compPrimitiveDipoleGD_ZZZZ_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_YZ)

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

                    mpol::compPrimitiveDipoleGD_ZZZZ_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_ZZ)

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

                    mpol::compPrimitiveDipoleGD_ZZZZ_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace mpol
