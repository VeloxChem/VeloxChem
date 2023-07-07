#include "OverlapGeom400RecPF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom400PF_X_XXX.hpp"
#include "PrimitiveOverlapGeom400PF_X_XXY.hpp"
#include "PrimitiveOverlapGeom400PF_X_XXZ.hpp"
#include "PrimitiveOverlapGeom400PF_X_XYY.hpp"
#include "PrimitiveOverlapGeom400PF_X_XYZ.hpp"
#include "PrimitiveOverlapGeom400PF_X_XZZ.hpp"
#include "PrimitiveOverlapGeom400PF_X_YYY.hpp"
#include "PrimitiveOverlapGeom400PF_X_YYZ.hpp"
#include "PrimitiveOverlapGeom400PF_X_YZZ.hpp"
#include "PrimitiveOverlapGeom400PF_X_ZZZ.hpp"
#include "PrimitiveOverlapGeom400PF_Y_XXX.hpp"
#include "PrimitiveOverlapGeom400PF_Y_XXY.hpp"
#include "PrimitiveOverlapGeom400PF_Y_XXZ.hpp"
#include "PrimitiveOverlapGeom400PF_Y_XYY.hpp"
#include "PrimitiveOverlapGeom400PF_Y_XYZ.hpp"
#include "PrimitiveOverlapGeom400PF_Y_XZZ.hpp"
#include "PrimitiveOverlapGeom400PF_Y_YYY.hpp"
#include "PrimitiveOverlapGeom400PF_Y_YYZ.hpp"
#include "PrimitiveOverlapGeom400PF_Y_YZZ.hpp"
#include "PrimitiveOverlapGeom400PF_Y_ZZZ.hpp"
#include "PrimitiveOverlapGeom400PF_Z_XXX.hpp"
#include "PrimitiveOverlapGeom400PF_Z_XXY.hpp"
#include "PrimitiveOverlapGeom400PF_Z_XXZ.hpp"
#include "PrimitiveOverlapGeom400PF_Z_XYY.hpp"
#include "PrimitiveOverlapGeom400PF_Z_XYZ.hpp"
#include "PrimitiveOverlapGeom400PF_Z_XZZ.hpp"
#include "PrimitiveOverlapGeom400PF_Z_YYY.hpp"
#include "PrimitiveOverlapGeom400PF_Z_YYZ.hpp"
#include "PrimitiveOverlapGeom400PF_Z_YZZ.hpp"
#include "PrimitiveOverlapGeom400PF_Z_ZZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom400PF(CSubMatrix*      matrix_xxxx,
                     CSubMatrix*      matrix_xxxy,
                     CSubMatrix*      matrix_xxxz,
                     CSubMatrix*      matrix_xxyy,
                     CSubMatrix*      matrix_xxyz,
                     CSubMatrix*      matrix_xxzz,
                     CSubMatrix*      matrix_xyyy,
                     CSubMatrix*      matrix_xyyz,
                     CSubMatrix*      matrix_xyzz,
                     CSubMatrix*      matrix_xzzz,
                     CSubMatrix*      matrix_yyyy,
                     CSubMatrix*      matrix_yyyz,
                     CSubMatrix*      matrix_yyzz,
                     CSubMatrix*      matrix_yzzz,
                     CSubMatrix*      matrix_zzzz,
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

    alignas(64) TDoubleArray buffer_xxxx;

    alignas(64) TDoubleArray buffer_xxxy;

    alignas(64) TDoubleArray buffer_xxxz;

    alignas(64) TDoubleArray buffer_xxyy;

    alignas(64) TDoubleArray buffer_xxyz;

    alignas(64) TDoubleArray buffer_xxzz;

    alignas(64) TDoubleArray buffer_xyyy;

    alignas(64) TDoubleArray buffer_xyyz;

    alignas(64) TDoubleArray buffer_xyzz;

    alignas(64) TDoubleArray buffer_xzzz;

    alignas(64) TDoubleArray buffer_yyyy;

    alignas(64) TDoubleArray buffer_yyyz;

    alignas(64) TDoubleArray buffer_yyzz;

    alignas(64) TDoubleArray buffer_yzzz;

    alignas(64) TDoubleArray buffer_zzzz;

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

            // compute primitive integrals block (X_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_XXX(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_XXY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_XXZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_XYY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_XYZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_XZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_YYY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_YYZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_YZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_X_ZZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_XXX(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_XXY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_XXZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_XYY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_XYZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_XZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_YYY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_YYZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_YZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Y_ZZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_XXX(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_XXY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_XXZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_XYY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_XYZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_XZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_YYY(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_YYZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_YZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400PF_Z_ZZZ(buffer_xxxx,
                                                                buffer_xxxy,
                                                                buffer_xxxz,
                                                                buffer_xxyy,
                                                                buffer_xxyz,
                                                                buffer_xxzz,
                                                                buffer_xyyy,
                                                                buffer_xyyz,
                                                                buffer_xyzz,
                                                                buffer_xzzz,
                                                                buffer_yyyy,
                                                                buffer_yyyz,
                                                                buffer_yyzz,
                                                                buffer_yzzz,
                                                                buffer_zzzz,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace ovlrec
