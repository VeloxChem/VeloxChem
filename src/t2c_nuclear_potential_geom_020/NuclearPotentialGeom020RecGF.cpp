#include "NuclearPotentialGeom020RecGF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "T2CDistributor.hpp"

#include "PrimitiveNuclearPotentialGeom020GF_XXXX_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXX_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXY_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXXZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYY_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXYZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XXZZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYY_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYYZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XYZZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_XZZZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYY_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYYZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YYZZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_YZZZ_ZZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_XXX.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_XXY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_XXZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_XYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_XYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_XZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_YYY.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_YYZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_YZZ.hpp"
#include "PrimitiveNuclearPotentialGeom020GF_ZZZZ_ZZZ.hpp"

namespace geom_npotrec { // geom_npotrec namespace

auto
compNuclearPotentialGeom020GF(      CSubMatrix* matrix_xx,
                                    CSubMatrix* matrix_xy,
                                    CSubMatrix* matrix_xz,
                                    CSubMatrix* matrix_yy,
                                    CSubMatrix* matrix_yz,
                                    CSubMatrix* matrix_zz,
                              const T2Tensor& quadrupole,
                              const TPoint3D& point,
                              const CGtoBlock&  bra_gto_block,
                              const CGtoBlock&  ket_gto_block,
                              const bool        ang_order,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void
{
    // spherical transformation factors

    const double f4_35 = 4.0 * std::sqrt(35);

    const double f4_17 = 4.0 * std::sqrt(17.5);

    const double f4_5 = 4.0 * std::sqrt(5.0);

    const double f4_2 = 4.0 * std::sqrt(2.5);

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

            // compute primitive integrals block (XXXX_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXX_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXX_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXY_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXXZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYY_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXYZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XXZZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYY_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                7, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYYZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XYZZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_XZZZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                8, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYY_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYYZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YYZZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_YZZZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XXX)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_XXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XXY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_XXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XXZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_XXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_XYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_XYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_XZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_XZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_YYY)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_YYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_YYZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_YYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_YZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_YZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ_ZZZ)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom020GF_ZZZZ_ZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, ang_order);

        }
    }
}

} // geom_npotrec namespace

