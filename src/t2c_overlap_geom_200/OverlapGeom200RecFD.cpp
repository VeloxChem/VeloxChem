#include "OverlapGeom200RecFD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom200FD_XXX_XX.hpp"
#include "PrimitiveOverlapGeom200FD_XXX_XY.hpp"
#include "PrimitiveOverlapGeom200FD_XXX_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXX_YY.hpp"
#include "PrimitiveOverlapGeom200FD_XXX_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXX_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXY_XX.hpp"
#include "PrimitiveOverlapGeom200FD_XXY_XY.hpp"
#include "PrimitiveOverlapGeom200FD_XXY_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXY_YY.hpp"
#include "PrimitiveOverlapGeom200FD_XXY_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXY_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXZ_XX.hpp"
#include "PrimitiveOverlapGeom200FD_XXZ_XY.hpp"
#include "PrimitiveOverlapGeom200FD_XXZ_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXZ_YY.hpp"
#include "PrimitiveOverlapGeom200FD_XXZ_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_XXZ_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_XYY_XX.hpp"
#include "PrimitiveOverlapGeom200FD_XYY_XY.hpp"
#include "PrimitiveOverlapGeom200FD_XYY_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_XYY_YY.hpp"
#include "PrimitiveOverlapGeom200FD_XYY_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_XYY_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_XYZ_XX.hpp"
#include "PrimitiveOverlapGeom200FD_XYZ_XY.hpp"
#include "PrimitiveOverlapGeom200FD_XYZ_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_XYZ_YY.hpp"
#include "PrimitiveOverlapGeom200FD_XYZ_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_XYZ_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_XZZ_XX.hpp"
#include "PrimitiveOverlapGeom200FD_XZZ_XY.hpp"
#include "PrimitiveOverlapGeom200FD_XZZ_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_XZZ_YY.hpp"
#include "PrimitiveOverlapGeom200FD_XZZ_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_XZZ_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_YYY_XX.hpp"
#include "PrimitiveOverlapGeom200FD_YYY_XY.hpp"
#include "PrimitiveOverlapGeom200FD_YYY_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_YYY_YY.hpp"
#include "PrimitiveOverlapGeom200FD_YYY_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_YYY_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_YYZ_XX.hpp"
#include "PrimitiveOverlapGeom200FD_YYZ_XY.hpp"
#include "PrimitiveOverlapGeom200FD_YYZ_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_YYZ_YY.hpp"
#include "PrimitiveOverlapGeom200FD_YYZ_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_YYZ_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_YZZ_XX.hpp"
#include "PrimitiveOverlapGeom200FD_YZZ_XY.hpp"
#include "PrimitiveOverlapGeom200FD_YZZ_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_YZZ_YY.hpp"
#include "PrimitiveOverlapGeom200FD_YZZ_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_YZZ_ZZ.hpp"
#include "PrimitiveOverlapGeom200FD_ZZZ_XX.hpp"
#include "PrimitiveOverlapGeom200FD_ZZZ_XY.hpp"
#include "PrimitiveOverlapGeom200FD_ZZZ_XZ.hpp"
#include "PrimitiveOverlapGeom200FD_ZZZ_YY.hpp"
#include "PrimitiveOverlapGeom200FD_ZZZ_YZ.hpp"
#include "PrimitiveOverlapGeom200FD_ZZZ_ZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom200FD(CSubMatrix*      matrix_xx,
                     CSubMatrix*      matrix_xy,
                     CSubMatrix*      matrix_xz,
                     CSubMatrix*      matrix_yy,
                     CSubMatrix*      matrix_yz,
                     CSubMatrix*      matrix_zz,
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

            // compute primitive integrals block (XXX_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXX_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXX_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXX_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXX_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXX_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXX_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXY_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXY_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXY_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXY_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXY_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXY_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXZ_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXZ_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXZ_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXZ_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXZ_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XXZ_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYY_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYY_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYY_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYY_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 3.0 * f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYY_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYY_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYZ_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYZ_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYZ_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYZ_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYZ_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XYZ_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XZZ_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XZZ_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XZZ_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XZZ_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XZZ_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_XZZ_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYY_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYY_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYY_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYY_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYY_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYY_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYZ_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYZ_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYZ_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYZ_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 0.5 * f3_15 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYZ_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f3_15 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YYZ_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YZZ_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, 4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YZZ_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YZZ_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YZZ_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xx, buffer_xx, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xy, buffer_xy, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xz, buffer_xz, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yy, buffer_yy, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yz, buffer_yz, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zz, buffer_zz, -4.0 * f3_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YZZ_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_YZZ_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_ZZZ_XX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_ZZZ_XY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_ZZZ_XZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_ZZZ_YY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_ZZZ_YZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom200FD_ZZZ_ZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace ovlrec
