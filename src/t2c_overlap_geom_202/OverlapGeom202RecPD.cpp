#include "OverlapGeom202RecPD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom202PD_X_XX.hpp"
#include "PrimitiveOverlapGeom202PD_X_XY.hpp"
#include "PrimitiveOverlapGeom202PD_X_XZ.hpp"
#include "PrimitiveOverlapGeom202PD_X_YY.hpp"
#include "PrimitiveOverlapGeom202PD_X_YZ.hpp"
#include "PrimitiveOverlapGeom202PD_X_ZZ.hpp"
#include "PrimitiveOverlapGeom202PD_Y_XX.hpp"
#include "PrimitiveOverlapGeom202PD_Y_XY.hpp"
#include "PrimitiveOverlapGeom202PD_Y_XZ.hpp"
#include "PrimitiveOverlapGeom202PD_Y_YY.hpp"
#include "PrimitiveOverlapGeom202PD_Y_YZ.hpp"
#include "PrimitiveOverlapGeom202PD_Y_ZZ.hpp"
#include "PrimitiveOverlapGeom202PD_Z_XX.hpp"
#include "PrimitiveOverlapGeom202PD_Z_XY.hpp"
#include "PrimitiveOverlapGeom202PD_Z_XZ.hpp"
#include "PrimitiveOverlapGeom202PD_Z_YY.hpp"
#include "PrimitiveOverlapGeom202PD_Z_YZ.hpp"
#include "PrimitiveOverlapGeom202PD_Z_ZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom202PD(CSubMatrix*      matrix_xx_xx,
                     CSubMatrix*      matrix_xx_xy,
                     CSubMatrix*      matrix_xx_xz,
                     CSubMatrix*      matrix_xx_yy,
                     CSubMatrix*      matrix_xx_yz,
                     CSubMatrix*      matrix_xx_zz,
                     CSubMatrix*      matrix_xy_xx,
                     CSubMatrix*      matrix_xy_xy,
                     CSubMatrix*      matrix_xy_xz,
                     CSubMatrix*      matrix_xy_yy,
                     CSubMatrix*      matrix_xy_yz,
                     CSubMatrix*      matrix_xy_zz,
                     CSubMatrix*      matrix_xz_xx,
                     CSubMatrix*      matrix_xz_xy,
                     CSubMatrix*      matrix_xz_xz,
                     CSubMatrix*      matrix_xz_yy,
                     CSubMatrix*      matrix_xz_yz,
                     CSubMatrix*      matrix_xz_zz,
                     CSubMatrix*      matrix_yy_xx,
                     CSubMatrix*      matrix_yy_xy,
                     CSubMatrix*      matrix_yy_xz,
                     CSubMatrix*      matrix_yy_yy,
                     CSubMatrix*      matrix_yy_yz,
                     CSubMatrix*      matrix_yy_zz,
                     CSubMatrix*      matrix_yz_xx,
                     CSubMatrix*      matrix_yz_xy,
                     CSubMatrix*      matrix_yz_xz,
                     CSubMatrix*      matrix_yz_yy,
                     CSubMatrix*      matrix_yz_yz,
                     CSubMatrix*      matrix_yz_zz,
                     CSubMatrix*      matrix_zz_xx,
                     CSubMatrix*      matrix_zz_xy,
                     CSubMatrix*      matrix_zz_xz,
                     CSubMatrix*      matrix_zz_yy,
                     CSubMatrix*      matrix_zz_yz,
                     CSubMatrix*      matrix_zz_zz,
                     const CGtoBlock& bra_gto_block,
                     const CGtoBlock& ket_gto_block,
                     const bool       ang_order,
                     const int64_t    bra_first,
                     const int64_t    bra_last) -> void
{
    // spherical transformation factors

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

    alignas(64) TDoubleArray buffer_xx_xx;

    alignas(64) TDoubleArray buffer_xx_xy;

    alignas(64) TDoubleArray buffer_xx_xz;

    alignas(64) TDoubleArray buffer_xx_yy;

    alignas(64) TDoubleArray buffer_xx_yz;

    alignas(64) TDoubleArray buffer_xx_zz;

    alignas(64) TDoubleArray buffer_xy_xx;

    alignas(64) TDoubleArray buffer_xy_xy;

    alignas(64) TDoubleArray buffer_xy_xz;

    alignas(64) TDoubleArray buffer_xy_yy;

    alignas(64) TDoubleArray buffer_xy_yz;

    alignas(64) TDoubleArray buffer_xy_zz;

    alignas(64) TDoubleArray buffer_xz_xx;

    alignas(64) TDoubleArray buffer_xz_xy;

    alignas(64) TDoubleArray buffer_xz_xz;

    alignas(64) TDoubleArray buffer_xz_yy;

    alignas(64) TDoubleArray buffer_xz_yz;

    alignas(64) TDoubleArray buffer_xz_zz;

    alignas(64) TDoubleArray buffer_yy_xx;

    alignas(64) TDoubleArray buffer_yy_xy;

    alignas(64) TDoubleArray buffer_yy_xz;

    alignas(64) TDoubleArray buffer_yy_yy;

    alignas(64) TDoubleArray buffer_yy_yz;

    alignas(64) TDoubleArray buffer_yy_zz;

    alignas(64) TDoubleArray buffer_yz_xx;

    alignas(64) TDoubleArray buffer_yz_xy;

    alignas(64) TDoubleArray buffer_yz_xz;

    alignas(64) TDoubleArray buffer_yz_yy;

    alignas(64) TDoubleArray buffer_yz_yz;

    alignas(64) TDoubleArray buffer_yz_zz;

    alignas(64) TDoubleArray buffer_zz_xx;

    alignas(64) TDoubleArray buffer_zz_xy;

    alignas(64) TDoubleArray buffer_zz_xz;

    alignas(64) TDoubleArray buffer_zz_yy;

    alignas(64) TDoubleArray buffer_zz_yz;

    alignas(64) TDoubleArray buffer_zz_zz;

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

            // compute primitive integrals block (X_XX)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_X_XX(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XY)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_X_XY(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_X_XZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YY)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_X_YY(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_X_YZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_ZZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_X_ZZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XX)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Y_XX(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XY)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Y_XY(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Y_XZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YY)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Y_YY(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Y_YZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_ZZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Y_ZZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XX)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Z_XX(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XY)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Z_XY(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Z_XZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YY)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Z_YY(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -1.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Z_YZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_ZZ)

            simd::zero(buffer_xx_xx);

            simd::zero(buffer_xx_xy);

            simd::zero(buffer_xx_xz);

            simd::zero(buffer_xx_yy);

            simd::zero(buffer_xx_yz);

            simd::zero(buffer_xx_zz);

            simd::zero(buffer_xy_xx);

            simd::zero(buffer_xy_xy);

            simd::zero(buffer_xy_xz);

            simd::zero(buffer_xy_yy);

            simd::zero(buffer_xy_yz);

            simd::zero(buffer_xy_zz);

            simd::zero(buffer_xz_xx);

            simd::zero(buffer_xz_xy);

            simd::zero(buffer_xz_xz);

            simd::zero(buffer_xz_yy);

            simd::zero(buffer_xz_yz);

            simd::zero(buffer_xz_zz);

            simd::zero(buffer_yy_xx);

            simd::zero(buffer_yy_xy);

            simd::zero(buffer_yy_xz);

            simd::zero(buffer_yy_yy);

            simd::zero(buffer_yy_yz);

            simd::zero(buffer_yy_zz);

            simd::zero(buffer_yz_xx);

            simd::zero(buffer_yz_xy);

            simd::zero(buffer_yz_xz);

            simd::zero(buffer_yz_yy);

            simd::zero(buffer_yz_yz);

            simd::zero(buffer_yz_zz);

            simd::zero(buffer_zz_xx);

            simd::zero(buffer_zz_xy);

            simd::zero(buffer_zz_xz);

            simd::zero(buffer_zz_yy);

            simd::zero(buffer_zz_yz);

            simd::zero(buffer_zz_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom202PD_Z_ZZ(buffer_xx_xx,
                                                               buffer_xx_xy,
                                                               buffer_xx_xz,
                                                               buffer_xx_yy,
                                                               buffer_xx_yz,
                                                               buffer_xx_zz,
                                                               buffer_xy_xx,
                                                               buffer_xy_xy,
                                                               buffer_xy_xz,
                                                               buffer_xy_yy,
                                                               buffer_xy_yz,
                                                               buffer_xy_zz,
                                                               buffer_xz_xx,
                                                               buffer_xz_xy,
                                                               buffer_xz_xz,
                                                               buffer_xz_yy,
                                                               buffer_xz_yz,
                                                               buffer_xz_zz,
                                                               buffer_yy_xx,
                                                               buffer_yy_xy,
                                                               buffer_yy_xz,
                                                               buffer_yy_yy,
                                                               buffer_yy_yz,
                                                               buffer_yy_zz,
                                                               buffer_yz_xx,
                                                               buffer_yz_xy,
                                                               buffer_yz_xz,
                                                               buffer_yz_yy,
                                                               buffer_yz_yz,
                                                               buffer_yz_zz,
                                                               buffer_zz_xx,
                                                               buffer_zz_xy,
                                                               buffer_zz_xz,
                                                               buffer_zz_yy,
                                                               buffer_zz_yz,
                                                               buffer_zz_zz,
                                                               bra_exp,
                                                               bra_norm,
                                                               bra_coord,
                                                               ket_exps,
                                                               ket_norms,
                                                               ket_coords_x,
                                                               ket_coords_y,
                                                               ket_coords_z,
                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xx_xx, buffer_xx_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xy, buffer_xx_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_xz, buffer_xx_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yy, buffer_xx_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_yz, buffer_xx_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_zz, buffer_xx_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xx, buffer_xy_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xy, buffer_xy_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_xz, buffer_xy_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yy, buffer_xy_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_yz, buffer_xy_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_zz, buffer_xy_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xx, buffer_xz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xy, buffer_xz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_xz, buffer_xz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yy, buffer_xz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_yz, buffer_xz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_zz, buffer_xz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xx, buffer_yy_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xy, buffer_yy_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_xz, buffer_yy_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yy, buffer_yy_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_yz, buffer_yy_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_zz, buffer_yy_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xx, buffer_yz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xy, buffer_yz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_xz, buffer_yz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yy, buffer_yz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_yz, buffer_yz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_zz, buffer_yz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xx, buffer_zz_xx, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xy, buffer_zz_xy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_xz, buffer_zz_xz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yy, buffer_zz_yy, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_yz, buffer_zz_yz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_zz, buffer_zz_zz, 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace ovlrec
