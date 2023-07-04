#include "QuadrupoleRecPG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveQuadrupolePG_X_XXXX.hpp"
#include "PrimitiveQuadrupolePG_X_XXXY.hpp"
#include "PrimitiveQuadrupolePG_X_XXXZ.hpp"
#include "PrimitiveQuadrupolePG_X_XXYY.hpp"
#include "PrimitiveQuadrupolePG_X_XXYZ.hpp"
#include "PrimitiveQuadrupolePG_X_XXZZ.hpp"
#include "PrimitiveQuadrupolePG_X_XYYY.hpp"
#include "PrimitiveQuadrupolePG_X_XYYZ.hpp"
#include "PrimitiveQuadrupolePG_X_XYZZ.hpp"
#include "PrimitiveQuadrupolePG_X_XZZZ.hpp"
#include "PrimitiveQuadrupolePG_X_YYYY.hpp"
#include "PrimitiveQuadrupolePG_X_YYYZ.hpp"
#include "PrimitiveQuadrupolePG_X_YYZZ.hpp"
#include "PrimitiveQuadrupolePG_X_YZZZ.hpp"
#include "PrimitiveQuadrupolePG_X_ZZZZ.hpp"
#include "PrimitiveQuadrupolePG_Y_XXXX.hpp"
#include "PrimitiveQuadrupolePG_Y_XXXY.hpp"
#include "PrimitiveQuadrupolePG_Y_XXXZ.hpp"
#include "PrimitiveQuadrupolePG_Y_XXYY.hpp"
#include "PrimitiveQuadrupolePG_Y_XXYZ.hpp"
#include "PrimitiveQuadrupolePG_Y_XXZZ.hpp"
#include "PrimitiveQuadrupolePG_Y_XYYY.hpp"
#include "PrimitiveQuadrupolePG_Y_XYYZ.hpp"
#include "PrimitiveQuadrupolePG_Y_XYZZ.hpp"
#include "PrimitiveQuadrupolePG_Y_XZZZ.hpp"
#include "PrimitiveQuadrupolePG_Y_YYYY.hpp"
#include "PrimitiveQuadrupolePG_Y_YYYZ.hpp"
#include "PrimitiveQuadrupolePG_Y_YYZZ.hpp"
#include "PrimitiveQuadrupolePG_Y_YZZZ.hpp"
#include "PrimitiveQuadrupolePG_Y_ZZZZ.hpp"
#include "PrimitiveQuadrupolePG_Z_XXXX.hpp"
#include "PrimitiveQuadrupolePG_Z_XXXY.hpp"
#include "PrimitiveQuadrupolePG_Z_XXXZ.hpp"
#include "PrimitiveQuadrupolePG_Z_XXYY.hpp"
#include "PrimitiveQuadrupolePG_Z_XXYZ.hpp"
#include "PrimitiveQuadrupolePG_Z_XXZZ.hpp"
#include "PrimitiveQuadrupolePG_Z_XYYY.hpp"
#include "PrimitiveQuadrupolePG_Z_XYYZ.hpp"
#include "PrimitiveQuadrupolePG_Z_XYZZ.hpp"
#include "PrimitiveQuadrupolePG_Z_XZZZ.hpp"
#include "PrimitiveQuadrupolePG_Z_YYYY.hpp"
#include "PrimitiveQuadrupolePG_Z_YYYZ.hpp"
#include "PrimitiveQuadrupolePG_Z_YYZZ.hpp"
#include "PrimitiveQuadrupolePG_Z_YZZZ.hpp"
#include "PrimitiveQuadrupolePG_Z_ZZZZ.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compQuadrupolePG(CSubMatrix*      matrix_xx,
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

            // compute primitive integrals block (X_XXXX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXXY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXXZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_ZZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_X_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXXX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXXY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXXZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_ZZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Y_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXXX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XXXX(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXXY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XXXY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXXZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XXXZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XXYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XXYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XXZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_XZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_YYYY(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_YYYZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_YYZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_YZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_ZZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveQuadrupolePG_Z_ZZZZ(buffer_xx,
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

            t2cfunc::distribute(matrix_xx, buffer_xx, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy, buffer_xy, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz, buffer_xz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy, buffer_yy, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz, buffer_yz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz, buffer_zz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace mpol
