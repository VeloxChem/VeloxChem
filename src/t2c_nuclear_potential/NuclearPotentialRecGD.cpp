#include "NuclearPotentialRecGD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveNuclearPotentialGD_XXXX_T.hpp"
#include "PrimitiveNuclearPotentialGD_XXXY_T.hpp"
#include "PrimitiveNuclearPotentialGD_XXXZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_XXYY_T.hpp"
#include "PrimitiveNuclearPotentialGD_XXYZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_XXZZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_XYYY_T.hpp"
#include "PrimitiveNuclearPotentialGD_XYYZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_XYZZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_XZZZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_YYYY_T.hpp"
#include "PrimitiveNuclearPotentialGD_YYYZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_YYZZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_YZZZ_T.hpp"
#include "PrimitiveNuclearPotentialGD_ZZZZ_T.hpp"
#include "T2CDistributor.hpp"

namespace npotrec {  // npotrec namespace

auto
compNuclearPotentialGD(CSubMatrix*      matrix,
                       const double     charge,
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

            // compute primitive integrals block (XXXX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XXXX_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XXXY_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XXXZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XXYY_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -6.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 1.50 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 6.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -1.50 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XXYZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XXZZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XYYY_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XYYZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 3.0 * f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XYZZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -6.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 6.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_XZZZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_YYYY_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -3.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -0.5 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -0.25 * f4_35 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 3.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 0.5 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 0.25 * f4_35 * f2_3, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_YYYZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, f4_17 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 3.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -f4_17 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_YYZZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, -24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, -3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, 24.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, 3.0 * f4_5 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -24.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, -3.0 * f4_5 * f2_3, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_YZZZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xx, 4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yy, -4.0 * f4_2 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 4.0 * f4_2 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialGD_ZZZZ_T(buffer_xx,
                                                                    buffer_xy,
                                                                    buffer_xz,
                                                                    buffer_yy,
                                                                    buffer_yz,
                                                                    buffer_zz,
                                                                    charge,
                                                                    point,
                                                                    bra_exp,
                                                                    bra_norm,
                                                                    bra_coord,
                                                                    ket_exps,
                                                                    ket_norms,
                                                                    ket_coords_x,
                                                                    ket_coords_y,
                                                                    ket_coords_z,
                                                                    ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xx, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xx, 8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xy, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yy, -8.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yz, 8.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace npotrec
