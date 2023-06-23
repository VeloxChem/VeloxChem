#include "OverlapRecGP.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGP(CSubMatrix*      matrix,
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

            // compute primitive integrals block (XXXX)

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

                    ovlrec::compPrimitiveOverlapGP_XXXX_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 3.0, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 3.0, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY)

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

                    ovlrec::compPrimitiveOverlapGP_XXXY_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ)

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

                    ovlrec::compPrimitiveOverlapGP_XXXZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, f4_17, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, f4_17, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY)

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

                    ovlrec::compPrimitiveOverlapGP_XXYY_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 6.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 6.0, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 6.0, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ)

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

                    ovlrec::compPrimitiveOverlapGP_XXYZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ)

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

                    ovlrec::compPrimitiveOverlapGP_XXZZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, -24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -24.0, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -24.0, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY)

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

                    ovlrec::compPrimitiveOverlapGP_XYYY_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ)

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

                    ovlrec::compPrimitiveOverlapGP_XYYZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ)

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

                    ovlrec::compPrimitiveOverlapGP_XYZZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ)

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

                    ovlrec::compPrimitiveOverlapGP_XZZZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY)

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

                    ovlrec::compPrimitiveOverlapGP_YYYY_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 3.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 3.0, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 3.0, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ)

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

                    ovlrec::compPrimitiveOverlapGP_YYYZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ)

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

                    ovlrec::compPrimitiveOverlapGP_YYZZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, -24.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -24.0, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -24.0, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ)

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

                    ovlrec::compPrimitiveOverlapGP_YZZZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ)

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

                    ovlrec::compPrimitiveOverlapGP_ZZZZ_T(buffer_x,
                                                          buffer_y,
                                                          buffer_z,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, 8.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 8.0, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 8.0, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);
        }
    }
}

auto
compPrimitiveOverlapGP_XXXX_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_x + 2.0 * fe_0 * rpa_x * rpa_x * rpa_x + 3.0 * fe_0 * fe_0 * rpa_x +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_x * rpa_x * rpa_x * rpa_x * rpb_x);

        fints_y[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_x * rpa_x * rpa_x * rpa_x * rpb_y);

        fints_z[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_x * rpa_x * rpa_x * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XXXY_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_y * rpa_x * rpa_x * rpa_x * rpb_x);

        fints_y[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_y * rpa_x * rpa_x * rpa_x * rpb_y);

        fints_z[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z + rpa_y * rpa_x * rpa_x * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XXXZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_x * rpa_x * rpa_x * rpb_x);

        fints_y[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y + rpa_z * rpa_x * rpa_x * rpa_x * rpb_y);

        fints_z[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_z * rpa_x * rpa_x * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XXYY_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * (fe_0 * rpa_y * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x +
                             (1.0 / 2.0) * fe_0 * fe_0 * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x);

        fints_x[i] += fss * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x;

        fints_y[i] += fss * (fe_0 * rpa_y * rpa_x * rpa_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y +
                             (1.0 / 2.0) * fe_0 * fe_0 * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y);

        fints_y[i] += fss * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y;

        fints_z[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_y * rpa_y * rpa_x * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XXYZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_x);

        fints_y[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_y * rpa_x * rpa_x * rpb_y);

        fints_z[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XXZZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * (fe_0 * rpa_z * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x +
                             (1.0 / 2.0) * fe_0 * fe_0 * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x);

        fints_x[i] += fss * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x;

        fints_y[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpa_x * rpa_x * rpb_y);

        fints_z[i] += fss * (fe_0 * rpa_z * rpa_x * rpa_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z +
                             (1.0 / 2.0) * fe_0 * fe_0 * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z);

        fints_z[i] += fss * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z;
    }
}

auto
compPrimitiveOverlapGP_XYYY_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_y * rpa_y * rpa_y * rpa_x * rpb_x);

        fints_y[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_y * rpa_y * rpa_y * rpa_x * rpb_y);

        fints_z[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z + rpa_y * rpa_y * rpa_y * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XYYZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_y * rpa_y * rpa_x * rpb_x);

        fints_y[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y + rpa_z * rpa_y * rpa_y * rpa_x * rpb_y);

        fints_z[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_z * rpa_y * rpa_y * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XYZZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y * rpa_x * rpb_x);

        fints_y[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_z * rpa_z * rpa_y * rpa_x * rpb_y);

        fints_z[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z + rpa_z * rpa_z * rpa_y * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_XZZZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z * rpa_x * rpb_x);

        fints_y[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y + rpa_z * rpa_z * rpa_z * rpa_x * rpb_y);

        fints_z[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_x + rpa_z * rpa_z * rpa_z * rpa_x * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_YYYY_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpa_y * rpa_y * rpb_x);

        fints_y[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y + 2.0 * fe_0 * rpa_y * rpa_y * rpa_y + 3.0 * fe_0 * fe_0 * rpa_y +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpa_y * rpb_y);

        fints_z[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_y * rpa_y * rpa_y * rpa_y * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_YYYZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_y * rpa_y * rpa_y * rpb_x);

        fints_y[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_y * rpa_y * rpa_y * rpb_y);

        fints_z[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_y * rpa_y * rpa_y * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_YYZZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x +
                             (1.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_x);

        fints_y[i] += fss * (fe_0 * rpa_z * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y +
                             (1.0 / 2.0) * fe_0 * fe_0 * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y);

        fints_y[i] += fss * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y;

        fints_z[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z +
                             (1.0 / 2.0) * fe_0 * fe_0 * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z);

        fints_z[i] += fss * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z;
    }
}

auto
compPrimitiveOverlapGP_YZZZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_y[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z * rpa_y * rpb_y);

        fints_z[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpa_y + rpa_z * rpa_z * rpa_z * rpa_y * rpb_z);
    }
}

auto
compPrimitiveOverlapGP_ZZZZ_T(TDoubleArray&       buffer_x,
                              TDoubleArray&       buffer_y,
                              TDoubleArray&       buffer_z,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_z * rpa_z * rpb_x);

        fints_y[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpa_z * rpa_z * rpb_y);

        fints_z[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_z + 2.0 * fe_0 * rpa_z * rpa_z * rpa_z + 3.0 * fe_0 * fe_0 * rpa_z +
                             (3.0 / 4.0) * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_z * rpa_z * rpb_z);
    }
}

}  // namespace ovlrec
