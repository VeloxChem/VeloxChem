#include "OverlapGeom101RecPP.hpp"

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom101PP_X_X.hpp"
#include "PrimitiveOverlapGeom101PP_X_Y.hpp"
#include "PrimitiveOverlapGeom101PP_X_Z.hpp"
#include "PrimitiveOverlapGeom101PP_Y_X.hpp"
#include "PrimitiveOverlapGeom101PP_Y_Y.hpp"
#include "PrimitiveOverlapGeom101PP_Y_Z.hpp"
#include "PrimitiveOverlapGeom101PP_Z_X.hpp"
#include "PrimitiveOverlapGeom101PP_Z_Y.hpp"
#include "PrimitiveOverlapGeom101PP_Z_Z.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom101PP(CSubMatrix*      matrix_x_x,
                     CSubMatrix*      matrix_x_y,
                     CSubMatrix*      matrix_x_z,
                     CSubMatrix*      matrix_y_x,
                     CSubMatrix*      matrix_y_y,
                     CSubMatrix*      matrix_y_z,
                     CSubMatrix*      matrix_z_x,
                     CSubMatrix*      matrix_z_y,
                     CSubMatrix*      matrix_z_z,
                     const CGtoBlock& bra_gto_block,
                     const CGtoBlock& ket_gto_block,
                     const int64_t    bra_first,
                     const int64_t    bra_last,
                     const mat_t      mat_type) -> void
{
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

    alignas(64) TDoubleArray buffer_x_x;

    alignas(64) TDoubleArray buffer_x_y;

    alignas(64) TDoubleArray buffer_x_z;

    alignas(64) TDoubleArray buffer_y_x;

    alignas(64) TDoubleArray buffer_y_y;

    alignas(64) TDoubleArray buffer_y_z;

    alignas(64) TDoubleArray buffer_z_x;

    alignas(64) TDoubleArray buffer_z_y;

    alignas(64) TDoubleArray buffer_z_z;

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

            // compute primitive integrals block (X_X)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_X_X(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (X_Y)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_X_Y(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (X_Z)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_X_Z(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Y_X)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_Y_X(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Y_Y)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_Y_Y(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Y_Z)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_Y_Z(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Z_X)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_Z_X(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Z_Y)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_Z_Y(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Z_Z)

            simd::zero(buffer_x_x);

            simd::zero(buffer_x_y);

            simd::zero(buffer_x_z);

            simd::zero(buffer_y_x);

            simd::zero(buffer_y_y);

            simd::zero(buffer_y_z);

            simd::zero(buffer_z_x);

            simd::zero(buffer_z_y);

            simd::zero(buffer_z_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom101PP_Z_Z(buffer_x_x,
                                                              buffer_x_y,
                                                              buffer_x_z,
                                                              buffer_y_x,
                                                              buffer_y_y,
                                                              buffer_y_z,
                                                              buffer_z_x,
                                                              buffer_z_y,
                                                              buffer_z_z,
                                                              bra_exp,
                                                              bra_norm,
                                                              bra_coord,
                                                              ket_exps,
                                                              ket_norms,
                                                              ket_coords_x,
                                                              ket_coords_y,
                                                              ket_coords_z,
                                                              ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x_x, buffer_x_x, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_y, buffer_x_y, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x_z, buffer_x_z, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_x, buffer_y_x, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_y, buffer_y_y, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y_z, buffer_y_z, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_x, buffer_z_x, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_y, buffer_z_y, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z_z, buffer_z_z, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace ovlrec
