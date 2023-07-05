#include "NuclearPotentialRecFP.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveNuclearPotentialFP_XXX_T.hpp"
#include "PrimitiveNuclearPotentialFP_XXY_T.hpp"
#include "PrimitiveNuclearPotentialFP_XXZ_T.hpp"
#include "PrimitiveNuclearPotentialFP_XYY_T.hpp"
#include "PrimitiveNuclearPotentialFP_XYZ_T.hpp"
#include "PrimitiveNuclearPotentialFP_XZZ_T.hpp"
#include "PrimitiveNuclearPotentialFP_YYY_T.hpp"
#include "PrimitiveNuclearPotentialFP_YYZ_T.hpp"
#include "PrimitiveNuclearPotentialFP_YZZ_T.hpp"
#include "PrimitiveNuclearPotentialFP_ZZZ_T.hpp"
#include "T2CDistributor.hpp"

namespace npotrec {  // npotrec namespace

auto
compNuclearPotentialFP(CSubMatrix*      matrix,
                       const double     charge,
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

            // compute primitive integrals block (XXX)

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

                    npotrec::compPrimitiveNuclearPotentialFP_XXX_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, f3_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY)

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

                    npotrec::compPrimitiveNuclearPotentialFP_XXY_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ)

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

                    npotrec::compPrimitiveNuclearPotentialFP_XXZ_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY)

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

                    npotrec::compPrimitiveNuclearPotentialFP_XYY_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ)

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

                    npotrec::compPrimitiveNuclearPotentialFP_XYZ_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ)

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

                    npotrec::compPrimitiveNuclearPotentialFP_XZZ_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY)

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

                    npotrec::compPrimitiveNuclearPotentialFP_YYY_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ)

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

                    npotrec::compPrimitiveNuclearPotentialFP_YYZ_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ)

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

                    npotrec::compPrimitiveNuclearPotentialFP_YZZ_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ)

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

                    npotrec::compPrimitiveNuclearPotentialFP_ZZZ_T(buffer_x,
                                                                   buffer_y,
                                                                   buffer_z,
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

            t2cfunc::distribute(matrix, buffer_x, 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_y, 2.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_z, 2.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace npotrec
