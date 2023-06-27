#include "NuclearPotentialGeom010RecFP.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "BoysFunc.hpp"
#include "T2CDistributor.hpp"

namespace geom_npotrec { // geom_npotrec namespace

auto
compNuclearPotentialGeom010FP(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  bra_gto_block,
                              const CGtoBlock&  ket_gto_block,
                              const bool        ang_order,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void
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

        simd::loadCoordinates(ket_coords_x,
                              ket_coords_y,
                              ket_coords_z,
                              ket_gto_coords,
                              ket_first,
                              ket_last);

        for (int64_t j = bra_first; j < bra_last; j++) 
        {
            const auto bra_coord = bra_gto_coords[j];

            // compute primitive integrals block (XXX_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXX_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXX_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXX_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXX_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXY_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXY_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXY_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXZ_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXZ_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XXZ_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XYY_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XYY_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XYY_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XYZ_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XYZ_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XYZ_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XZZ_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XZZ_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_XZZ_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YYY_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YYY_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YYY_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YYZ_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YYZ_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YYZ_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YZZ_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YZZ_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_YZZ_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_X)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_ZZZ_X(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_Y)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_ZZZ_Y(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_Z)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FP_ZZZ_Z(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, ang_order);

        }
    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXX_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x + (3.0 / 2.0) * fe_0 * rpb_x + 3.0 * rpa_x * rpa_x * rpb_x + rpa_x * rpa_x * rpa_x);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpb_x - 6.0 * fe_0 * rpc_x - 6.0 * rpa_x * rpb_x * rpc_x - 6.0 * rpa_x * rpa_x * rpc_x);

        fints_x[i] += dip_x * fss * b3_vals[i] * (6.0 * fe_0 * rpc_x + 9.0 * rpa_x * rpc_x * rpc_x + 3.0 * rpb_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-4.0 * rpc_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-3.0 * rpa_x * rpa_x * rpb_x * rpc_x - rpa_x * rpa_x * rpa_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_x + 3.0 * fe_0 * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * 3.0 * rpa_x * rpa_x * rpc_x * rpc_x;

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-3.0 * fe_0 * rpc_x * rpc_x - 3.0 * rpa_x * rpc_x * rpc_x * rpc_x - rpb_x * rpc_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-3.0 * rpa_x * rpa_x * rpb_x * rpc_x - rpa_x * rpa_x * rpa_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_x + 3.0 * fe_0 * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * 3.0 * rpa_x * rpa_x * rpc_x * rpc_x;

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-3.0 * fe_0 * rpc_x * rpc_x - 3.0 * rpa_x * rpc_x * rpc_x * rpc_x - rpb_x * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-3.0 * rpa_x * rpa_x * rpb_x * rpc_x - rpa_x * rpa_x * rpa_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_x + 3.0 * fe_0 * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * 3.0 * rpa_x * rpa_x * rpc_x * rpc_x;

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-3.0 * fe_0 * rpc_x * rpc_x - 3.0 * rpa_x * rpc_x * rpc_x * rpc_x - rpb_x * rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_x * rpc_x * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXX_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y + 3.0 * rpa_x * rpa_x * rpb_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y - (3.0 / 2.0) * fe_0 * rpc_y - 6.0 * rpa_x * rpb_y * rpc_x - 3.0 * rpa_x * rpa_x * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 6.0 * rpa_x * rpc_y * rpc_x + 3.0 * rpb_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-3.0 * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_x * rpa_x * rpa_x * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 3.0 * rpa_x * rpc_y * rpc_x * rpc_x - rpb_y * rpc_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x + rpa_x * rpa_x * rpa_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpc_x - 3.0 * rpa_x * rpa_x * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 3.0 * rpa_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_x * rpa_x * rpa_x * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 3.0 * rpa_x * rpc_y * rpc_x * rpc_x - rpb_y * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_x * rpa_x * rpa_x * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 3.0 * rpa_x * rpc_y * rpc_x * rpc_x - rpb_y * rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXX_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z + 3.0 * rpa_x * rpa_x * rpb_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z - (3.0 / 2.0) * fe_0 * rpc_z - 6.0 * rpa_x * rpb_z * rpc_x - 3.0 * rpa_x * rpa_x * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 6.0 * rpa_x * rpc_z * rpc_x + 3.0 * rpb_z * rpc_x * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_x * rpa_x * rpa_x * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 3.0 * rpa_x * rpc_z * rpc_x * rpc_x - rpb_z * rpc_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_x * rpa_x * rpa_x * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 3.0 * rpa_x * rpc_z * rpc_x * rpc_x - rpb_z * rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x + rpa_x * rpa_x * rpa_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpc_x - 3.0 * rpa_x * rpa_x * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 3.0 * rpa_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_x * rpa_x * rpa_x * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 3.0 * rpa_x * rpc_z * rpc_x * rpc_x - rpb_z * rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXY_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y + 2.0 * rpa_y * rpa_x * rpb_x + rpa_y * rpa_x * rpa_x);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpc_y - 4.0 * rpa_y * rpa_x * rpc_x - 2.0 * rpa_y * rpb_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-rpa_x * rpa_x * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 3.0 * rpa_y * rpc_x * rpc_x + 4.0 * rpa_x * rpc_y * rpc_x + 2.0 * rpb_x * rpc_y * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-3.0 * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * (fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_y * rpa_x * rpa_x * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_x - fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_x - rpa_x * rpa_x * rpb_x * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_x + fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - rpa_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_x * rpc_x - rpb_x * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * (fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * rpb_x + rpa_x * rpa_x * rpb_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpb_x - (3.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_x * rpb_x * rpc_x - rpa_x * rpa_x * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_x * rpc_x * rpc_x + rpb_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * (fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_y * rpa_x * rpa_x * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_x - fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_x - rpa_x * rpa_x * rpb_x * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_x + fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - rpa_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_x * rpc_x - rpb_x * rpc_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * (fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_y * rpa_x * rpa_x * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_x - fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_x - rpa_x * rpa_x * rpb_x * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_x + fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - rpa_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_x * rpc_x - rpb_x * rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXY_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * (fe_0 * rpa_x + 2.0 * rpa_y * rpa_x * rpb_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-fe_0 * rpa_x - fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpc_y - 2.0 * rpa_y * rpb_y * rpc_x - 2.0 * rpa_x * rpb_y * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * (fe_0 * rpc_x + 2.0 * rpa_y * rpc_y * rpc_x + 2.0 * rpa_x * rpc_y * rpc_y + 2.0 * rpb_y * rpc_y * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-2.0 * rpc_y * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y - fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_y - rpa_x * rpa_x * rpb_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_y + fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_y * rpa_x * rpc_y * rpc_x + rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - rpa_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_y * rpc_x - rpb_y * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * rpb_y + rpa_y * rpa_x * rpa_x + rpa_x * rpa_x * rpb_y);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpb_y - fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpc_x - 2.0 * rpa_x * rpb_y * rpc_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-2.0 * rpa_x * rpa_x * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * (fe_0 * rpc_y + rpa_y * rpc_x * rpc_x + 4.0 * rpa_x * rpc_y * rpc_x + rpb_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-2.0 * rpc_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y - fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_y - rpa_x * rpa_x * rpb_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_y + fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_y * rpa_x * rpc_y * rpc_x + rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - rpa_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_y * rpc_x - rpb_y * rpc_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y - fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_y - rpa_x * rpa_x * rpb_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_y + fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_y * rpa_x * rpc_y * rpc_x + rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - rpa_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_y * rpc_x - rpb_y * rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXY_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * 2.0 * rpa_y * rpa_x * rpb_z;

        fints_x[i] += dip_x * fss * b2_vals[i] * (-2.0 * rpa_y * rpa_x * rpc_z - 2.0 * rpa_y * rpb_z * rpc_x - 2.0 * rpa_x * rpb_z * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * (2.0 * rpa_y * rpc_z * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y + 2.0 * rpb_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_y * rpa_x * rpa_x * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_x * rpa_x * rpb_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpc_z * rpc_x + rpa_y * rpb_z * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (2.0 * rpa_x * rpb_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x - rpb_z * rpc_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z + rpa_x * rpa_x * rpb_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_x * rpb_z * rpc_x - rpa_x * rpa_x * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + 2.0 * rpa_x * rpc_z * rpc_x + rpb_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_y * rpa_x * rpa_x * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_x * rpa_x * rpb_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpc_z * rpc_x + rpa_y * rpb_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (2.0 * rpa_x * rpb_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x - rpb_z * rpc_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_x * rpa_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpc_x - rpa_x * rpa_x * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + rpa_y * rpc_x * rpc_x + 2.0 * rpa_x * rpc_y * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_y * rpa_x * rpa_x * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_x * rpa_x * rpb_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpc_z * rpc_x + rpa_y * rpb_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (2.0 * rpa_x * rpb_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x - rpb_z * rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXZ_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z + 2.0 * rpa_z * rpa_x * rpb_x + rpa_z * rpa_x * rpa_x);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpc_z - 4.0 * rpa_z * rpa_x * rpc_x - 2.0 * rpa_z * rpb_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-rpa_x * rpa_x * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 3.0 * rpa_z * rpc_x * rpc_x + 4.0 * rpa_x * rpc_z * rpc_x + 2.0 * rpb_x * rpc_z * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * (fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_x * rpa_x * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_x - fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_x + fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_x * rpc_x - rpb_x * rpc_z * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * (fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_x * rpa_x * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_x - fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_x + fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_x * rpc_x - rpb_x * rpc_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * (fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * rpb_x + rpa_x * rpa_x * rpb_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpb_x - (3.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_x * rpb_x * rpc_x - rpa_x * rpa_x * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_x * rpc_x * rpc_x + rpb_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * (fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_x * rpa_x * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_x - fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_x + fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_x * rpc_x - rpb_x * rpc_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXZ_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * 2.0 * rpa_z * rpa_x * rpb_y;

        fints_x[i] += dip_x * fss * b2_vals[i] * (-2.0 * rpa_z * rpa_x * rpc_y - 2.0 * rpa_z * rpb_y * rpc_x - 2.0 * rpa_x * rpb_y * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * (2.0 * rpa_z * rpc_y * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y + 2.0 * rpb_y * rpc_z * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_x * rpa_x * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_x * rpb_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_x * rpa_x * rpb_y * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpc_y * rpc_x + rpa_z * rpb_y * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (2.0 * rpa_x * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x - rpb_y * rpc_z * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_x * rpa_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpc_x - rpa_x * rpa_x * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_x * rpa_x * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_x * rpb_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_x * rpa_x * rpb_y * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpc_y * rpc_x + rpa_z * rpb_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (2.0 * rpa_x * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x - rpb_y * rpc_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y + rpa_x * rpa_x * rpb_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_x * rpb_y * rpc_x - rpa_x * rpa_x * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_x * rpc_y * rpc_x + rpb_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_x * rpa_x * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_x * rpb_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_x * rpa_x * rpb_y * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpc_y * rpc_x + rpa_z * rpb_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (2.0 * rpa_x * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x - rpb_y * rpc_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XXZ_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * (fe_0 * rpa_x + 2.0 * rpa_z * rpa_x * rpb_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-fe_0 * rpa_x - fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpc_z - 2.0 * rpa_z * rpb_z * rpc_x - 2.0 * rpa_x * rpb_z * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * (fe_0 * rpc_x + 2.0 * rpa_z * rpc_z * rpc_x + 2.0 * rpa_x * rpc_z * rpc_z + 2.0 * rpb_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z - fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z - rpa_x * rpa_x * rpb_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z + fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_x * rpc_z * rpc_x + rpa_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - rpa_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_x - rpb_z * rpc_z * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z - fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z - rpa_x * rpa_x * rpb_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z + fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_x * rpc_z * rpc_x + rpa_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - rpa_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_x - rpb_z * rpc_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * rpb_z + rpa_z * rpa_x * rpa_x + rpa_x * rpa_x * rpb_z);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpb_z - fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpc_x - 2.0 * rpa_x * rpb_z * rpc_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_x * rpa_x * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * (fe_0 * rpc_z + rpa_z * rpc_x * rpc_x + 4.0 * rpa_x * rpc_z * rpc_x + rpb_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z - fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z - rpa_x * rpa_x * rpb_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z + fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_x * rpc_z * rpc_x + rpa_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - rpa_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_x - rpb_z * rpc_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XYY_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * rpb_x + rpa_y * rpa_y * rpa_x + rpa_y * rpa_y * rpb_x);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpb_x - fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpc_y - 2.0 * rpa_y * rpb_x * rpc_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-2.0 * rpa_y * rpa_y * rpc_x);

        fints_x[i] += dip_x * fss * b3_vals[i] * (fe_0 * rpc_x + 4.0 * rpa_y * rpc_y * rpc_x + rpa_x * rpc_y * rpc_y + rpb_x * rpc_y * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-2.0 * rpc_y * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_x * rpc_y - rpa_y * rpa_y * rpa_x * rpc_x - rpa_y * rpa_y * rpb_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_y * rpa_x * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpc_x * rpc_x + rpa_x * rpb_x * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpa_y * rpc_y * rpc_x * rpc_x - rpa_x * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * (fe_0 * rpa_y + 2.0 * rpa_y * rpa_x * rpb_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-fe_0 * rpa_y - fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpc_x - 2.0 * rpa_y * rpb_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * (fe_0 * rpc_y + 2.0 * rpa_y * rpc_x * rpc_x + 2.0 * rpa_x * rpc_y * rpc_x + 2.0 * rpb_x * rpc_y * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-2.0 * rpc_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_x * rpc_y - rpa_y * rpa_y * rpa_x * rpc_x - rpa_y * rpa_y * rpb_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_y * rpa_x * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpc_x * rpc_x + rpa_x * rpb_x * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpa_y * rpc_y * rpc_x * rpc_x - rpa_x * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_x * rpc_y - rpa_y * rpa_y * rpa_x * rpc_x - rpa_y * rpa_y * rpb_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_y * rpa_x * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpc_x * rpc_x + rpa_x * rpb_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpa_y * rpc_y * rpc_x * rpc_x - rpa_x * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XYY_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * (fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * rpb_y + rpa_y * rpa_y * rpb_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpb_y - (3.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_y * rpb_y * rpc_y - rpa_y * rpa_y * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_y * rpc_y * rpc_y + rpb_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * (fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_y * rpa_y * rpa_x * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_y * rpa_x - fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_y * rpc_y - rpa_y * rpa_y * rpa_x * rpc_y - rpa_y * rpa_y * rpb_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (2.0 * rpa_y * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpc_y * rpc_x + rpa_x * rpb_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_y * rpc_y * rpc_y * rpc_x - rpa_x * rpc_y * rpc_y * rpc_y - rpb_y * rpc_y * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x + 2.0 * rpa_y * rpa_x * rpb_y + rpa_y * rpa_y * rpa_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpc_x - 4.0 * rpa_y * rpa_x * rpc_y - 2.0 * rpa_y * rpb_y * rpc_x - rpa_y * rpa_y * rpc_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-2.0 * rpa_x * rpb_y * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 4.0 * rpa_y * rpc_y * rpc_x + 3.0 * rpa_x * rpc_y * rpc_y + 2.0 * rpb_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-3.0 * rpc_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * (fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_y * rpa_y * rpa_x * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_y * rpa_x - fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_y * rpc_y - rpa_y * rpa_y * rpa_x * rpc_y - rpa_y * rpa_y * rpb_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (2.0 * rpa_y * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpc_y * rpc_x + rpa_x * rpb_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_y * rpc_y * rpc_y * rpc_x - rpa_x * rpc_y * rpc_y * rpc_y - rpb_y * rpc_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * (fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_y * rpa_y * rpa_x * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_y * rpa_x - fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_y * rpc_y - rpa_y * rpa_y * rpa_x * rpc_y - rpa_y * rpa_y * rpb_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (2.0 * rpa_y * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpc_y * rpc_x + rpa_x * rpb_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_y * rpc_y * rpc_y * rpc_x - rpa_x * rpc_y * rpc_y * rpc_y - rpb_y * rpc_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XYY_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z + rpa_y * rpa_y * rpb_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_y * rpb_z * rpc_y - rpa_y * rpa_y * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + 2.0 * rpa_y * rpc_z * rpc_y + rpb_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_y * rpa_y * rpa_x * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpc_y - rpa_y * rpa_y * rpa_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_y * rpa_y * rpb_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_x + rpa_x * rpb_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * 2.0 * rpa_y * rpa_x * rpb_z;

        fints_y[i] += dip_y * fss * b2_vals[i] * (-2.0 * rpa_y * rpa_x * rpc_z - 2.0 * rpa_y * rpb_z * rpc_x - 2.0 * rpa_x * rpb_z * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * (2.0 * rpa_y * rpc_z * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y + 2.0 * rpb_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_y * rpa_y * rpa_x * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpc_y - rpa_y * rpa_y * rpa_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_y * rpa_y * rpb_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_x + rpa_x * rpb_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_y * rpa_y * rpa_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpc_y - rpa_y * rpa_y * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_y * rpc_y * rpc_x + rpa_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_y * rpa_y * rpa_x * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpc_y - rpa_y * rpa_y * rpa_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_y * rpa_y * rpb_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_x + rpa_x * rpb_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XYZ_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * (rpa_z * rpa_y * rpa_x + rpa_z * rpa_y * rpb_x);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-2.0 * rpa_z * rpa_y * rpc_x - rpa_z * rpa_x * rpc_y - rpa_z * rpb_x * rpc_y - rpa_y * rpa_x * rpc_z - rpa_y * rpb_x * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * (2.0 * rpa_z * rpc_y * rpc_x + 2.0 * rpa_y * rpc_z * rpc_x + rpa_x * rpc_z * rpc_y + rpb_x * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y + rpa_z * rpa_y * rpa_x * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z - rpa_z * rpa_y * rpa_x * rpc_x - rpa_z * rpa_y * rpb_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_z * rpa_x * rpb_x * rpc_y - rpa_y * rpa_x * rpb_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + rpa_z * rpa_y * rpc_x * rpc_x + rpa_z * rpa_x * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpb_x * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_x + rpa_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_x * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_x * rpb_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - rpa_z * rpa_x * rpc_x - rpa_z * rpb_x * rpc_x - rpa_x * rpb_x * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_x * rpc_x + rpa_x * rpc_z * rpc_x + rpb_x * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y + rpa_z * rpa_y * rpa_x * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z - rpa_z * rpa_y * rpa_x * rpc_x - rpa_z * rpa_y * rpb_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_z * rpa_x * rpb_x * rpc_y - rpa_y * rpa_x * rpb_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + rpa_z * rpa_y * rpc_x * rpc_x + rpa_z * rpa_x * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpb_x * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_x + rpa_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_x * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_x * rpb_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - rpa_y * rpa_x * rpc_x - rpa_y * rpb_x * rpc_x - rpa_x * rpb_x * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + rpa_y * rpc_x * rpc_x + rpa_x * rpc_y * rpc_x + rpb_x * rpc_y * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y + rpa_z * rpa_y * rpa_x * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z - rpa_z * rpa_y * rpa_x * rpc_x - rpa_z * rpa_y * rpb_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_z * rpa_x * rpb_x * rpc_y - rpa_y * rpa_x * rpb_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + rpa_z * rpa_y * rpc_x * rpc_x + rpa_z * rpa_x * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpb_x * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_x + rpa_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XYZ_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_y * rpb_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - rpa_z * rpa_y * rpc_y - rpa_z * rpb_y * rpc_y - rpa_y * rpb_y * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_y * rpc_y + rpa_y * rpc_z * rpc_y + rpb_y * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x + rpa_z * rpa_y * rpa_x * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z - rpa_z * rpa_y * rpa_x * rpc_y - rpa_z * rpa_y * rpb_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_z * rpa_x * rpb_y * rpc_y - rpa_y * rpa_x * rpb_y * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + rpa_z * rpa_y * rpc_y * rpc_x + rpa_z * rpa_x * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpb_y * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_y + rpa_y * rpb_y * rpc_z * rpc_x + rpa_x * rpb_y * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y - rpb_y * rpc_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * (rpa_z * rpa_y * rpa_x + rpa_z * rpa_x * rpb_y);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-rpa_z * rpa_y * rpc_x - 2.0 * rpa_z * rpa_x * rpc_y - rpa_z * rpb_y * rpc_x - rpa_y * rpa_x * rpc_z - rpa_x * rpb_y * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * (2.0 * rpa_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y + rpb_y * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x + rpa_z * rpa_y * rpa_x * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z - rpa_z * rpa_y * rpa_x * rpc_y - rpa_z * rpa_y * rpb_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_z * rpa_x * rpb_y * rpc_y - rpa_y * rpa_x * rpb_y * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + rpa_z * rpa_y * rpc_y * rpc_x + rpa_z * rpa_x * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpb_y * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_y + rpa_y * rpb_y * rpc_z * rpc_x + rpa_x * rpb_y * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y - rpb_y * rpc_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_y * rpa_x * rpb_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - rpa_y * rpa_x * rpc_y - rpa_y * rpb_y * rpc_x - rpa_x * rpb_y * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + rpa_y * rpc_y * rpc_x + rpa_x * rpc_y * rpc_y + rpb_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x + rpa_z * rpa_y * rpa_x * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z - rpa_z * rpa_y * rpa_x * rpc_y - rpa_z * rpa_y * rpb_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_z * rpa_x * rpb_y * rpc_y - rpa_y * rpa_x * rpb_y * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + rpa_z * rpa_y * rpc_y * rpc_x + rpa_z * rpa_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpb_y * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_y + rpa_y * rpb_y * rpc_z * rpc_x + rpa_x * rpb_y * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y - rpb_y * rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XYZ_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_z * rpa_y * rpb_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - rpa_z * rpa_y * rpc_z - rpa_z * rpb_z * rpc_y - rpa_y * rpb_z * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + rpa_z * rpc_z * rpc_y + rpa_y * rpc_z * rpc_z + rpb_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x + rpa_z * rpa_y * rpa_x * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y - rpa_z * rpa_y * rpa_x * rpc_z - rpa_z * rpa_y * rpb_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_z * rpa_x * rpb_z * rpc_y - rpa_y * rpa_x * rpb_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpc_z * rpc_x + rpa_z * rpa_x * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpb_z * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_z + rpa_y * rpb_z * rpc_z * rpc_x + rpa_x * rpb_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - rpa_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y - rpb_z * rpc_z * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_z * rpa_x * rpb_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - rpa_z * rpa_x * rpc_z - rpa_z * rpb_z * rpc_x - rpa_x * rpb_z * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + rpa_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z + rpb_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x + rpa_z * rpa_y * rpa_x * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y - rpa_z * rpa_y * rpa_x * rpc_z - rpa_z * rpa_y * rpb_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_z * rpa_x * rpb_z * rpc_y - rpa_y * rpa_x * rpb_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpc_z * rpc_x + rpa_z * rpa_x * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpb_z * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_z + rpa_y * rpb_z * rpc_z * rpc_x + rpa_x * rpb_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - rpa_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y - rpb_z * rpc_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * (rpa_z * rpa_y * rpa_x + rpa_y * rpa_x * rpb_z);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-rpa_z * rpa_y * rpc_x - rpa_z * rpa_x * rpc_y - 2.0 * rpa_y * rpa_x * rpc_z - rpa_y * rpb_z * rpc_x - rpa_x * rpb_z * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * (rpa_z * rpc_y * rpc_x + 2.0 * rpa_y * rpc_z * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y + rpb_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x + rpa_z * rpa_y * rpa_x * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y - rpa_z * rpa_y * rpa_x * rpc_z - rpa_z * rpa_y * rpb_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_z * rpa_x * rpb_z * rpc_y - rpa_y * rpa_x * rpb_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpc_z * rpc_x + rpa_z * rpa_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpb_z * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_z + rpa_y * rpb_z * rpc_z * rpc_x + rpa_x * rpb_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - rpa_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y - rpb_z * rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XZZ_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + (1.0 / 2.0) * fe_0 * rpb_x + rpa_z * rpa_z * rpa_x + rpa_z * rpa_z * rpb_x);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpb_x - fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpc_z - 2.0 * rpa_z * rpb_x * rpc_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-2.0 * rpa_z * rpa_z * rpc_x);

        fints_x[i] += dip_x * fss * b3_vals[i] * (fe_0 * rpc_x + 4.0 * rpa_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z + rpb_x * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z - rpa_z * rpa_z * rpa_x * rpc_x - rpa_z * rpa_z * rpb_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (fe_0 * rpa_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_x + rpa_z * rpa_z * rpc_x * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z - rpa_z * rpa_z * rpa_x * rpc_x - rpa_z * rpa_z * rpb_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (fe_0 * rpa_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_x + rpa_z * rpa_z * rpc_x * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * (fe_0 * rpa_z + 2.0 * rpa_z * rpa_x * rpb_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-fe_0 * rpa_z - fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpc_x - 2.0 * rpa_z * rpb_x * rpc_x - 2.0 * rpa_x * rpb_x * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * (fe_0 * rpc_z + 2.0 * rpa_z * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_x + 2.0 * rpb_x * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z - rpa_z * rpa_z * rpa_x * rpc_x - rpa_z * rpa_z * rpb_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (fe_0 * rpa_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_x + rpa_z * rpa_z * rpc_x * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XZZ_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y + rpa_z * rpa_z * rpb_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_z * rpb_y * rpc_z - rpa_z * rpa_z * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_z * rpc_z * rpc_y + rpb_y * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_z * rpa_z * rpa_x * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z - rpa_z * rpa_z * rpa_x * rpc_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_z * rpa_z * rpb_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_x + rpa_x * rpb_y * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_z * rpa_z * rpa_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpc_z - rpa_z * rpa_z * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_z * rpa_z * rpa_x * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z - rpa_z * rpa_z * rpa_x * rpc_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_z * rpa_z * rpb_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_x + rpa_x * rpb_y * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * 2.0 * rpa_z * rpa_x * rpb_y;

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_z * rpa_x * rpc_y - 2.0 * rpa_z * rpb_y * rpc_x - 2.0 * rpa_x * rpb_y * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * (2.0 * rpa_z * rpc_y * rpc_x + 2.0 * rpa_x * rpc_z * rpc_y + 2.0 * rpb_y * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y + rpa_z * rpa_z * rpa_x * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z - rpa_z * rpa_z * rpa_x * rpc_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_z * rpa_z * rpb_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_x + rpa_x * rpb_y * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_XZZ_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * (fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * rpb_z + rpa_z * rpa_z * rpb_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpb_z - (3.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpb_z * rpc_z - rpa_z * rpa_z * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 2.0 * rpa_z * rpc_z * rpc_z + rpb_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * (fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_z * rpa_z * rpa_x * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_z * rpa_x - fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_z * rpc_z - rpa_z * rpa_z * rpa_x * rpc_z - rpa_z * rpa_z * rpb_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (2.0 * rpa_z * rpb_z * rpc_z * rpc_x + rpa_z * rpa_z * rpc_z * rpc_x + rpa_x * rpb_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * (fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_z * rpa_z * rpa_x * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_z * rpa_x - fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_z * rpc_z - rpa_z * rpa_z * rpa_x * rpc_z - rpa_z * rpa_z * rpb_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (2.0 * rpa_z * rpb_z * rpc_z * rpc_x + rpa_z * rpa_z * rpc_z * rpc_x + rpa_x * rpb_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x + 2.0 * rpa_z * rpa_x * rpb_z + rpa_z * rpa_z * rpa_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpc_x - 4.0 * rpa_z * rpa_x * rpc_z - 2.0 * rpa_z * rpb_z * rpc_x - rpa_z * rpa_z * rpc_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_x * rpb_z * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 4.0 * rpa_z * rpc_z * rpc_x + 3.0 * rpa_x * rpc_z * rpc_z + 2.0 * rpb_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * (fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z + rpa_z * rpa_z * rpa_x * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_z * rpa_x - fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_z * rpc_z - rpa_z * rpa_z * rpa_x * rpc_z - rpa_z * rpa_z * rpb_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (2.0 * rpa_z * rpb_z * rpc_z * rpc_x + rpa_z * rpa_z * rpc_z * rpc_x + rpa_x * rpb_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YYY_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_y * rpa_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 3.0 * rpa_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_y * rpa_y * rpa_y * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpb_x * rpc_y - rpa_y * rpa_y * rpa_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 3.0 * rpa_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x + 3.0 * rpa_y * rpa_y * rpb_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x - (3.0 / 2.0) * fe_0 * rpc_x - 6.0 * rpa_y * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 6.0 * rpa_y * rpc_y * rpc_x + 3.0 * rpb_x * rpc_y * rpc_y);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-3.0 * rpc_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_y * rpa_y * rpa_y * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpb_x * rpc_y - rpa_y * rpa_y * rpa_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 3.0 * rpa_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_y * rpa_y * rpa_y * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpb_x * rpc_y - rpa_y * rpa_y * rpa_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 3.0 * rpa_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YYY_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y - (9.0 / 2.0) * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpc_y - rpa_y * rpa_y * rpa_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y + 3.0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * 3.0 * rpa_y * rpa_y * rpc_y * rpc_y;

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y - 3.0 * rpa_y * rpc_y * rpc_y * rpc_y - rpb_y * rpc_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y + (3.0 / 2.0) * fe_0 * rpb_y + 3.0 * rpa_y * rpa_y * rpb_y + rpa_y * rpa_y * rpa_y);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpb_y - 6.0 * fe_0 * rpc_y - 6.0 * rpa_y * rpb_y * rpc_y - 6.0 * rpa_y * rpa_y * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * (6.0 * fe_0 * rpc_y + 9.0 * rpa_y * rpc_y * rpc_y + 3.0 * rpb_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-4.0 * rpc_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y - (9.0 / 2.0) * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpc_y - rpa_y * rpa_y * rpa_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y + 3.0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * 3.0 * rpa_y * rpa_y * rpc_y * rpc_y;

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y - 3.0 * rpa_y * rpc_y * rpc_y * rpc_y - rpb_y * rpc_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y - (9.0 / 2.0) * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpc_y - rpa_y * rpa_y * rpa_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y + 3.0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * 3.0 * rpa_y * rpa_y * rpc_y * rpc_y;

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y - 3.0 * rpa_y * rpc_y * rpc_y * rpc_y - rpb_y * rpc_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YYY_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_y * rpa_y * rpa_y * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpc_y - rpa_y * rpa_y * rpa_y * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 3.0 * rpa_y * rpc_z * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z + 3.0 * rpa_y * rpa_y * rpb_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z - (3.0 / 2.0) * fe_0 * rpc_z - 6.0 * rpa_y * rpb_z * rpc_y - 3.0 * rpa_y * rpa_y * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 6.0 * rpa_y * rpc_z * rpc_y + 3.0 * rpb_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_y * rpa_y * rpa_y * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpc_y - rpa_y * rpa_y * rpa_y * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 3.0 * rpa_y * rpc_z * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_y * rpa_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 3.0 * rpa_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_y * rpa_y * rpa_y * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpc_y - rpa_y * rpa_y * rpa_y * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 3.0 * rpa_y * rpc_z * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YYZ_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_y * rpa_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpc_y - rpa_y * rpa_y * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_y * rpc_y + 2.0 * rpa_y * rpc_z * rpc_y);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_y * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_y * rpa_y * rpb_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpc_y * rpc_x + rpa_z * rpb_x * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (2.0 * rpa_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * 2.0 * rpa_z * rpa_y * rpb_x;

        fints_y[i] += dip_y * fss * b2_vals[i] * (-2.0 * rpa_z * rpa_y * rpc_x - 2.0 * rpa_z * rpb_x * rpc_y - 2.0 * rpa_y * rpb_x * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * (2.0 * rpa_z * rpc_y * rpc_x + 2.0 * rpa_y * rpc_z * rpc_x + 2.0 * rpb_x * rpc_z * rpc_y);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_y * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_y * rpa_y * rpb_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpc_y * rpc_x + rpa_z * rpb_x * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (2.0 * rpa_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x + rpa_y * rpa_y * rpb_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_y * rpb_x * rpc_y - rpa_y * rpa_y * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_y * rpc_y * rpc_x + rpb_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_y * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_y * rpa_y * rpb_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpc_y * rpc_x + rpa_z * rpb_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (2.0 * rpa_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YYZ_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * (fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_y * rpa_y * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y - fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_y * rpc_y - rpa_z * rpa_y * rpa_y * rpc_y - rpa_y * rpa_y * rpb_y * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y + fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y - rpb_y * rpc_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z + 2.0 * rpa_z * rpa_y * rpb_y + rpa_z * rpa_y * rpa_y);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpc_z - 4.0 * rpa_z * rpa_y * rpc_y - 2.0 * rpa_z * rpb_y * rpc_y - 2.0 * rpa_y * rpb_y * rpc_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-rpa_y * rpa_y * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 3.0 * rpa_z * rpc_y * rpc_y + 4.0 * rpa_y * rpc_z * rpc_y + 2.0 * rpb_y * rpc_z * rpc_y);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * (fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_y * rpa_y * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y - fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_y * rpc_y - rpa_z * rpa_y * rpa_y * rpc_y - rpa_y * rpa_y * rpb_y * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y + fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y - rpb_y * rpc_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_z[i] += dip_z * fss * b1_vals[i] * (fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * rpb_y + rpa_y * rpa_y * rpb_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpb_y - (3.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_y * rpb_y * rpc_y - rpa_y * rpa_y * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_y * rpc_y * rpc_y + rpb_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-rpc_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * (fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_y * rpa_y * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y - fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_y * rpc_y - rpa_z * rpa_y * rpa_y * rpc_y - rpa_y * rpa_y * rpb_y * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y + fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpa_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y - rpb_y * rpc_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YYZ_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z - fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z - rpa_y * rpa_y * rpb_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z + fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y + rpa_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y - rpa_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y - rpb_z * rpc_z * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_y[i] += dip_y * fss * b1_vals[i] * (fe_0 * rpa_y + 2.0 * rpa_z * rpa_y * rpb_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-fe_0 * rpa_y - fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpc_z - 2.0 * rpa_z * rpb_z * rpc_y - 2.0 * rpa_y * rpb_z * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * (fe_0 * rpc_y + 2.0 * rpa_z * rpc_z * rpc_y + 2.0 * rpa_y * rpc_z * rpc_z + 2.0 * rpb_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z - fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z - rpa_y * rpa_y * rpb_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z + fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y + rpa_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y - rpa_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y - rpb_z * rpc_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * rpb_z + rpa_z * rpa_y * rpa_y + rpa_y * rpa_y * rpb_z);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpb_z - fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpc_y - 2.0 * rpa_y * rpb_z * rpc_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_y * rpa_y * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * (fe_0 * rpc_z + rpa_z * rpc_y * rpc_y + 4.0 * rpa_y * rpc_z * rpc_y + rpb_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z - fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z - rpa_y * rpa_y * rpb_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z + fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y + rpa_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y - rpa_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y - rpb_z * rpc_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YZZ_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpc_z - rpa_z * rpa_z * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_z * rpc_z * rpc_y + rpa_y * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_z * rpa_z * rpa_y * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z - rpa_z * rpa_z * rpa_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_z * rpa_z * rpb_x * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x + rpa_z * rpa_z * rpb_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_z * rpb_x * rpc_z - rpa_z * rpa_z * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_z * rpc_z * rpc_x + rpb_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_z * rpa_z * rpa_y * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z - rpa_z * rpa_z * rpa_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_z * rpa_z * rpb_x * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * 2.0 * rpa_z * rpa_y * rpb_x;

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_z * rpa_y * rpc_x - 2.0 * rpa_z * rpb_x * rpc_y - 2.0 * rpa_y * rpb_x * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * (2.0 * rpa_z * rpc_y * rpc_x + 2.0 * rpa_y * rpc_z * rpc_x + 2.0 * rpb_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_x + rpa_z * rpa_z * rpa_y * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z - rpa_z * rpa_z * rpa_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_z * rpa_z * rpb_x * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YZZ_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z - rpa_z * rpa_z * rpa_y * rpc_y - rpa_z * rpa_z * rpb_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (fe_0 * rpa_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_y + rpa_z * rpa_z * rpc_y * rpc_y + rpa_y * rpb_y * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y - rpa_y * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + (1.0 / 2.0) * fe_0 * rpb_y + rpa_z * rpa_z * rpa_y + rpa_z * rpa_z * rpb_y);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpb_y - fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpc_z - 2.0 * rpa_z * rpb_y * rpc_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-2.0 * rpa_z * rpa_z * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * (fe_0 * rpc_y + 4.0 * rpa_z * rpc_z * rpc_y + rpa_y * rpc_z * rpc_z + rpb_y * rpc_z * rpc_z);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z - rpa_z * rpa_z * rpa_y * rpc_y - rpa_z * rpa_z * rpb_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (fe_0 * rpa_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_y + rpa_z * rpa_z * rpc_y * rpc_y + rpa_y * rpb_y * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y - rpa_y * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_z[i] += dip_z * fss * b1_vals[i] * (fe_0 * rpa_z + 2.0 * rpa_z * rpa_y * rpb_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-fe_0 * rpa_z - fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpc_y - 2.0 * rpa_z * rpb_y * rpc_y - 2.0 * rpa_y * rpb_y * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * (fe_0 * rpc_z + 2.0 * rpa_z * rpc_y * rpc_y + 2.0 * rpa_y * rpc_z * rpc_y + 2.0 * rpb_y * rpc_z * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-2.0 * rpc_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z - rpa_z * rpa_z * rpa_y * rpc_y - rpa_z * rpa_z * rpb_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (fe_0 * rpa_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_y + rpa_z * rpa_z * rpc_y * rpc_y + rpa_y * rpb_y * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y - rpa_y * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_YZZ_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * (fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_z * rpa_z * rpa_y * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-fe_0 * rpa_z * rpa_y - fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_z * rpc_z - rpa_z * rpa_z * rpa_y * rpc_z - rpa_z * rpa_z * rpb_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (2.0 * rpa_z * rpb_z * rpc_z * rpc_y + rpa_z * rpa_z * rpc_z * rpc_y + rpa_y * rpb_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y - rpa_y * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_y[i] += dip_y * fss * b1_vals[i] * (fe_0 * rpa_z + (1.0 / 2.0) * fe_0 * rpb_z + rpa_z * rpa_z * rpb_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpb_z - (3.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpb_z * rpc_z - rpa_z * rpa_z * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 2.0 * rpa_z * rpc_z * rpc_z + rpb_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * (fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_z * rpa_z * rpa_y * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-fe_0 * rpa_z * rpa_y - fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_z * rpc_z - rpa_z * rpa_z * rpa_y * rpc_z - rpa_z * rpa_z * rpb_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (2.0 * rpa_z * rpb_z * rpc_z * rpc_y + rpa_z * rpa_z * rpc_z * rpc_y + rpa_y * rpb_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y - rpa_y * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y + 2.0 * rpa_z * rpa_y * rpb_z + rpa_z * rpa_z * rpa_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpc_y - 4.0 * rpa_z * rpa_y * rpc_z - 2.0 * rpa_z * rpb_z * rpc_y - rpa_z * rpa_z * rpc_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_y * rpb_z * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 4.0 * rpa_z * rpc_z * rpc_y + 3.0 * rpa_y * rpc_z * rpc_z + 2.0 * rpb_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * (fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z + rpa_z * rpa_z * rpa_y * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-fe_0 * rpa_z * rpa_y - fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_z * rpc_z - rpa_z * rpa_z * rpa_y * rpc_z - rpa_z * rpa_z * rpb_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (2.0 * rpa_z * rpb_z * rpc_z * rpc_y + rpa_z * rpa_z * rpc_z * rpc_y + rpa_y * rpb_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y - rpa_y * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_ZZZ_X(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpc_z - 3.0 * rpa_z * rpa_z * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 3.0 * rpa_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_z * rpa_z * rpb_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpb_x * rpc_z - rpa_z * rpa_z * rpa_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpb_x * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 3.0 * rpa_z * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_z * rpa_z * rpb_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpb_x * rpc_z - rpa_z * rpa_z * rpa_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpb_x * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 3.0 * rpa_z * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x + 3.0 * rpa_z * rpa_z * rpb_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x - (3.0 / 2.0) * fe_0 * rpc_x - 6.0 * rpa_z * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 6.0 * rpa_z * rpc_z * rpc_x + 3.0 * rpb_x * rpc_z * rpc_z);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x + rpa_z * rpa_z * rpa_z * rpb_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpb_x * rpc_z - rpa_z * rpa_z * rpa_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpb_x * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 3.0 * rpa_z * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_ZZZ_Y(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_z * rpa_z * rpb_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpb_y * rpc_z - rpa_z * rpa_z * rpa_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpb_y * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 3.0 * rpa_z * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_y[i] += dip_y * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpc_z - 3.0 * rpa_z * rpa_z * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 3.0 * rpa_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * fss * b4_vals[i] * (-rpc_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_z * rpa_z * rpb_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpb_y * rpc_z - rpa_z * rpa_z * rpa_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpb_y * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 3.0 * rpa_z * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y + 3.0 * rpa_z * rpa_z * rpb_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y - (3.0 / 2.0) * fe_0 * rpc_y - 6.0 * rpa_z * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 6.0 * rpa_z * rpc_z * rpc_y + 3.0 * rpb_y * rpc_z * rpc_z);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-3.0 * rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y + rpa_z * rpa_z * rpa_z * rpb_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpb_y * rpc_z - rpa_z * rpa_z * rpa_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpb_y * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 3.0 * rpa_z * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialGeom010FP_ZZZ_Z(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<4> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<5> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<5>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_z * rpc_z - rpa_z * rpa_z * rpa_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_z * rpb_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * 3.0 * rpa_z * rpa_z * rpc_z * rpc_z;

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_z * rpc_z - rpa_z * rpa_z * rpa_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_z * rpb_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * 3.0 * rpa_z * rpa_z * rpc_z * rpc_z;

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z;

        fints_z[i] += dip_z * fss * b1_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z + (3.0 / 2.0) * fe_0 * rpb_z + 3.0 * rpa_z * rpa_z * rpb_z + rpa_z * rpa_z * rpa_z);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpb_z - 6.0 * fe_0 * rpc_z - 6.0 * rpa_z * rpb_z * rpc_z - 6.0 * rpa_z * rpa_z * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * (6.0 * fe_0 * rpc_z + 9.0 * rpa_z * rpc_z * rpc_z + 3.0 * rpb_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * fss * b4_vals[i] * (-4.0 * rpc_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_z * rpc_z - rpa_z * rpa_z * rpa_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 + 3.0 * rpa_z * rpb_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * 3.0 * rpa_z * rpa_z * rpc_z * rpc_z;

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z - rpb_z * rpc_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z;

    }
}

} // geom_npotrec namespace

