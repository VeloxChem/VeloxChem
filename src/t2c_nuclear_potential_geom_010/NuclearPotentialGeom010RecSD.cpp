#include "NuclearPotentialGeom010RecSD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveNuclearPotentialGeom010SD_0_XX.hpp"
#include "PrimitiveNuclearPotentialGeom010SD_0_XY.hpp"
#include "PrimitiveNuclearPotentialGeom010SD_0_XZ.hpp"
#include "PrimitiveNuclearPotentialGeom010SD_0_YY.hpp"
#include "PrimitiveNuclearPotentialGeom010SD_0_YZ.hpp"
#include "PrimitiveNuclearPotentialGeom010SD_0_ZZ.hpp"
#include "T2CDistributor.hpp"

namespace npotg010rec {  // npotg010rec namespace

auto
compNuclearPotentialGeom010SD(CSubMatrix*      matrix_x,
                              CSubMatrix*      matrix_y,
                              CSubMatrix*      matrix_z,
                              const TPoint3D&  dipole,
                              const TPoint3D&  point,
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

            // compute primitive integrals block (0_XX)

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

                    npotg010rec::compPrimitiveNuclearPotentialGeom010SD_0_XX(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_XY)

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

                    npotg010rec::compPrimitiveNuclearPotentialGeom010SD_0_XY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_XZ)

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

                    npotg010rec::compPrimitiveNuclearPotentialGeom010SD_0_XZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_YY)

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

                    npotg010rec::compPrimitiveNuclearPotentialGeom010SD_0_YY(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -1.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_YZ)

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

                    npotg010rec::compPrimitiveNuclearPotentialGeom010SD_0_YZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_ZZ)

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

                    npotg010rec::compPrimitiveNuclearPotentialGeom010SD_0_ZZ(buffer_x,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace npotg010rec
