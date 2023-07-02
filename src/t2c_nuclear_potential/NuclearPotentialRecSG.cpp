#include "NuclearPotentialRecSG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveNuclearPotentialSG.hpp"
#include "T2CDistributor.hpp"

namespace npotrec {  // npotrec namespace

auto
compNuclearPotentialSG(CSubMatrix*      matrix,
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

    alignas(64) TDoubleArray buffer_xxxx;

    alignas(64) TDoubleArray buffer_xxxy;

    alignas(64) TDoubleArray buffer_xxxz;

    alignas(64) TDoubleArray buffer_xxyy;

    alignas(64) TDoubleArray buffer_xxyz;

    alignas(64) TDoubleArray buffer_xxzz;

    alignas(64) TDoubleArray buffer_xyyy;

    alignas(64) TDoubleArray buffer_xyyz;

    alignas(64) TDoubleArray buffer_xyzz;

    alignas(64) TDoubleArray buffer_xzzz;

    alignas(64) TDoubleArray buffer_yyyy;

    alignas(64) TDoubleArray buffer_yyyz;

    alignas(64) TDoubleArray buffer_yyzz;

    alignas(64) TDoubleArray buffer_yzzz;

    alignas(64) TDoubleArray buffer_zzzz;

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

            // compute primitive integrals block

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialSG(buffer_xxxx,
                                                             buffer_xxxy,
                                                             buffer_xxxz,
                                                             buffer_xxyy,
                                                             buffer_xxyz,
                                                             buffer_xxzz,
                                                             buffer_xyyy,
                                                             buffer_xyyz,
                                                             buffer_xyzz,
                                                             buffer_xzzz,
                                                             buffer_yyyy,
                                                             buffer_yyyz,
                                                             buffer_yyzz,
                                                             buffer_yzzz,
                                                             buffer_zzzz,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyy, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzzz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace npotrec
