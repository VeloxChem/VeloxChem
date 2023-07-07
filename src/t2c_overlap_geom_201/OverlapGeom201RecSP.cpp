#include "OverlapGeom201RecSP.hpp"

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom201SP_0_X.hpp"
#include "PrimitiveOverlapGeom201SP_0_Y.hpp"
#include "PrimitiveOverlapGeom201SP_0_Z.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom201SP(CSubMatrix*      matrix_xx_x,
                     CSubMatrix*      matrix_xx_y,
                     CSubMatrix*      matrix_xx_z,
                     CSubMatrix*      matrix_xy_x,
                     CSubMatrix*      matrix_xy_y,
                     CSubMatrix*      matrix_xy_z,
                     CSubMatrix*      matrix_xz_x,
                     CSubMatrix*      matrix_xz_y,
                     CSubMatrix*      matrix_xz_z,
                     CSubMatrix*      matrix_yy_x,
                     CSubMatrix*      matrix_yy_y,
                     CSubMatrix*      matrix_yy_z,
                     CSubMatrix*      matrix_yz_x,
                     CSubMatrix*      matrix_yz_y,
                     CSubMatrix*      matrix_yz_z,
                     CSubMatrix*      matrix_zz_x,
                     CSubMatrix*      matrix_zz_y,
                     CSubMatrix*      matrix_zz_z,
                     const CGtoBlock& bra_gto_block,
                     const CGtoBlock& ket_gto_block,
                     const bool       ang_order,
                     const int64_t    bra_first,
                     const int64_t    bra_last) -> void
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

    alignas(64) TDoubleArray buffer_xx_x;

    alignas(64) TDoubleArray buffer_xx_y;

    alignas(64) TDoubleArray buffer_xx_z;

    alignas(64) TDoubleArray buffer_xy_x;

    alignas(64) TDoubleArray buffer_xy_y;

    alignas(64) TDoubleArray buffer_xy_z;

    alignas(64) TDoubleArray buffer_xz_x;

    alignas(64) TDoubleArray buffer_xz_y;

    alignas(64) TDoubleArray buffer_xz_z;

    alignas(64) TDoubleArray buffer_yy_x;

    alignas(64) TDoubleArray buffer_yy_y;

    alignas(64) TDoubleArray buffer_yy_z;

    alignas(64) TDoubleArray buffer_yz_x;

    alignas(64) TDoubleArray buffer_yz_y;

    alignas(64) TDoubleArray buffer_yz_z;

    alignas(64) TDoubleArray buffer_zz_x;

    alignas(64) TDoubleArray buffer_zz_y;

    alignas(64) TDoubleArray buffer_zz_z;

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

            // compute primitive integrals block (0_X)

            simd::zero(buffer_xx_x);

            simd::zero(buffer_xx_y);

            simd::zero(buffer_xx_z);

            simd::zero(buffer_xy_x);

            simd::zero(buffer_xy_y);

            simd::zero(buffer_xy_z);

            simd::zero(buffer_xz_x);

            simd::zero(buffer_xz_y);

            simd::zero(buffer_xz_z);

            simd::zero(buffer_yy_x);

            simd::zero(buffer_yy_y);

            simd::zero(buffer_yy_z);

            simd::zero(buffer_yz_x);

            simd::zero(buffer_yz_y);

            simd::zero(buffer_yz_z);

            simd::zero(buffer_zz_x);

            simd::zero(buffer_zz_y);

            simd::zero(buffer_zz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom201SP_0_X(buffer_xx_x,
                                                              buffer_xx_y,
                                                              buffer_xx_z,
                                                              buffer_xy_x,
                                                              buffer_xy_y,
                                                              buffer_xy_z,
                                                              buffer_xz_x,
                                                              buffer_xz_y,
                                                              buffer_xz_z,
                                                              buffer_yy_x,
                                                              buffer_yy_y,
                                                              buffer_yy_z,
                                                              buffer_yz_x,
                                                              buffer_yz_y,
                                                              buffer_yz_z,
                                                              buffer_zz_x,
                                                              buffer_zz_y,
                                                              buffer_zz_z,
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

            t2cfunc::distribute(matrix_xx_x, buffer_xx_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_y, buffer_xx_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_z, buffer_xx_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_x, buffer_xy_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_y, buffer_xy_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_z, buffer_xy_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_x, buffer_xz_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_y, buffer_xz_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_z, buffer_xz_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_x, buffer_yy_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_y, buffer_yy_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_z, buffer_yy_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_x, buffer_yz_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_y, buffer_yz_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_z, buffer_yz_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_x, buffer_zz_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_y, buffer_zz_y, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_z, buffer_zz_z, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_Y)

            simd::zero(buffer_xx_x);

            simd::zero(buffer_xx_y);

            simd::zero(buffer_xx_z);

            simd::zero(buffer_xy_x);

            simd::zero(buffer_xy_y);

            simd::zero(buffer_xy_z);

            simd::zero(buffer_xz_x);

            simd::zero(buffer_xz_y);

            simd::zero(buffer_xz_z);

            simd::zero(buffer_yy_x);

            simd::zero(buffer_yy_y);

            simd::zero(buffer_yy_z);

            simd::zero(buffer_yz_x);

            simd::zero(buffer_yz_y);

            simd::zero(buffer_yz_z);

            simd::zero(buffer_zz_x);

            simd::zero(buffer_zz_y);

            simd::zero(buffer_zz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom201SP_0_Y(buffer_xx_x,
                                                              buffer_xx_y,
                                                              buffer_xx_z,
                                                              buffer_xy_x,
                                                              buffer_xy_y,
                                                              buffer_xy_z,
                                                              buffer_xz_x,
                                                              buffer_xz_y,
                                                              buffer_xz_z,
                                                              buffer_yy_x,
                                                              buffer_yy_y,
                                                              buffer_yy_z,
                                                              buffer_yz_x,
                                                              buffer_yz_y,
                                                              buffer_yz_z,
                                                              buffer_zz_x,
                                                              buffer_zz_y,
                                                              buffer_zz_z,
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

            t2cfunc::distribute(matrix_xx_x, buffer_xx_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_y, buffer_xx_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_z, buffer_xx_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_x, buffer_xy_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_y, buffer_xy_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_z, buffer_xy_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_x, buffer_xz_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_y, buffer_xz_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_z, buffer_xz_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_x, buffer_yy_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_y, buffer_yy_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_z, buffer_yy_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_x, buffer_yz_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_y, buffer_yz_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_z, buffer_yz_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_x, buffer_zz_x, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_y, buffer_zz_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_z, buffer_zz_z, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_Z)

            simd::zero(buffer_xx_x);

            simd::zero(buffer_xx_y);

            simd::zero(buffer_xx_z);

            simd::zero(buffer_xy_x);

            simd::zero(buffer_xy_y);

            simd::zero(buffer_xy_z);

            simd::zero(buffer_xz_x);

            simd::zero(buffer_xz_y);

            simd::zero(buffer_xz_z);

            simd::zero(buffer_yy_x);

            simd::zero(buffer_yy_y);

            simd::zero(buffer_yy_z);

            simd::zero(buffer_yz_x);

            simd::zero(buffer_yz_y);

            simd::zero(buffer_yz_z);

            simd::zero(buffer_zz_x);

            simd::zero(buffer_zz_y);

            simd::zero(buffer_zz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom201SP_0_Z(buffer_xx_x,
                                                              buffer_xx_y,
                                                              buffer_xx_z,
                                                              buffer_xy_x,
                                                              buffer_xy_y,
                                                              buffer_xy_z,
                                                              buffer_xz_x,
                                                              buffer_xz_y,
                                                              buffer_xz_z,
                                                              buffer_yy_x,
                                                              buffer_yy_y,
                                                              buffer_yy_z,
                                                              buffer_yz_x,
                                                              buffer_yz_y,
                                                              buffer_yz_z,
                                                              buffer_zz_x,
                                                              buffer_zz_y,
                                                              buffer_zz_z,
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

            t2cfunc::distribute(matrix_xx_x, buffer_xx_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_y, buffer_xx_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xx_z, buffer_xx_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_x, buffer_xy_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_y, buffer_xy_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xy_z, buffer_xy_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_x, buffer_xz_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_y, buffer_xz_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xz_z, buffer_xz_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_x, buffer_yy_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_y, buffer_yy_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yy_z, buffer_yy_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_x, buffer_yz_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_y, buffer_yz_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yz_z, buffer_yz_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_x, buffer_zz_x, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_y, buffer_zz_y, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zz_z, buffer_zz_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace ovlrec
