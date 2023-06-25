#include "NuclearPotentialRecSG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "BoysFunc.hpp"
#include "T2CDistributor.hpp"

namespace npotrec { // npotrec namespace

auto
compNuclearPotentialSG(      CSubMatrix* matrix,
                       const double charge,
                       const TPoint3D& point,
                       const CGtoBlock&  bra_gto_block,
                       const CGtoBlock&  ket_gto_block,
                       const bool        ang_order,
                       const int64_t     bra_first,
                       const int64_t     bra_last) -> void
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

        simd::loadCoordinates(ket_coords_x,
                              ket_coords_y,
                              ket_coords_z,
                              ket_gto_coords,
                              ket_first,
                              ket_last);

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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes,
                                0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_35, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_17, bra_gto_indexes, ket_gto_indexes,
                                0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyy, 6.0, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes,
                                0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxzz, -24.0, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_35, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes,
                                0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes,
                                0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_17, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyzz, -24.0, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzzz, 8.0, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

        }
    }
}

auto
compPrimitiveNuclearPotentialSG(      TDoubleArray& buffer_xxxx,
                                      TDoubleArray& buffer_xxxy,
                                      TDoubleArray& buffer_xxxz,
                                      TDoubleArray& buffer_xxyy,
                                      TDoubleArray& buffer_xxyz,
                                      TDoubleArray& buffer_xxzz,
                                      TDoubleArray& buffer_xyyy,
                                      TDoubleArray& buffer_xyyz,
                                      TDoubleArray& buffer_xyzz,
                                      TDoubleArray& buffer_xzzz,
                                      TDoubleArray& buffer_yyyy,
                                      TDoubleArray& buffer_yyyz,
                                      TDoubleArray& buffer_yyzz,
                                      TDoubleArray& buffer_yzzz,
                                      TDoubleArray& buffer_zzzz,
                       const double charge,
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

    // set up pointer to integrals buffer(s)

    auto fints_xxxx = buffer_xxxx.data();

    auto fints_xxxy = buffer_xxxy.data();

    auto fints_xxxz = buffer_xxxz.data();

    auto fints_xxyy = buffer_xxyy.data();

    auto fints_xxyz = buffer_xxyz.data();

    auto fints_xxzz = buffer_xxzz.data();

    auto fints_xyyy = buffer_xyyy.data();

    auto fints_xyyz = buffer_xyyz.data();

    auto fints_xyzz = buffer_xyzz.data();

    auto fints_xzzz = buffer_xzzz.data();

    auto fints_yyyy = buffer_yyyy.data();

    auto fints_yyyz = buffer_yyyz.data();

    auto fints_yyzz = buffer_yyzz.data();

    auto fints_yzzz = buffer_yzzz.data();

    auto fints_zzzz = buffer_zzzz.data();

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

    #pragma omp simd aligned(fints_xxxx,\
                             fints_xxxy,\
                             fints_xxxz,\
                             fints_xxyy,\
                             fints_xxyz,\
                             fints_xxzz,\
                             fints_xyyy,\
                             fints_xyyz,\
                             fints_xyzz,\
                             fints_xzzz,\
                             fints_yyyy,\
                             fints_yyyz,\
                             fints_yyzz,\
                             fints_yzzz,\
                             fints_zzzz,\
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

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxxx[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 + rpb_x * rpb_x * rpb_x * rpb_x);

        fints_xxxx[i] += fss * b1_vals[i] * (-6.0 * fe_0 * rpb_x * rpc_x - 3.0 * fe_0 * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 - 4.0 * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xxxx[i] += fss * b2_vals[i] * (6.0 * fe_0 * rpb_x * rpc_x + 3.0 * fe_0 * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 + 6.0 * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpc_x * rpc_x - 4.0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxxx[i] += fss * b4_vals[i] * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxxy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_x + rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xxxy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y - 3.0 * rpb_y * rpb_x * rpb_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_y);

        fints_xxxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpb_y * rpb_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - rpb_y * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxxy[i] += fss * b4_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xxxz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_x + rpb_z * rpb_x * rpb_x * rpb_x);

        fints_xxxz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z - 3.0 * rpb_z * rpb_x * rpb_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpb_z * rpb_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpb_z * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxxz[i] += fss * b4_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_xxyy[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xxyy[i] += fss * b1_vals[i] * (-fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y - fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * fe_0);

        fints_xxyy[i] += fss * b1_vals[i] * (-2.0 * rpb_y * rpb_x * rpb_x * rpc_y - 2.0 * rpb_y * rpb_y * rpb_x * rpc_x);

        fints_xxyy[i] += fss * b2_vals[i] * (fe_0 * rpb_y * rpc_y + fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0);

        fints_xxyy[i] += fss * b2_vals[i] * (4.0 * rpb_y * rpb_x * rpc_y * rpc_x + rpb_y * rpb_y * rpc_x * rpc_x + rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxyy[i] += fss * b4_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xxyz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_y + rpb_z * rpb_y * rpb_x * rpb_x);

        fints_xxyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z - 2.0 * rpb_z * rpb_y * rpb_x * rpc_x - rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxyz[i] += fss * b1_vals[i] * (-rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y + rpb_z * rpb_y * rpc_x * rpc_x + 2.0 * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b2_vals[i] * (2.0 * rpb_y * rpb_x * rpc_z * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpb_z * rpc_y * rpc_x * rpc_x - rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxyz[i] += fss * b4_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xxzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 + rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xxzz[i] += fss * b1_vals[i] * (-fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z - fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * fe_0);

        fints_xxzz[i] += fss * b1_vals[i] * (-2.0 * rpb_z * rpb_x * rpb_x * rpc_z - 2.0 * rpb_z * rpb_z * rpb_x * rpc_x);

        fints_xxzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpc_z + fe_0 * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0);

        fints_xxzz[i] += fss * b2_vals[i] * (4.0 * rpb_z * rpb_x * rpc_z * rpc_x + rpb_z * rpb_z * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x - 2.0 * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxzz[i] += fss * b4_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_xyyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_x + rpb_y * rpb_y * rpb_y * rpb_x);

        fints_xyyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y - 3.0 * rpb_y * rpb_y * rpb_x * rpc_y - rpb_y * rpb_y * rpb_y * rpc_x);

        fints_xyyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x + 3.0 * rpb_y * rpb_x * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xyyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x - 3.0 * rpb_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyyy[i] += fss * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xyyz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_x + rpb_z * rpb_y * rpb_y * rpb_x);

        fints_xyyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z - 2.0 * rpb_z * rpb_y * rpb_x * rpc_y - rpb_z * rpb_y * rpb_y * rpc_x);

        fints_xyyz[i] += fss * b1_vals[i] * (-rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x + 2.0 * rpb_z * rpb_y * rpc_y * rpc_x + rpb_z * rpb_x * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b2_vals[i] * (2.0 * rpb_y * rpb_x * rpc_z * rpc_y + rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_x - rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyyz[i] += fss * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_xyzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_x + rpb_z * rpb_z * rpb_y * rpb_x);

        fints_xyzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y - 2.0 * rpb_z * rpb_y * rpb_x * rpc_z - rpb_z * rpb_z * rpb_y * rpc_x);

        fints_xyzz[i] += fss * b1_vals[i] * (-rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xyzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x + 2.0 * rpb_z * rpb_y * rpc_z * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b2_vals[i] * (rpb_z * rpb_z * rpc_y * rpc_x + rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_y * rpc_x - 2.0 * rpb_z * rpc_z * rpc_y * rpc_x - rpb_y * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyzz[i] += fss * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_xzzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_x + rpb_z * rpb_z * rpb_z * rpb_x);

        fints_xzzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z - 3.0 * rpb_z * rpb_z * rpb_x * rpc_z - rpb_z * rpb_z * rpb_z * rpc_x);

        fints_xzzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x + 3.0 * rpb_z * rpb_x * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x - 3.0 * rpb_z * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xzzz[i] += fss * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_yyyy[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 + rpb_y * rpb_y * rpb_y * rpb_y);

        fints_yyyy[i] += fss * b1_vals[i] * (-6.0 * fe_0 * rpb_y * rpc_y - 3.0 * fe_0 * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 - 4.0 * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yyyy[i] += fss * b2_vals[i] * (6.0 * fe_0 * rpb_y * rpc_y + 3.0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 + 6.0 * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y - 4.0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyyy[i] += fss * b4_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yyyz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_y + rpb_z * rpb_y * rpb_y * rpb_y);

        fints_yyyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z - 3.0 * rpb_z * rpb_y * rpb_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpb_z * rpb_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - rpb_z * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyyz[i] += fss * b4_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_yyzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 + rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yyzz[i] += fss * b1_vals[i] * (-fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z - fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * fe_0);

        fints_yyzz[i] += fss * b1_vals[i] * (-2.0 * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_yyzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpc_z + fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0);

        fints_yyzz[i] += fss * b2_vals[i] * (4.0 * rpb_z * rpb_y * rpc_z * rpc_y + rpb_z * rpb_z * rpc_y * rpc_y + rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y - 2.0 * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyzz[i] += fss * b4_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_yzzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_y + rpb_z * rpb_z * rpb_z * rpb_y);

        fints_yzzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z - 3.0 * rpb_z * rpb_z * rpb_y * rpc_z - rpb_z * rpb_z * rpb_z * rpc_y);

        fints_yzzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y + 3.0 * rpb_z * rpb_y * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yzzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y - 3.0 * rpb_z * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzzz[i] += fss * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_zzzz[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 + rpb_z * rpb_z * rpb_z * rpb_z);

        fints_zzzz[i] += fss * b1_vals[i] * (-6.0 * fe_0 * rpb_z * rpc_z - 3.0 * fe_0 * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * fe_0 - 4.0 * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzzz[i] += fss * b2_vals[i] * (6.0 * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 + 6.0 * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z - 4.0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzzz[i] += fss * b4_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z;

    }
}

} // npotrec namespace

