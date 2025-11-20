#include "ElectronRepulsionGeom0010ContrRecXXGS.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxgs(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxgs,
                                            const size_t idx_xxfs,
                                            const size_t idx_geom_10_xxfs,
                                            const size_t idx_geom_10_xxfp,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSFS

            const auto fs_off = idx_xxfs + (i * bcomps + j) * 10;

            auto g_xxx_0 = cbuffer.data(fs_off + 0);

            auto g_yyy_0 = cbuffer.data(fs_off + 6);

            auto g_zzz_0 = cbuffer.data(fs_off + 9);

            /// Set up components of auxilary buffer : SSFS

            const auto fs_geom_10_off = idx_geom_10_xxfs + (i * bcomps + j) * 10;

            auto g_x_0_xxx_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_yyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_yyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_yzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_zzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_y_0_xxx_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 0);

            auto g_y_0_xxy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 1);

            auto g_y_0_xxz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 2);

            auto g_y_0_xyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 3);

            auto g_y_0_xyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 4);

            auto g_y_0_xzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 5);

            auto g_y_0_yyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 6);

            auto g_y_0_yyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 7);

            auto g_y_0_yzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 8);

            auto g_y_0_zzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 9);

            auto g_z_0_xxx_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 0);

            auto g_z_0_xxy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 1);

            auto g_z_0_xxz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 2);

            auto g_z_0_xyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 3);

            auto g_z_0_xyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 4);

            auto g_z_0_xzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 5);

            auto g_z_0_yyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 6);

            auto g_z_0_yyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 7);

            auto g_z_0_yzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 8);

            auto g_z_0_zzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 9);

            /// Set up components of auxilary buffer : SSFP

            const auto fp_geom_10_off = idx_geom_10_xxfp + (i * bcomps + j) * 30;

            auto g_x_0_xxx_x = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxy_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_yyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_yyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_yzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_zzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_zzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_y_0_xxx_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 0);

            auto g_y_0_xxy_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 3);

            auto g_y_0_xxz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 6);

            auto g_y_0_xyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 9);

            auto g_y_0_xyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 12);

            auto g_y_0_xzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 15);

            auto g_y_0_yyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 18);

            auto g_y_0_yyy_y = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 19);

            auto g_y_0_yyy_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 20);

            auto g_y_0_yyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 21);

            auto g_y_0_yyz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 23);

            auto g_y_0_yzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 24);

            auto g_y_0_yzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 26);

            auto g_y_0_zzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 27);

            auto g_y_0_zzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps * bcomps + 29);

            auto g_z_0_xxx_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 0);

            auto g_z_0_xxy_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 3);

            auto g_z_0_xxz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 6);

            auto g_z_0_xyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 9);

            auto g_z_0_xyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 12);

            auto g_z_0_xzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 15);

            auto g_z_0_yyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 18);

            auto g_z_0_yyy_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 19);

            auto g_z_0_yyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 21);

            auto g_z_0_yyz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 22);

            auto g_z_0_yzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 24);

            auto g_z_0_yzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 25);

            auto g_z_0_zzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 27);

            auto g_z_0_zzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 28);

            auto g_z_0_zzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps * bcomps + 29);

            /// set up bra offset for contr_buffer_xxgs

            const auto gs_geom_10_off = idx_geom_10_xxgs + (i * bcomps + j) * 15;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_x_0_xxx_0, g_x_0_xxx_x, g_x_0_xxxx_0, g_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_0[k] = -g_xxx_0[k] - g_x_0_xxx_0[k] * cd_x[k] + g_x_0_xxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_x_0_xxx_0, g_x_0_xxx_y, g_x_0_xxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_0[k] = -g_x_0_xxx_0[k] * cd_y[k] + g_x_0_xxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_x_0_xxx_0, g_x_0_xxx_z, g_x_0_xxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_0[k] = -g_x_0_xxx_0[k] * cd_z[k] + g_x_0_xxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_y, g_x_0_xxy_0, g_x_0_xxy_y, g_x_0_xxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_0[k] = -g_x_0_xxy_0[k] * cd_y[k] + g_x_0_xxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_y, g_x_0_xxyz_0, g_x_0_xxz_0, g_x_0_xxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_0[k] = -g_x_0_xxz_0[k] * cd_y[k] + g_x_0_xxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_x_0_xxz_0, g_x_0_xxz_z, g_x_0_xxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_0[k] = -g_x_0_xxz_0[k] * cd_z[k] + g_x_0_xxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_y, g_x_0_xyy_0, g_x_0_xyy_y, g_x_0_xyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_0[k] = -g_x_0_xyy_0[k] * cd_y[k] + g_x_0_xyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_y, g_x_0_xyyz_0, g_x_0_xyz_0, g_x_0_xyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_0[k] = -g_x_0_xyz_0[k] * cd_y[k] + g_x_0_xyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_y, g_x_0_xyzz_0, g_x_0_xzz_0, g_x_0_xzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_0[k] = -g_x_0_xzz_0[k] * cd_y[k] + g_x_0_xzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_z, g_x_0_xzz_0, g_x_0_xzz_z, g_x_0_xzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_0[k] = -g_x_0_xzz_0[k] * cd_z[k] + g_x_0_xzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_y, g_x_0_yyy_0, g_x_0_yyy_y, g_x_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_0[k] = -g_x_0_yyy_0[k] * cd_y[k] + g_x_0_yyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_yyyz_0, g_x_0_yyz_0, g_x_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_0[k] = -g_x_0_yyz_0[k] * cd_y[k] + g_x_0_yyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_y, g_x_0_yyzz_0, g_x_0_yzz_0, g_x_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_0[k] = -g_x_0_yzz_0[k] * cd_y[k] + g_x_0_yzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_y, g_x_0_yzzz_0, g_x_0_zzz_0, g_x_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_0[k] = -g_x_0_zzz_0[k] * cd_y[k] + g_x_0_zzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_z, g_x_0_zzz_0, g_x_0_zzz_z, g_x_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_0[k] = -g_x_0_zzz_0[k] * cd_z[k] + g_x_0_zzz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_y_0_xxx_0, g_y_0_xxx_x, g_y_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_0[k] = -g_y_0_xxx_0[k] * cd_x[k] + g_y_0_xxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_y_0_xxxy_0, g_y_0_xxy_0, g_y_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_0[k] = -g_y_0_xxy_0[k] * cd_x[k] + g_y_0_xxy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xxxz_0, g_y_0_xxz_0, g_y_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_0[k] = -g_y_0_xxz_0[k] * cd_x[k] + g_y_0_xxz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_y_0_xxyy_0, g_y_0_xyy_0, g_y_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_0[k] = -g_y_0_xyy_0[k] * cd_x[k] + g_y_0_xyy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_y_0_xxyz_0, g_y_0_xyz_0, g_y_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_0[k] = -g_y_0_xyz_0[k] * cd_x[k] + g_y_0_xyz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxzz_0, g_y_0_xzz_0, g_y_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_0[k] = -g_y_0_xzz_0[k] * cd_x[k] + g_y_0_xzz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_x, g_y_0_xyyy_0, g_y_0_yyy_0, g_y_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_0[k] = -g_y_0_yyy_0[k] * cd_x[k] + g_y_0_yyy_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_x, g_y_0_xyyz_0, g_y_0_yyz_0, g_y_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_0[k] = -g_y_0_yyz_0[k] * cd_x[k] + g_y_0_yyz_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_y_0_xyzz_0, g_y_0_yzz_0, g_y_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_0[k] = -g_y_0_yzz_0[k] * cd_x[k] + g_y_0_yzz_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_y_0_xzzz_0, g_y_0_zzz_0, g_y_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_0[k] = -g_y_0_zzz_0[k] * cd_x[k] + g_y_0_zzz_x[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_y, g_y_0_yyy_0, g_y_0_yyy_y, g_y_0_yyyy_0, g_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_0[k] = -g_yyy_0[k] - g_y_0_yyy_0[k] * cd_y[k] + g_y_0_yyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_z, g_y_0_yyy_0, g_y_0_yyy_z, g_y_0_yyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_0[k] = -g_y_0_yyy_0[k] * cd_z[k] + g_y_0_yyy_z[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_z, g_y_0_yyz_0, g_y_0_yyz_z, g_y_0_yyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_0[k] = -g_y_0_yyz_0[k] * cd_z[k] + g_y_0_yyz_z[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_z, g_y_0_yzz_0, g_y_0_yzz_z, g_y_0_yzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_0[k] = -g_y_0_yzz_0[k] * cd_z[k] + g_y_0_yzz_z[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_z, g_y_0_zzz_0, g_y_0_zzz_z, g_y_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_0[k] = -g_y_0_zzz_0[k] * cd_z[k] + g_y_0_zzz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_z_0_xxx_0, g_z_0_xxx_x, g_z_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_0[k] = -g_z_0_xxx_0[k] * cd_x[k] + g_z_0_xxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_z_0_xxxy_0, g_z_0_xxy_0, g_z_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_0[k] = -g_z_0_xxy_0[k] * cd_x[k] + g_z_0_xxy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xxxz_0, g_z_0_xxz_0, g_z_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_0[k] = -g_z_0_xxz_0[k] * cd_x[k] + g_z_0_xxz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_z_0_xxyy_0, g_z_0_xyy_0, g_z_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_0[k] = -g_z_0_xyy_0[k] * cd_x[k] + g_z_0_xyy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_z_0_xxyz_0, g_z_0_xyz_0, g_z_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_0[k] = -g_z_0_xyz_0[k] * cd_x[k] + g_z_0_xyz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxzz_0, g_z_0_xzz_0, g_z_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_0[k] = -g_z_0_xzz_0[k] * cd_x[k] + g_z_0_xzz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_x, g_z_0_xyyy_0, g_z_0_yyy_0, g_z_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_0[k] = -g_z_0_yyy_0[k] * cd_x[k] + g_z_0_yyy_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_x, g_z_0_xyyz_0, g_z_0_yyz_0, g_z_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_0[k] = -g_z_0_yyz_0[k] * cd_x[k] + g_z_0_yyz_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_z_0_xyzz_0, g_z_0_yzz_0, g_z_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_0[k] = -g_z_0_yzz_0[k] * cd_x[k] + g_z_0_yzz_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_z_0_xzzz_0, g_z_0_zzz_0, g_z_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_0[k] = -g_z_0_zzz_0[k] * cd_x[k] + g_z_0_zzz_x[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_y, g_z_0_yyy_0, g_z_0_yyy_y, g_z_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_0[k] = -g_z_0_yyy_0[k] * cd_y[k] + g_z_0_yyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_z_0_yyyz_0, g_z_0_yyz_0, g_z_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_0[k] = -g_z_0_yyz_0[k] * cd_y[k] + g_z_0_yyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_y, g_z_0_yyzz_0, g_z_0_yzz_0, g_z_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_0[k] = -g_z_0_yzz_0[k] * cd_y[k] + g_z_0_yzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_y, g_z_0_yzzz_0, g_z_0_zzz_0, g_z_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_0[k] = -g_z_0_zzz_0[k] * cd_y[k] + g_z_0_zzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_z, g_z_0_zzz_0, g_z_0_zzz_z, g_z_0_zzzz_0, g_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_0[k] = -g_zzz_0[k] - g_z_0_zzz_0[k] * cd_z[k] + g_z_0_zzz_z[k];
            }
        }
    }
}

} // erirec namespace

