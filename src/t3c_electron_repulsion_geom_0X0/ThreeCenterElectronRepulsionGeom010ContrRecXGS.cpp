#include "ThreeCenterElectronRepulsionGeom010ContrRecXGS.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xgs(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xgs,
                                        const size_t idx_xfs,
                                        const size_t idx_geom_10_xfs,
                                        const size_t idx_geom_10_xfp,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SFS

        const auto fs_off = idx_xfs + i * 10;

        auto g_xxx_0 = cbuffer.data(fs_off + 0);

        auto g_xxy_0 = cbuffer.data(fs_off + 1);

        auto g_xxz_0 = cbuffer.data(fs_off + 2);

        auto g_xyy_0 = cbuffer.data(fs_off + 3);

        auto g_xyz_0 = cbuffer.data(fs_off + 4);

        auto g_xzz_0 = cbuffer.data(fs_off + 5);

        auto g_yyy_0 = cbuffer.data(fs_off + 6);

        auto g_yyz_0 = cbuffer.data(fs_off + 7);

        auto g_yzz_0 = cbuffer.data(fs_off + 8);

        auto g_zzz_0 = cbuffer.data(fs_off + 9);

        /// Set up components of auxilary buffer : SFS

        const auto fs_geom_10_off = idx_geom_10_xfs + i * 10;

        auto g_x_0_xxx_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xxy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xxz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 5);

        auto g_x_0_yyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 6);

        auto g_x_0_yyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 7);

        auto g_x_0_yzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 8);

        auto g_x_0_zzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps + 9);

        auto g_y_0_xxx_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 0);

        auto g_y_0_xxy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 1);

        auto g_y_0_xxz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 2);

        auto g_y_0_xyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 3);

        auto g_y_0_xyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 4);

        auto g_y_0_xzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 5);

        auto g_y_0_yyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 6);

        auto g_y_0_yyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 7);

        auto g_y_0_yzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 8);

        auto g_y_0_zzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps + 9);

        auto g_z_0_xxx_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 0);

        auto g_z_0_xxy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 1);

        auto g_z_0_xxz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 2);

        auto g_z_0_xyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 3);

        auto g_z_0_xyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 4);

        auto g_z_0_xzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 5);

        auto g_z_0_yyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 6);

        auto g_z_0_yyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 7);

        auto g_z_0_yzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 8);

        auto g_z_0_zzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps + 9);

        /// Set up components of auxilary buffer : SFP

        const auto fp_geom_10_off = idx_geom_10_xfp + i * 30;

        auto g_x_0_xxx_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xxx_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xxx_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xxy_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xxy_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xxy_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xxz_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xxz_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xxz_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xyy_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xyy_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xyz_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xyz_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 17);

        auto g_x_0_yyy_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 18);

        auto g_x_0_yyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 19);

        auto g_x_0_yyy_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 20);

        auto g_x_0_yyz_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 21);

        auto g_x_0_yyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 22);

        auto g_x_0_yyz_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 23);

        auto g_x_0_yzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 24);

        auto g_x_0_yzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 25);

        auto g_x_0_yzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 26);

        auto g_x_0_zzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps + 27);

        auto g_x_0_zzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps + 28);

        auto g_x_0_zzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps + 29);

        auto g_y_0_xxx_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 0);

        auto g_y_0_xxx_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 1);

        auto g_y_0_xxx_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 2);

        auto g_y_0_xxy_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 3);

        auto g_y_0_xxy_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 4);

        auto g_y_0_xxy_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 5);

        auto g_y_0_xxz_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 6);

        auto g_y_0_xxz_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 7);

        auto g_y_0_xxz_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 8);

        auto g_y_0_xyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 9);

        auto g_y_0_xyy_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 10);

        auto g_y_0_xyy_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 11);

        auto g_y_0_xyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 12);

        auto g_y_0_xyz_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 13);

        auto g_y_0_xyz_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 14);

        auto g_y_0_xzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 15);

        auto g_y_0_xzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 16);

        auto g_y_0_xzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 17);

        auto g_y_0_yyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 18);

        auto g_y_0_yyy_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 19);

        auto g_y_0_yyy_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 20);

        auto g_y_0_yyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 21);

        auto g_y_0_yyz_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 22);

        auto g_y_0_yyz_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 23);

        auto g_y_0_yzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 24);

        auto g_y_0_yzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 25);

        auto g_y_0_yzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 26);

        auto g_y_0_zzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps + 27);

        auto g_y_0_zzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps + 28);

        auto g_y_0_zzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps + 29);

        auto g_z_0_xxx_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 0);

        auto g_z_0_xxx_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 1);

        auto g_z_0_xxx_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 2);

        auto g_z_0_xxy_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 3);

        auto g_z_0_xxy_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 4);

        auto g_z_0_xxy_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 5);

        auto g_z_0_xxz_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 6);

        auto g_z_0_xxz_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 7);

        auto g_z_0_xxz_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 8);

        auto g_z_0_xyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 9);

        auto g_z_0_xyy_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 10);

        auto g_z_0_xyy_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 11);

        auto g_z_0_xyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 12);

        auto g_z_0_xyz_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 13);

        auto g_z_0_xyz_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 14);

        auto g_z_0_xzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 15);

        auto g_z_0_xzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 16);

        auto g_z_0_xzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 17);

        auto g_z_0_yyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 18);

        auto g_z_0_yyy_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 19);

        auto g_z_0_yyy_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 20);

        auto g_z_0_yyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 21);

        auto g_z_0_yyz_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 22);

        auto g_z_0_yyz_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 23);

        auto g_z_0_yzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 24);

        auto g_z_0_yzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 25);

        auto g_z_0_yzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 26);

        auto g_z_0_zzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps + 27);

        auto g_z_0_zzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps + 28);

        auto g_z_0_zzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps + 29);

        /// set up bra offset for contr_buffer_xxgs

        const auto gs_geom_10_off = idx_geom_10_xgs + i * 15;

        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 0);

        #pragma omp simd aligned(cd_x, g_x_0_xxx_0, g_x_0_xxx_x, g_x_0_xxxx_0, g_xxx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxx_0[k] = -g_xxx_0[k] - g_x_0_xxx_0[k] * cd_x[k] + g_x_0_xxx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 1);

        #pragma omp simd aligned(cd_y, g_x_0_xxx_0, g_x_0_xxx_y, g_x_0_xxxy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxy_0[k] = -g_x_0_xxx_0[k] * cd_y[k] + g_x_0_xxx_y[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 2);

        #pragma omp simd aligned(cd_z, g_x_0_xxx_0, g_x_0_xxx_z, g_x_0_xxxz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxz_0[k] = -g_x_0_xxx_0[k] * cd_z[k] + g_x_0_xxx_z[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 3);

        #pragma omp simd aligned(cd_y, g_x_0_xxy_0, g_x_0_xxy_y, g_x_0_xxyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxyy_0[k] = -g_x_0_xxy_0[k] * cd_y[k] + g_x_0_xxy_y[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 4);

        #pragma omp simd aligned(cd_y, g_x_0_xxyz_0, g_x_0_xxz_0, g_x_0_xxz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxyz_0[k] = -g_x_0_xxz_0[k] * cd_y[k] + g_x_0_xxz_y[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_z, g_x_0_xxz_0, g_x_0_xxz_z, g_x_0_xxzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxzz_0[k] = -g_x_0_xxz_0[k] * cd_z[k] + g_x_0_xxz_z[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 6);

        #pragma omp simd aligned(cd_y, g_x_0_xyy_0, g_x_0_xyy_y, g_x_0_xyyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyyy_0[k] = -g_x_0_xyy_0[k] * cd_y[k] + g_x_0_xyy_y[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 7);

        #pragma omp simd aligned(cd_y, g_x_0_xyyz_0, g_x_0_xyz_0, g_x_0_xyz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyyz_0[k] = -g_x_0_xyz_0[k] * cd_y[k] + g_x_0_xyz_y[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 8);

        #pragma omp simd aligned(cd_y, g_x_0_xyzz_0, g_x_0_xzz_0, g_x_0_xzz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyzz_0[k] = -g_x_0_xzz_0[k] * cd_y[k] + g_x_0_xzz_y[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 9);

        #pragma omp simd aligned(cd_z, g_x_0_xzz_0, g_x_0_xzz_z, g_x_0_xzzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzzz_0[k] = -g_x_0_xzz_0[k] * cd_z[k] + g_x_0_xzz_z[k];
        }

        /// Set up 10-11 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 10);

        #pragma omp simd aligned(cd_y, g_x_0_yyy_0, g_x_0_yyy_y, g_x_0_yyyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyyy_0[k] = -g_x_0_yyy_0[k] * cd_y[k] + g_x_0_yyy_y[k];
        }

        /// Set up 11-12 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_x_0_yyyz_0, g_x_0_yyz_0, g_x_0_yyz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyyz_0[k] = -g_x_0_yyz_0[k] * cd_y[k] + g_x_0_yyz_y[k];
        }

        /// Set up 12-13 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 12);

        #pragma omp simd aligned(cd_y, g_x_0_yyzz_0, g_x_0_yzz_0, g_x_0_yzz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyzz_0[k] = -g_x_0_yzz_0[k] * cd_y[k] + g_x_0_yzz_y[k];
        }

        /// Set up 13-14 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 13);

        #pragma omp simd aligned(cd_y, g_x_0_yzzz_0, g_x_0_zzz_0, g_x_0_zzz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzzz_0[k] = -g_x_0_zzz_0[k] * cd_y[k] + g_x_0_zzz_y[k];
        }

        /// Set up 14-15 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps  + 14);

        #pragma omp simd aligned(cd_z, g_x_0_zzz_0, g_x_0_zzz_z, g_x_0_zzzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzzz_0[k] = -g_x_0_zzz_0[k] * cd_z[k] + g_x_0_zzz_z[k];
        }
        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 0);

        #pragma omp simd aligned(cd_x, g_y_0_xxx_0, g_y_0_xxx_x, g_y_0_xxxx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxx_0[k] = -g_y_0_xxx_0[k] * cd_x[k] + g_y_0_xxx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 1);

        #pragma omp simd aligned(cd_x, g_y_0_xxxy_0, g_y_0_xxy_0, g_y_0_xxy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxy_0[k] = -g_y_0_xxy_0[k] * cd_x[k] + g_y_0_xxy_x[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_y_0_xxxz_0, g_y_0_xxz_0, g_y_0_xxz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxz_0[k] = -g_y_0_xxz_0[k] * cd_x[k] + g_y_0_xxz_x[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 3);

        #pragma omp simd aligned(cd_x, g_y_0_xxyy_0, g_y_0_xyy_0, g_y_0_xyy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxyy_0[k] = -g_y_0_xyy_0[k] * cd_x[k] + g_y_0_xyy_x[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 4);

        #pragma omp simd aligned(cd_x, g_y_0_xxyz_0, g_y_0_xyz_0, g_y_0_xyz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxyz_0[k] = -g_y_0_xyz_0[k] * cd_x[k] + g_y_0_xyz_x[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_xxzz_0, g_y_0_xzz_0, g_y_0_xzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxzz_0[k] = -g_y_0_xzz_0[k] * cd_x[k] + g_y_0_xzz_x[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 6);

        #pragma omp simd aligned(cd_x, g_y_0_xyyy_0, g_y_0_yyy_0, g_y_0_yyy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyyy_0[k] = -g_y_0_yyy_0[k] * cd_x[k] + g_y_0_yyy_x[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 7);

        #pragma omp simd aligned(cd_x, g_y_0_xyyz_0, g_y_0_yyz_0, g_y_0_yyz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyyz_0[k] = -g_y_0_yyz_0[k] * cd_x[k] + g_y_0_yyz_x[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 8);

        #pragma omp simd aligned(cd_x, g_y_0_xyzz_0, g_y_0_yzz_0, g_y_0_yzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyzz_0[k] = -g_y_0_yzz_0[k] * cd_x[k] + g_y_0_yzz_x[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_y_0_xzzz_0, g_y_0_zzz_0, g_y_0_zzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzzz_0[k] = -g_y_0_zzz_0[k] * cd_x[k] + g_y_0_zzz_x[k];
        }

        /// Set up 10-11 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 10);

        #pragma omp simd aligned(cd_y, g_y_0_yyy_0, g_y_0_yyy_y, g_y_0_yyyy_0, g_yyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyyy_0[k] = -g_yyy_0[k] - g_y_0_yyy_0[k] * cd_y[k] + g_y_0_yyy_y[k];
        }

        /// Set up 11-12 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 11);

        #pragma omp simd aligned(cd_z, g_y_0_yyy_0, g_y_0_yyy_z, g_y_0_yyyz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyyz_0[k] = -g_y_0_yyy_0[k] * cd_z[k] + g_y_0_yyy_z[k];
        }

        /// Set up 12-13 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 12);

        #pragma omp simd aligned(cd_z, g_y_0_yyz_0, g_y_0_yyz_z, g_y_0_yyzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyzz_0[k] = -g_y_0_yyz_0[k] * cd_z[k] + g_y_0_yyz_z[k];
        }

        /// Set up 13-14 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 13);

        #pragma omp simd aligned(cd_z, g_y_0_yzz_0, g_y_0_yzz_z, g_y_0_yzzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzzz_0[k] = -g_y_0_yzz_0[k] * cd_z[k] + g_y_0_yzz_z[k];
        }

        /// Set up 14-15 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps  + 14);

        #pragma omp simd aligned(cd_z, g_y_0_zzz_0, g_y_0_zzz_z, g_y_0_zzzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzzz_0[k] = -g_y_0_zzz_0[k] * cd_z[k] + g_y_0_zzz_z[k];
        }
        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 0);

        #pragma omp simd aligned(cd_x, g_z_0_xxx_0, g_z_0_xxx_x, g_z_0_xxxx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxx_0[k] = -g_z_0_xxx_0[k] * cd_x[k] + g_z_0_xxx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 1);

        #pragma omp simd aligned(cd_x, g_z_0_xxxy_0, g_z_0_xxy_0, g_z_0_xxy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxy_0[k] = -g_z_0_xxy_0[k] * cd_x[k] + g_z_0_xxy_x[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_z_0_xxxz_0, g_z_0_xxz_0, g_z_0_xxz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxz_0[k] = -g_z_0_xxz_0[k] * cd_x[k] + g_z_0_xxz_x[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 3);

        #pragma omp simd aligned(cd_x, g_z_0_xxyy_0, g_z_0_xyy_0, g_z_0_xyy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxyy_0[k] = -g_z_0_xyy_0[k] * cd_x[k] + g_z_0_xyy_x[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 4);

        #pragma omp simd aligned(cd_x, g_z_0_xxyz_0, g_z_0_xyz_0, g_z_0_xyz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxyz_0[k] = -g_z_0_xyz_0[k] * cd_x[k] + g_z_0_xyz_x[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_xxzz_0, g_z_0_xzz_0, g_z_0_xzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxzz_0[k] = -g_z_0_xzz_0[k] * cd_x[k] + g_z_0_xzz_x[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 6);

        #pragma omp simd aligned(cd_x, g_z_0_xyyy_0, g_z_0_yyy_0, g_z_0_yyy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyyy_0[k] = -g_z_0_yyy_0[k] * cd_x[k] + g_z_0_yyy_x[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 7);

        #pragma omp simd aligned(cd_x, g_z_0_xyyz_0, g_z_0_yyz_0, g_z_0_yyz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyyz_0[k] = -g_z_0_yyz_0[k] * cd_x[k] + g_z_0_yyz_x[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 8);

        #pragma omp simd aligned(cd_x, g_z_0_xyzz_0, g_z_0_yzz_0, g_z_0_yzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyzz_0[k] = -g_z_0_yzz_0[k] * cd_x[k] + g_z_0_yzz_x[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_z_0_xzzz_0, g_z_0_zzz_0, g_z_0_zzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzzz_0[k] = -g_z_0_zzz_0[k] * cd_x[k] + g_z_0_zzz_x[k];
        }

        /// Set up 10-11 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 10);

        #pragma omp simd aligned(cd_y, g_z_0_yyy_0, g_z_0_yyy_y, g_z_0_yyyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyyy_0[k] = -g_z_0_yyy_0[k] * cd_y[k] + g_z_0_yyy_y[k];
        }

        /// Set up 11-12 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_z_0_yyyz_0, g_z_0_yyz_0, g_z_0_yyz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyyz_0[k] = -g_z_0_yyz_0[k] * cd_y[k] + g_z_0_yyz_y[k];
        }

        /// Set up 12-13 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 12);

        #pragma omp simd aligned(cd_y, g_z_0_yyzz_0, g_z_0_yzz_0, g_z_0_yzz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyzz_0[k] = -g_z_0_yzz_0[k] * cd_y[k] + g_z_0_yzz_y[k];
        }

        /// Set up 13-14 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 13);

        #pragma omp simd aligned(cd_y, g_z_0_yzzz_0, g_z_0_zzz_0, g_z_0_zzz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzzz_0[k] = -g_z_0_zzz_0[k] * cd_y[k] + g_z_0_zzz_y[k];
        }

        /// Set up 14-15 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps  + 14);

        #pragma omp simd aligned(cd_z, g_z_0_zzz_0, g_z_0_zzz_z, g_z_0_zzzz_0, g_zzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzzz_0[k] = -g_zzz_0[k] - g_z_0_zzz_0[k] * cd_z[k] + g_z_0_zzz_z[k];
        }
    }
}

} // t3ceri namespace

