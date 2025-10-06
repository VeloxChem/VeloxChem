#include "ElectronRepulsionGeom0100ContrRecGSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_gsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_gsxx,
                                            const size_t idx_fsxx,
                                            const size_t idx_geom_01_fsxx,
                                            const size_t idx_geom_01_fpxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : FSSS

            const auto fs_off = idx_fsxx + i * dcomps + j;

            auto g_xxx_0 = cbuffer.data(fs_off + 0 * ccomps * dcomps);

            auto g_xxy_0 = cbuffer.data(fs_off + 1 * ccomps * dcomps);

            auto g_xxz_0 = cbuffer.data(fs_off + 2 * ccomps * dcomps);

            auto g_xyy_0 = cbuffer.data(fs_off + 3 * ccomps * dcomps);

            auto g_xyz_0 = cbuffer.data(fs_off + 4 * ccomps * dcomps);

            auto g_xzz_0 = cbuffer.data(fs_off + 5 * ccomps * dcomps);

            auto g_yyy_0 = cbuffer.data(fs_off + 6 * ccomps * dcomps);

            auto g_yyz_0 = cbuffer.data(fs_off + 7 * ccomps * dcomps);

            auto g_yzz_0 = cbuffer.data(fs_off + 8 * ccomps * dcomps);

            auto g_zzz_0 = cbuffer.data(fs_off + 9 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FSSS

            const auto fs_geom_01_off = idx_geom_01_fsxx + i * dcomps + j;

            auto g_0_x_xxx_0 = cbuffer.data(fs_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxy_0 = cbuffer.data(fs_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxz_0 = cbuffer.data(fs_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xyy_0 = cbuffer.data(fs_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xyz_0 = cbuffer.data(fs_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xzz_0 = cbuffer.data(fs_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_yyy_0 = cbuffer.data(fs_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_yyz_0 = cbuffer.data(fs_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_yzz_0 = cbuffer.data(fs_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_zzz_0 = cbuffer.data(fs_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_xxx_0 = cbuffer.data(fs_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_xxy_0 = cbuffer.data(fs_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_y_xxz_0 = cbuffer.data(fs_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_y_xyy_0 = cbuffer.data(fs_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_y_xyz_0 = cbuffer.data(fs_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_y_xzz_0 = cbuffer.data(fs_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_yyy_0 = cbuffer.data(fs_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_yyz_0 = cbuffer.data(fs_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_yzz_0 = cbuffer.data(fs_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_zzz_0 = cbuffer.data(fs_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_z_xxx_0 = cbuffer.data(fs_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_z_xxy_0 = cbuffer.data(fs_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_z_xxz_0 = cbuffer.data(fs_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_z_xyy_0 = cbuffer.data(fs_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_z_xyz_0 = cbuffer.data(fs_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_z_xzz_0 = cbuffer.data(fs_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_z_yyy_0 = cbuffer.data(fs_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_z_yyz_0 = cbuffer.data(fs_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_z_yzz_0 = cbuffer.data(fs_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_z_zzz_0 = cbuffer.data(fs_geom_01_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FPSS

            const auto fp_geom_01_off = idx_geom_01_fpxx + i * dcomps + j;

            auto g_0_x_xxx_x = cbuffer.data(fp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_y = cbuffer.data(fp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_z = cbuffer.data(fp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxy_x = cbuffer.data(fp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxy_y = cbuffer.data(fp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxy_z = cbuffer.data(fp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxz_x = cbuffer.data(fp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxz_y = cbuffer.data(fp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxz_z = cbuffer.data(fp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xyy_x = cbuffer.data(fp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xyy_y = cbuffer.data(fp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xyy_z = cbuffer.data(fp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xyz_x = cbuffer.data(fp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xyz_y = cbuffer.data(fp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xyz_z = cbuffer.data(fp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xzz_x = cbuffer.data(fp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xzz_y = cbuffer.data(fp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xzz_z = cbuffer.data(fp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_yyy_x = cbuffer.data(fp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_yyy_y = cbuffer.data(fp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_yyy_z = cbuffer.data(fp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_yyz_x = cbuffer.data(fp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_yyz_y = cbuffer.data(fp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_yyz_z = cbuffer.data(fp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_yzz_x = cbuffer.data(fp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_yzz_y = cbuffer.data(fp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_yzz_z = cbuffer.data(fp_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_zzz_x = cbuffer.data(fp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_zzz_y = cbuffer.data(fp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_zzz_z = cbuffer.data(fp_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_xxx_x = cbuffer.data(fp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_xxx_y = cbuffer.data(fp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_xxx_z = cbuffer.data(fp_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_xxy_x = cbuffer.data(fp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_xxy_y = cbuffer.data(fp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_xxy_z = cbuffer.data(fp_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_y_xxz_x = cbuffer.data(fp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_xxz_y = cbuffer.data(fp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_xxz_z = cbuffer.data(fp_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_xyy_x = cbuffer.data(fp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_xyy_y = cbuffer.data(fp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_xyy_z = cbuffer.data(fp_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_y_xyz_x = cbuffer.data(fp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_y_xyz_y = cbuffer.data(fp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_y_xyz_z = cbuffer.data(fp_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_xzz_x = cbuffer.data(fp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_xzz_y = cbuffer.data(fp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_xzz_z = cbuffer.data(fp_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_yyy_x = cbuffer.data(fp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_yyy_y = cbuffer.data(fp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_yyy_z = cbuffer.data(fp_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_yyz_x = cbuffer.data(fp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_yyz_y = cbuffer.data(fp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_yyz_z = cbuffer.data(fp_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_yzz_x = cbuffer.data(fp_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_yzz_y = cbuffer.data(fp_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_yzz_z = cbuffer.data(fp_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_y_zzz_x = cbuffer.data(fp_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_zzz_y = cbuffer.data(fp_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_zzz_z = cbuffer.data(fp_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_z_xxx_x = cbuffer.data(fp_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_z_xxx_y = cbuffer.data(fp_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_z_xxx_z = cbuffer.data(fp_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_z_xxy_x = cbuffer.data(fp_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_z_xxy_y = cbuffer.data(fp_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_z_xxy_z = cbuffer.data(fp_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_z_xxz_x = cbuffer.data(fp_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_z_xxz_y = cbuffer.data(fp_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_z_xxz_z = cbuffer.data(fp_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_z_xyy_x = cbuffer.data(fp_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_z_xyy_y = cbuffer.data(fp_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_z_xyy_z = cbuffer.data(fp_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_z_xyz_x = cbuffer.data(fp_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_z_xyz_y = cbuffer.data(fp_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_z_xyz_z = cbuffer.data(fp_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_z_xzz_x = cbuffer.data(fp_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_z_xzz_y = cbuffer.data(fp_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_z_xzz_z = cbuffer.data(fp_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_z_yyy_x = cbuffer.data(fp_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_z_yyy_y = cbuffer.data(fp_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_z_yyy_z = cbuffer.data(fp_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_z_yyz_x = cbuffer.data(fp_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_z_yyz_y = cbuffer.data(fp_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_z_yyz_z = cbuffer.data(fp_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_z_yzz_x = cbuffer.data(fp_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_z_yzz_y = cbuffer.data(fp_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_z_yzz_z = cbuffer.data(fp_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_z_zzz_x = cbuffer.data(fp_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_z_zzz_y = cbuffer.data(fp_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_z_zzz_z = cbuffer.data(fp_geom_01_off + 89 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gsxx

            const auto gs_geom_01_off = idx_geom_01_gsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxx_0 = cbuffer.data(gs_geom_01_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_0, g_0_x_xxx_x, g_0_x_xxxx_0, g_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxx_0[k] = g_xxx_0[k] - g_0_x_xxx_0[k] * ab_x + g_0_x_xxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxy_0 = cbuffer.data(gs_geom_01_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_0, g_0_x_xxx_y, g_0_x_xxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxy_0[k] = -g_0_x_xxx_0[k] * ab_y + g_0_x_xxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxz_0 = cbuffer.data(gs_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_0, g_0_x_xxx_z, g_0_x_xxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxz_0[k] = -g_0_x_xxx_0[k] * ab_z + g_0_x_xxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyy_0 = cbuffer.data(gs_geom_01_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxy_0, g_0_x_xxy_y, g_0_x_xxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyy_0[k] = -g_0_x_xxy_0[k] * ab_y + g_0_x_xxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyz_0 = cbuffer.data(gs_geom_01_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyz_0, g_0_x_xxz_0, g_0_x_xxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyz_0[k] = -g_0_x_xxz_0[k] * ab_y + g_0_x_xxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzz_0 = cbuffer.data(gs_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxz_0, g_0_x_xxz_z, g_0_x_xxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzz_0[k] = -g_0_x_xxz_0[k] * ab_z + g_0_x_xxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyy_0 = cbuffer.data(gs_geom_01_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyy_0, g_0_x_xyy_y, g_0_x_xyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyy_0[k] = -g_0_x_xyy_0[k] * ab_y + g_0_x_xyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyz_0 = cbuffer.data(gs_geom_01_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyz_0, g_0_x_xyz_0, g_0_x_xyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyz_0[k] = -g_0_x_xyz_0[k] * ab_y + g_0_x_xyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzz_0 = cbuffer.data(gs_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzz_0, g_0_x_xzz_0, g_0_x_xzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzz_0[k] = -g_0_x_xzz_0[k] * ab_y + g_0_x_xzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzz_0 = cbuffer.data(gs_geom_01_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzz_0, g_0_x_xzz_z, g_0_x_xzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzz_0[k] = -g_0_x_xzz_0[k] * ab_z + g_0_x_xzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyy_0 = cbuffer.data(gs_geom_01_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyy_0, g_0_x_yyy_y, g_0_x_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyy_0[k] = -g_0_x_yyy_0[k] * ab_y + g_0_x_yyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyz_0 = cbuffer.data(gs_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyz_0, g_0_x_yyz_0, g_0_x_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyz_0[k] = -g_0_x_yyz_0[k] * ab_y + g_0_x_yyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzz_0 = cbuffer.data(gs_geom_01_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzz_0, g_0_x_yzz_0, g_0_x_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzz_0[k] = -g_0_x_yzz_0[k] * ab_y + g_0_x_yzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzz_0 = cbuffer.data(gs_geom_01_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzz_0, g_0_x_zzz_0, g_0_x_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzz_0[k] = -g_0_x_zzz_0[k] * ab_y + g_0_x_zzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzz_0 = cbuffer.data(gs_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzz_0, g_0_x_zzz_z, g_0_x_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzz_0[k] = -g_0_x_zzz_0[k] * ab_z + g_0_x_zzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxx_0 = cbuffer.data(gs_geom_01_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxx_0, g_0_y_xxx_x, g_0_y_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxx_0[k] = -g_0_y_xxx_0[k] * ab_x + g_0_y_xxx_x[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxy_0 = cbuffer.data(gs_geom_01_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxy_0, g_0_y_xxy_0, g_0_y_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxy_0[k] = -g_0_y_xxy_0[k] * ab_x + g_0_y_xxy_x[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxz_0 = cbuffer.data(gs_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxz_0, g_0_y_xxz_0, g_0_y_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxz_0[k] = -g_0_y_xxz_0[k] * ab_x + g_0_y_xxz_x[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyy_0 = cbuffer.data(gs_geom_01_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyy_0, g_0_y_xyy_0, g_0_y_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyy_0[k] = -g_0_y_xyy_0[k] * ab_x + g_0_y_xyy_x[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyz_0 = cbuffer.data(gs_geom_01_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyz_0, g_0_y_xyz_0, g_0_y_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyz_0[k] = -g_0_y_xyz_0[k] * ab_x + g_0_y_xyz_x[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzz_0 = cbuffer.data(gs_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzz_0, g_0_y_xzz_0, g_0_y_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzz_0[k] = -g_0_y_xzz_0[k] * ab_x + g_0_y_xzz_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyy_0 = cbuffer.data(gs_geom_01_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyy_0, g_0_y_yyy_0, g_0_y_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyy_0[k] = -g_0_y_yyy_0[k] * ab_x + g_0_y_yyy_x[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyz_0 = cbuffer.data(gs_geom_01_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyz_0, g_0_y_yyz_0, g_0_y_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyz_0[k] = -g_0_y_yyz_0[k] * ab_x + g_0_y_yyz_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzz_0 = cbuffer.data(gs_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzz_0, g_0_y_yzz_0, g_0_y_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzz_0[k] = -g_0_y_yzz_0[k] * ab_x + g_0_y_yzz_x[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzz_0 = cbuffer.data(gs_geom_01_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzz_0, g_0_y_zzz_0, g_0_y_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzz_0[k] = -g_0_y_zzz_0[k] * ab_x + g_0_y_zzz_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyy_0 = cbuffer.data(gs_geom_01_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_0, g_0_y_yyy_y, g_0_y_yyyy_0, g_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyy_0[k] = g_yyy_0[k] - g_0_y_yyy_0[k] * ab_y + g_0_y_yyy_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyz_0 = cbuffer.data(gs_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_0, g_0_y_yyy_z, g_0_y_yyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyz_0[k] = -g_0_y_yyy_0[k] * ab_z + g_0_y_yyy_z[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzz_0 = cbuffer.data(gs_geom_01_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyz_0, g_0_y_yyz_z, g_0_y_yyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzz_0[k] = -g_0_y_yyz_0[k] * ab_z + g_0_y_yyz_z[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzz_0 = cbuffer.data(gs_geom_01_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzz_0, g_0_y_yzz_z, g_0_y_yzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzz_0[k] = -g_0_y_yzz_0[k] * ab_z + g_0_y_yzz_z[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzz_0 = cbuffer.data(gs_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzz_0, g_0_y_zzz_z, g_0_y_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzz_0[k] = -g_0_y_zzz_0[k] * ab_z + g_0_y_zzz_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxx_0 = cbuffer.data(gs_geom_01_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxx_0, g_0_z_xxx_x, g_0_z_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxx_0[k] = -g_0_z_xxx_0[k] * ab_x + g_0_z_xxx_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxy_0 = cbuffer.data(gs_geom_01_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxy_0, g_0_z_xxy_0, g_0_z_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxy_0[k] = -g_0_z_xxy_0[k] * ab_x + g_0_z_xxy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxz_0 = cbuffer.data(gs_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxz_0, g_0_z_xxz_0, g_0_z_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxz_0[k] = -g_0_z_xxz_0[k] * ab_x + g_0_z_xxz_x[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyy_0 = cbuffer.data(gs_geom_01_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyy_0, g_0_z_xyy_0, g_0_z_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyy_0[k] = -g_0_z_xyy_0[k] * ab_x + g_0_z_xyy_x[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyz_0 = cbuffer.data(gs_geom_01_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyz_0, g_0_z_xyz_0, g_0_z_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyz_0[k] = -g_0_z_xyz_0[k] * ab_x + g_0_z_xyz_x[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzz_0 = cbuffer.data(gs_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzz_0, g_0_z_xzz_0, g_0_z_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzz_0[k] = -g_0_z_xzz_0[k] * ab_x + g_0_z_xzz_x[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyy_0 = cbuffer.data(gs_geom_01_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyy_0, g_0_z_yyy_0, g_0_z_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyy_0[k] = -g_0_z_yyy_0[k] * ab_x + g_0_z_yyy_x[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyz_0 = cbuffer.data(gs_geom_01_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyz_0, g_0_z_yyz_0, g_0_z_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyz_0[k] = -g_0_z_yyz_0[k] * ab_x + g_0_z_yyz_x[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzz_0 = cbuffer.data(gs_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzz_0, g_0_z_yzz_0, g_0_z_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzz_0[k] = -g_0_z_yzz_0[k] * ab_x + g_0_z_yzz_x[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzz_0 = cbuffer.data(gs_geom_01_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzz_0, g_0_z_zzz_0, g_0_z_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzz_0[k] = -g_0_z_zzz_0[k] * ab_x + g_0_z_zzz_x[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyy_0 = cbuffer.data(gs_geom_01_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyy_0, g_0_z_yyy_y, g_0_z_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyy_0[k] = -g_0_z_yyy_0[k] * ab_y + g_0_z_yyy_y[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyz_0 = cbuffer.data(gs_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyz_0, g_0_z_yyz_0, g_0_z_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyz_0[k] = -g_0_z_yyz_0[k] * ab_y + g_0_z_yyz_y[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzz_0 = cbuffer.data(gs_geom_01_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzz_0, g_0_z_yzz_0, g_0_z_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzz_0[k] = -g_0_z_yzz_0[k] * ab_y + g_0_z_yzz_y[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzz_0 = cbuffer.data(gs_geom_01_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzz_0, g_0_z_zzz_0, g_0_z_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzz_0[k] = -g_0_z_zzz_0[k] * ab_y + g_0_z_zzz_y[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzz_0 = cbuffer.data(gs_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_0, g_0_z_zzz_z, g_0_z_zzzz_0, g_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzz_0[k] = g_zzz_0[k] - g_0_z_zzz_0[k] * ab_z + g_0_z_zzz_z[k];
            }
        }
    }
}

} // erirec namespace

