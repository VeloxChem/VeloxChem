#include "ElectronRepulsionContrRecDFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_dfxx(CSimdArray<double>& contr_buffer_dfxx,
                                     const CSimdArray<double>& contr_buffer_pfxx,
                                     const CSimdArray<double>& contr_buffer_pgxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto ndims = contr_buffer_dfxx.number_of_columns();

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_pfxx

            const auto pf_off = i * dcomps + j;

            auto g_x_xxx = contr_buffer_pfxx[pf_off + 0 * ccomps * dcomps];

            auto g_x_xxy = contr_buffer_pfxx[pf_off + 1 * ccomps * dcomps];

            auto g_x_xxz = contr_buffer_pfxx[pf_off + 2 * ccomps * dcomps];

            auto g_x_xyy = contr_buffer_pfxx[pf_off + 3 * ccomps * dcomps];

            auto g_x_xyz = contr_buffer_pfxx[pf_off + 4 * ccomps * dcomps];

            auto g_x_xzz = contr_buffer_pfxx[pf_off + 5 * ccomps * dcomps];

            auto g_x_yyy = contr_buffer_pfxx[pf_off + 6 * ccomps * dcomps];

            auto g_x_yyz = contr_buffer_pfxx[pf_off + 7 * ccomps * dcomps];

            auto g_x_yzz = contr_buffer_pfxx[pf_off + 8 * ccomps * dcomps];

            auto g_x_zzz = contr_buffer_pfxx[pf_off + 9 * ccomps * dcomps];

            auto g_y_xxx = contr_buffer_pfxx[pf_off + 10 * ccomps * dcomps];

            auto g_y_xxy = contr_buffer_pfxx[pf_off + 11 * ccomps * dcomps];

            auto g_y_xxz = contr_buffer_pfxx[pf_off + 12 * ccomps * dcomps];

            auto g_y_xyy = contr_buffer_pfxx[pf_off + 13 * ccomps * dcomps];

            auto g_y_xyz = contr_buffer_pfxx[pf_off + 14 * ccomps * dcomps];

            auto g_y_xzz = contr_buffer_pfxx[pf_off + 15 * ccomps * dcomps];

            auto g_y_yyy = contr_buffer_pfxx[pf_off + 16 * ccomps * dcomps];

            auto g_y_yyz = contr_buffer_pfxx[pf_off + 17 * ccomps * dcomps];

            auto g_y_yzz = contr_buffer_pfxx[pf_off + 18 * ccomps * dcomps];

            auto g_y_zzz = contr_buffer_pfxx[pf_off + 19 * ccomps * dcomps];

            auto g_z_xxx = contr_buffer_pfxx[pf_off + 20 * ccomps * dcomps];

            auto g_z_xxy = contr_buffer_pfxx[pf_off + 21 * ccomps * dcomps];

            auto g_z_xxz = contr_buffer_pfxx[pf_off + 22 * ccomps * dcomps];

            auto g_z_xyy = contr_buffer_pfxx[pf_off + 23 * ccomps * dcomps];

            auto g_z_xyz = contr_buffer_pfxx[pf_off + 24 * ccomps * dcomps];

            auto g_z_xzz = contr_buffer_pfxx[pf_off + 25 * ccomps * dcomps];

            auto g_z_yyy = contr_buffer_pfxx[pf_off + 26 * ccomps * dcomps];

            auto g_z_yyz = contr_buffer_pfxx[pf_off + 27 * ccomps * dcomps];

            auto g_z_yzz = contr_buffer_pfxx[pf_off + 28 * ccomps * dcomps];

            auto g_z_zzz = contr_buffer_pfxx[pf_off + 29 * ccomps * dcomps];

            /// Set up components of auxilary buffer : contr_buffer_pgxx

            const auto pg_off = i * dcomps + j;

            auto g_x_xxxx = contr_buffer_pgxx[pg_off + 0 * ccomps * dcomps];

            auto g_x_xxxy = contr_buffer_pgxx[pg_off + 1 * ccomps * dcomps];

            auto g_x_xxxz = contr_buffer_pgxx[pg_off + 2 * ccomps * dcomps];

            auto g_x_xxyy = contr_buffer_pgxx[pg_off + 3 * ccomps * dcomps];

            auto g_x_xxyz = contr_buffer_pgxx[pg_off + 4 * ccomps * dcomps];

            auto g_x_xxzz = contr_buffer_pgxx[pg_off + 5 * ccomps * dcomps];

            auto g_x_xyyy = contr_buffer_pgxx[pg_off + 6 * ccomps * dcomps];

            auto g_x_xyyz = contr_buffer_pgxx[pg_off + 7 * ccomps * dcomps];

            auto g_x_xyzz = contr_buffer_pgxx[pg_off + 8 * ccomps * dcomps];

            auto g_x_xzzz = contr_buffer_pgxx[pg_off + 9 * ccomps * dcomps];

            auto g_y_xxxx = contr_buffer_pgxx[pg_off + 15 * ccomps * dcomps];

            auto g_y_xxxy = contr_buffer_pgxx[pg_off + 16 * ccomps * dcomps];

            auto g_y_xxxz = contr_buffer_pgxx[pg_off + 17 * ccomps * dcomps];

            auto g_y_xxyy = contr_buffer_pgxx[pg_off + 18 * ccomps * dcomps];

            auto g_y_xxyz = contr_buffer_pgxx[pg_off + 19 * ccomps * dcomps];

            auto g_y_xxzz = contr_buffer_pgxx[pg_off + 20 * ccomps * dcomps];

            auto g_y_xyyy = contr_buffer_pgxx[pg_off + 21 * ccomps * dcomps];

            auto g_y_xyyz = contr_buffer_pgxx[pg_off + 22 * ccomps * dcomps];

            auto g_y_xyzz = contr_buffer_pgxx[pg_off + 23 * ccomps * dcomps];

            auto g_y_xzzz = contr_buffer_pgxx[pg_off + 24 * ccomps * dcomps];

            auto g_y_yyyy = contr_buffer_pgxx[pg_off + 25 * ccomps * dcomps];

            auto g_y_yyyz = contr_buffer_pgxx[pg_off + 26 * ccomps * dcomps];

            auto g_y_yyzz = contr_buffer_pgxx[pg_off + 27 * ccomps * dcomps];

            auto g_y_yzzz = contr_buffer_pgxx[pg_off + 28 * ccomps * dcomps];

            auto g_z_xxxx = contr_buffer_pgxx[pg_off + 30 * ccomps * dcomps];

            auto g_z_xxxy = contr_buffer_pgxx[pg_off + 31 * ccomps * dcomps];

            auto g_z_xxxz = contr_buffer_pgxx[pg_off + 32 * ccomps * dcomps];

            auto g_z_xxyy = contr_buffer_pgxx[pg_off + 33 * ccomps * dcomps];

            auto g_z_xxyz = contr_buffer_pgxx[pg_off + 34 * ccomps * dcomps];

            auto g_z_xxzz = contr_buffer_pgxx[pg_off + 35 * ccomps * dcomps];

            auto g_z_xyyy = contr_buffer_pgxx[pg_off + 36 * ccomps * dcomps];

            auto g_z_xyyz = contr_buffer_pgxx[pg_off + 37 * ccomps * dcomps];

            auto g_z_xyzz = contr_buffer_pgxx[pg_off + 38 * ccomps * dcomps];

            auto g_z_xzzz = contr_buffer_pgxx[pg_off + 39 * ccomps * dcomps];

            auto g_z_yyyy = contr_buffer_pgxx[pg_off + 40 * ccomps * dcomps];

            auto g_z_yyyz = contr_buffer_pgxx[pg_off + 41 * ccomps * dcomps];

            auto g_z_yyzz = contr_buffer_pgxx[pg_off + 42 * ccomps * dcomps];

            auto g_z_yzzz = contr_buffer_pgxx[pg_off + 43 * ccomps * dcomps];

            auto g_z_zzzz = contr_buffer_pgxx[pg_off + 44 * ccomps * dcomps];

            /// set up bra offset for contr_buffer_dfxx

            const auto df_off = i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : contr_buffer_dfxx

            auto g_xx_xxx = contr_buffer_dfxx[df_off + 0 * ccomps * dcomps];

            auto g_xx_xxy = contr_buffer_dfxx[df_off + 1 * ccomps * dcomps];

            auto g_xx_xxz = contr_buffer_dfxx[df_off + 2 * ccomps * dcomps];

            auto g_xx_xyy = contr_buffer_dfxx[df_off + 3 * ccomps * dcomps];

            auto g_xx_xyz = contr_buffer_dfxx[df_off + 4 * ccomps * dcomps];

            auto g_xx_xzz = contr_buffer_dfxx[df_off + 5 * ccomps * dcomps];

            auto g_xx_yyy = contr_buffer_dfxx[df_off + 6 * ccomps * dcomps];

            auto g_xx_yyz = contr_buffer_dfxx[df_off + 7 * ccomps * dcomps];

            auto g_xx_yzz = contr_buffer_dfxx[df_off + 8 * ccomps * dcomps];

            auto g_xx_zzz = contr_buffer_dfxx[df_off + 9 * ccomps * dcomps];

            #pragma omp simd aligned(g_x_xxx, g_x_xxxx, g_x_xxxy, g_x_xxxz, g_x_xxy, g_x_xxyy, g_x_xxyz, g_x_xxz, g_x_xxzz, g_x_xyy, g_x_xyyy, g_x_xyyz, g_x_xyz, g_x_xyzz, g_x_xzz, g_x_xzzz, g_x_yyy, g_x_yyz, g_x_yzz, g_x_zzz, g_xx_xxx, g_xx_xxy, g_xx_xxz, g_xx_xyy, g_xx_xyz, g_xx_xzz, g_xx_yyy, g_xx_yyz, g_xx_yzz, g_xx_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xx_xxx[k] = -g_x_xxx[k] * ab_x + g_x_xxxx[k];

                g_xx_xxy[k] = -g_x_xxy[k] * ab_x + g_x_xxxy[k];

                g_xx_xxz[k] = -g_x_xxz[k] * ab_x + g_x_xxxz[k];

                g_xx_xyy[k] = -g_x_xyy[k] * ab_x + g_x_xxyy[k];

                g_xx_xyz[k] = -g_x_xyz[k] * ab_x + g_x_xxyz[k];

                g_xx_xzz[k] = -g_x_xzz[k] * ab_x + g_x_xxzz[k];

                g_xx_yyy[k] = -g_x_yyy[k] * ab_x + g_x_xyyy[k];

                g_xx_yyz[k] = -g_x_yyz[k] * ab_x + g_x_xyyz[k];

                g_xx_yzz[k] = -g_x_yzz[k] * ab_x + g_x_xyzz[k];

                g_xx_zzz[k] = -g_x_zzz[k] * ab_x + g_x_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : contr_buffer_dfxx

            auto g_xy_xxx = contr_buffer_dfxx[df_off + 10 * ccomps * dcomps];

            auto g_xy_xxy = contr_buffer_dfxx[df_off + 11 * ccomps * dcomps];

            auto g_xy_xxz = contr_buffer_dfxx[df_off + 12 * ccomps * dcomps];

            auto g_xy_xyy = contr_buffer_dfxx[df_off + 13 * ccomps * dcomps];

            auto g_xy_xyz = contr_buffer_dfxx[df_off + 14 * ccomps * dcomps];

            auto g_xy_xzz = contr_buffer_dfxx[df_off + 15 * ccomps * dcomps];

            auto g_xy_yyy = contr_buffer_dfxx[df_off + 16 * ccomps * dcomps];

            auto g_xy_yyz = contr_buffer_dfxx[df_off + 17 * ccomps * dcomps];

            auto g_xy_yzz = contr_buffer_dfxx[df_off + 18 * ccomps * dcomps];

            auto g_xy_zzz = contr_buffer_dfxx[df_off + 19 * ccomps * dcomps];

            #pragma omp simd aligned(g_xy_xxx, g_xy_xxy, g_xy_xxz, g_xy_xyy, g_xy_xyz, g_xy_xzz, g_xy_yyy, g_xy_yyz, g_xy_yzz, g_xy_zzz, g_y_xxx, g_y_xxxx, g_y_xxxy, g_y_xxxz, g_y_xxy, g_y_xxyy, g_y_xxyz, g_y_xxz, g_y_xxzz, g_y_xyy, g_y_xyyy, g_y_xyyz, g_y_xyz, g_y_xyzz, g_y_xzz, g_y_xzzz, g_y_yyy, g_y_yyz, g_y_yzz, g_y_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xy_xxx[k] = -g_y_xxx[k] * ab_x + g_y_xxxx[k];

                g_xy_xxy[k] = -g_y_xxy[k] * ab_x + g_y_xxxy[k];

                g_xy_xxz[k] = -g_y_xxz[k] * ab_x + g_y_xxxz[k];

                g_xy_xyy[k] = -g_y_xyy[k] * ab_x + g_y_xxyy[k];

                g_xy_xyz[k] = -g_y_xyz[k] * ab_x + g_y_xxyz[k];

                g_xy_xzz[k] = -g_y_xzz[k] * ab_x + g_y_xxzz[k];

                g_xy_yyy[k] = -g_y_yyy[k] * ab_x + g_y_xyyy[k];

                g_xy_yyz[k] = -g_y_yyz[k] * ab_x + g_y_xyyz[k];

                g_xy_yzz[k] = -g_y_yzz[k] * ab_x + g_y_xyzz[k];

                g_xy_zzz[k] = -g_y_zzz[k] * ab_x + g_y_xzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : contr_buffer_dfxx

            auto g_xz_xxx = contr_buffer_dfxx[df_off + 20 * ccomps * dcomps];

            auto g_xz_xxy = contr_buffer_dfxx[df_off + 21 * ccomps * dcomps];

            auto g_xz_xxz = contr_buffer_dfxx[df_off + 22 * ccomps * dcomps];

            auto g_xz_xyy = contr_buffer_dfxx[df_off + 23 * ccomps * dcomps];

            auto g_xz_xyz = contr_buffer_dfxx[df_off + 24 * ccomps * dcomps];

            auto g_xz_xzz = contr_buffer_dfxx[df_off + 25 * ccomps * dcomps];

            auto g_xz_yyy = contr_buffer_dfxx[df_off + 26 * ccomps * dcomps];

            auto g_xz_yyz = contr_buffer_dfxx[df_off + 27 * ccomps * dcomps];

            auto g_xz_yzz = contr_buffer_dfxx[df_off + 28 * ccomps * dcomps];

            auto g_xz_zzz = contr_buffer_dfxx[df_off + 29 * ccomps * dcomps];

            #pragma omp simd aligned(g_xz_xxx, g_xz_xxy, g_xz_xxz, g_xz_xyy, g_xz_xyz, g_xz_xzz, g_xz_yyy, g_xz_yyz, g_xz_yzz, g_xz_zzz, g_z_xxx, g_z_xxxx, g_z_xxxy, g_z_xxxz, g_z_xxy, g_z_xxyy, g_z_xxyz, g_z_xxz, g_z_xxzz, g_z_xyy, g_z_xyyy, g_z_xyyz, g_z_xyz, g_z_xyzz, g_z_xzz, g_z_xzzz, g_z_yyy, g_z_yyz, g_z_yzz, g_z_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xz_xxx[k] = -g_z_xxx[k] * ab_x + g_z_xxxx[k];

                g_xz_xxy[k] = -g_z_xxy[k] * ab_x + g_z_xxxy[k];

                g_xz_xxz[k] = -g_z_xxz[k] * ab_x + g_z_xxxz[k];

                g_xz_xyy[k] = -g_z_xyy[k] * ab_x + g_z_xxyy[k];

                g_xz_xyz[k] = -g_z_xyz[k] * ab_x + g_z_xxyz[k];

                g_xz_xzz[k] = -g_z_xzz[k] * ab_x + g_z_xxzz[k];

                g_xz_yyy[k] = -g_z_yyy[k] * ab_x + g_z_xyyy[k];

                g_xz_yyz[k] = -g_z_yyz[k] * ab_x + g_z_xyyz[k];

                g_xz_yzz[k] = -g_z_yzz[k] * ab_x + g_z_xyzz[k];

                g_xz_zzz[k] = -g_z_zzz[k] * ab_x + g_z_xzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : contr_buffer_dfxx

            auto g_yy_xxx = contr_buffer_dfxx[df_off + 30 * ccomps * dcomps];

            auto g_yy_xxy = contr_buffer_dfxx[df_off + 31 * ccomps * dcomps];

            auto g_yy_xxz = contr_buffer_dfxx[df_off + 32 * ccomps * dcomps];

            auto g_yy_xyy = contr_buffer_dfxx[df_off + 33 * ccomps * dcomps];

            auto g_yy_xyz = contr_buffer_dfxx[df_off + 34 * ccomps * dcomps];

            auto g_yy_xzz = contr_buffer_dfxx[df_off + 35 * ccomps * dcomps];

            auto g_yy_yyy = contr_buffer_dfxx[df_off + 36 * ccomps * dcomps];

            auto g_yy_yyz = contr_buffer_dfxx[df_off + 37 * ccomps * dcomps];

            auto g_yy_yzz = contr_buffer_dfxx[df_off + 38 * ccomps * dcomps];

            auto g_yy_zzz = contr_buffer_dfxx[df_off + 39 * ccomps * dcomps];

            #pragma omp simd aligned(g_y_xxx, g_y_xxxy, g_y_xxy, g_y_xxyy, g_y_xxyz, g_y_xxz, g_y_xyy, g_y_xyyy, g_y_xyyz, g_y_xyz, g_y_xyzz, g_y_xzz, g_y_yyy, g_y_yyyy, g_y_yyyz, g_y_yyz, g_y_yyzz, g_y_yzz, g_y_yzzz, g_y_zzz, g_yy_xxx, g_yy_xxy, g_yy_xxz, g_yy_xyy, g_yy_xyz, g_yy_xzz, g_yy_yyy, g_yy_yyz, g_yy_yzz, g_yy_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yy_xxx[k] = -g_y_xxx[k] * ab_y + g_y_xxxy[k];

                g_yy_xxy[k] = -g_y_xxy[k] * ab_y + g_y_xxyy[k];

                g_yy_xxz[k] = -g_y_xxz[k] * ab_y + g_y_xxyz[k];

                g_yy_xyy[k] = -g_y_xyy[k] * ab_y + g_y_xyyy[k];

                g_yy_xyz[k] = -g_y_xyz[k] * ab_y + g_y_xyyz[k];

                g_yy_xzz[k] = -g_y_xzz[k] * ab_y + g_y_xyzz[k];

                g_yy_yyy[k] = -g_y_yyy[k] * ab_y + g_y_yyyy[k];

                g_yy_yyz[k] = -g_y_yyz[k] * ab_y + g_y_yyyz[k];

                g_yy_yzz[k] = -g_y_yzz[k] * ab_y + g_y_yyzz[k];

                g_yy_zzz[k] = -g_y_zzz[k] * ab_y + g_y_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : contr_buffer_dfxx

            auto g_yz_xxx = contr_buffer_dfxx[df_off + 40 * ccomps * dcomps];

            auto g_yz_xxy = contr_buffer_dfxx[df_off + 41 * ccomps * dcomps];

            auto g_yz_xxz = contr_buffer_dfxx[df_off + 42 * ccomps * dcomps];

            auto g_yz_xyy = contr_buffer_dfxx[df_off + 43 * ccomps * dcomps];

            auto g_yz_xyz = contr_buffer_dfxx[df_off + 44 * ccomps * dcomps];

            auto g_yz_xzz = contr_buffer_dfxx[df_off + 45 * ccomps * dcomps];

            auto g_yz_yyy = contr_buffer_dfxx[df_off + 46 * ccomps * dcomps];

            auto g_yz_yyz = contr_buffer_dfxx[df_off + 47 * ccomps * dcomps];

            auto g_yz_yzz = contr_buffer_dfxx[df_off + 48 * ccomps * dcomps];

            auto g_yz_zzz = contr_buffer_dfxx[df_off + 49 * ccomps * dcomps];

            #pragma omp simd aligned(g_yz_xxx, g_yz_xxy, g_yz_xxz, g_yz_xyy, g_yz_xyz, g_yz_xzz, g_yz_yyy, g_yz_yyz, g_yz_yzz, g_yz_zzz, g_z_xxx, g_z_xxxy, g_z_xxy, g_z_xxyy, g_z_xxyz, g_z_xxz, g_z_xyy, g_z_xyyy, g_z_xyyz, g_z_xyz, g_z_xyzz, g_z_xzz, g_z_yyy, g_z_yyyy, g_z_yyyz, g_z_yyz, g_z_yyzz, g_z_yzz, g_z_yzzz, g_z_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yz_xxx[k] = -g_z_xxx[k] * ab_y + g_z_xxxy[k];

                g_yz_xxy[k] = -g_z_xxy[k] * ab_y + g_z_xxyy[k];

                g_yz_xxz[k] = -g_z_xxz[k] * ab_y + g_z_xxyz[k];

                g_yz_xyy[k] = -g_z_xyy[k] * ab_y + g_z_xyyy[k];

                g_yz_xyz[k] = -g_z_xyz[k] * ab_y + g_z_xyyz[k];

                g_yz_xzz[k] = -g_z_xzz[k] * ab_y + g_z_xyzz[k];

                g_yz_yyy[k] = -g_z_yyy[k] * ab_y + g_z_yyyy[k];

                g_yz_yyz[k] = -g_z_yyz[k] * ab_y + g_z_yyyz[k];

                g_yz_yzz[k] = -g_z_yzz[k] * ab_y + g_z_yyzz[k];

                g_yz_zzz[k] = -g_z_zzz[k] * ab_y + g_z_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : contr_buffer_dfxx

            auto g_zz_xxx = contr_buffer_dfxx[df_off + 50 * ccomps * dcomps];

            auto g_zz_xxy = contr_buffer_dfxx[df_off + 51 * ccomps * dcomps];

            auto g_zz_xxz = contr_buffer_dfxx[df_off + 52 * ccomps * dcomps];

            auto g_zz_xyy = contr_buffer_dfxx[df_off + 53 * ccomps * dcomps];

            auto g_zz_xyz = contr_buffer_dfxx[df_off + 54 * ccomps * dcomps];

            auto g_zz_xzz = contr_buffer_dfxx[df_off + 55 * ccomps * dcomps];

            auto g_zz_yyy = contr_buffer_dfxx[df_off + 56 * ccomps * dcomps];

            auto g_zz_yyz = contr_buffer_dfxx[df_off + 57 * ccomps * dcomps];

            auto g_zz_yzz = contr_buffer_dfxx[df_off + 58 * ccomps * dcomps];

            auto g_zz_zzz = contr_buffer_dfxx[df_off + 59 * ccomps * dcomps];

            #pragma omp simd aligned(g_z_xxx, g_z_xxxz, g_z_xxy, g_z_xxyz, g_z_xxz, g_z_xxzz, g_z_xyy, g_z_xyyz, g_z_xyz, g_z_xyzz, g_z_xzz, g_z_xzzz, g_z_yyy, g_z_yyyz, g_z_yyz, g_z_yyzz, g_z_yzz, g_z_yzzz, g_z_zzz, g_z_zzzz, g_zz_xxx, g_zz_xxy, g_zz_xxz, g_zz_xyy, g_zz_xyz, g_zz_xzz, g_zz_yyy, g_zz_yyz, g_zz_yzz, g_zz_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_zz_xxx[k] = -g_z_xxx[k] * ab_z + g_z_xxxz[k];

                g_zz_xxy[k] = -g_z_xxy[k] * ab_z + g_z_xxyz[k];

                g_zz_xxz[k] = -g_z_xxz[k] * ab_z + g_z_xxzz[k];

                g_zz_xyy[k] = -g_z_xyy[k] * ab_z + g_z_xyyz[k];

                g_zz_xyz[k] = -g_z_xyz[k] * ab_z + g_z_xyzz[k];

                g_zz_xzz[k] = -g_z_xzz[k] * ab_z + g_z_xzzz[k];

                g_zz_yyy[k] = -g_z_yyy[k] * ab_z + g_z_yyyz[k];

                g_zz_yyz[k] = -g_z_yyz[k] * ab_z + g_z_yyzz[k];

                g_zz_yzz[k] = -g_z_yzz[k] * ab_z + g_z_yzzz[k];

                g_zz_zzz[k] = -g_z_zzz[k] * ab_z + g_z_zzzz[k];
            }
        }
    }
}

} // erirec namespace

