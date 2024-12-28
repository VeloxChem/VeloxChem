#include "ElectronRepulsionContrRecDFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto
comp_bra_hrr_electron_repulsion_dfxx(CSimdArray<double>&   cbuffer,
                                     const size_t          idx_dfxx,
                                     const size_t          idx_pfxx,
                                     const size_t          idx_pgxx,
                                     const TPoint<double>& r_ab,
                                     const int             c_angmom,
                                     const int             d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{
        c_angmom,
    });

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{
        d_angmom,
    });

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : PFSS

            const auto pf_off = idx_pfxx + i * dcomps + j;

            auto g_x_xxx = cbuffer.data(pf_off + 0 * ccomps * dcomps);

            auto g_x_xxy = cbuffer.data(pf_off + 1 * ccomps * dcomps);

            auto g_x_xxz = cbuffer.data(pf_off + 2 * ccomps * dcomps);

            auto g_x_xyy = cbuffer.data(pf_off + 3 * ccomps * dcomps);

            auto g_x_xyz = cbuffer.data(pf_off + 4 * ccomps * dcomps);

            auto g_x_xzz = cbuffer.data(pf_off + 5 * ccomps * dcomps);

            auto g_x_yyy = cbuffer.data(pf_off + 6 * ccomps * dcomps);

            auto g_x_yyz = cbuffer.data(pf_off + 7 * ccomps * dcomps);

            auto g_x_yzz = cbuffer.data(pf_off + 8 * ccomps * dcomps);

            auto g_x_zzz = cbuffer.data(pf_off + 9 * ccomps * dcomps);

            auto g_y_xxx = cbuffer.data(pf_off + 10 * ccomps * dcomps);

            auto g_y_xxy = cbuffer.data(pf_off + 11 * ccomps * dcomps);

            auto g_y_xxz = cbuffer.data(pf_off + 12 * ccomps * dcomps);

            auto g_y_xyy = cbuffer.data(pf_off + 13 * ccomps * dcomps);

            auto g_y_xyz = cbuffer.data(pf_off + 14 * ccomps * dcomps);

            auto g_y_xzz = cbuffer.data(pf_off + 15 * ccomps * dcomps);

            auto g_y_yyy = cbuffer.data(pf_off + 16 * ccomps * dcomps);

            auto g_y_yyz = cbuffer.data(pf_off + 17 * ccomps * dcomps);

            auto g_y_yzz = cbuffer.data(pf_off + 18 * ccomps * dcomps);

            auto g_y_zzz = cbuffer.data(pf_off + 19 * ccomps * dcomps);

            auto g_z_xxx = cbuffer.data(pf_off + 20 * ccomps * dcomps);

            auto g_z_xxy = cbuffer.data(pf_off + 21 * ccomps * dcomps);

            auto g_z_xxz = cbuffer.data(pf_off + 22 * ccomps * dcomps);

            auto g_z_xyy = cbuffer.data(pf_off + 23 * ccomps * dcomps);

            auto g_z_xyz = cbuffer.data(pf_off + 24 * ccomps * dcomps);

            auto g_z_xzz = cbuffer.data(pf_off + 25 * ccomps * dcomps);

            auto g_z_yyy = cbuffer.data(pf_off + 26 * ccomps * dcomps);

            auto g_z_yyz = cbuffer.data(pf_off + 27 * ccomps * dcomps);

            auto g_z_yzz = cbuffer.data(pf_off + 28 * ccomps * dcomps);

            auto g_z_zzz = cbuffer.data(pf_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PGSS

            const auto pg_off = idx_pgxx + i * dcomps + j;

            auto g_x_xxxx = cbuffer.data(pg_off + 0 * ccomps * dcomps);

            auto g_x_xxxy = cbuffer.data(pg_off + 1 * ccomps * dcomps);

            auto g_x_xxxz = cbuffer.data(pg_off + 2 * ccomps * dcomps);

            auto g_x_xxyy = cbuffer.data(pg_off + 3 * ccomps * dcomps);

            auto g_x_xxyz = cbuffer.data(pg_off + 4 * ccomps * dcomps);

            auto g_x_xxzz = cbuffer.data(pg_off + 5 * ccomps * dcomps);

            auto g_x_xyyy = cbuffer.data(pg_off + 6 * ccomps * dcomps);

            auto g_x_xyyz = cbuffer.data(pg_off + 7 * ccomps * dcomps);

            auto g_x_xyzz = cbuffer.data(pg_off + 8 * ccomps * dcomps);

            auto g_x_xzzz = cbuffer.data(pg_off + 9 * ccomps * dcomps);

            auto g_y_xxxx = cbuffer.data(pg_off + 15 * ccomps * dcomps);

            auto g_y_xxxy = cbuffer.data(pg_off + 16 * ccomps * dcomps);

            auto g_y_xxxz = cbuffer.data(pg_off + 17 * ccomps * dcomps);

            auto g_y_xxyy = cbuffer.data(pg_off + 18 * ccomps * dcomps);

            auto g_y_xxyz = cbuffer.data(pg_off + 19 * ccomps * dcomps);

            auto g_y_xxzz = cbuffer.data(pg_off + 20 * ccomps * dcomps);

            auto g_y_xyyy = cbuffer.data(pg_off + 21 * ccomps * dcomps);

            auto g_y_xyyz = cbuffer.data(pg_off + 22 * ccomps * dcomps);

            auto g_y_xyzz = cbuffer.data(pg_off + 23 * ccomps * dcomps);

            auto g_y_xzzz = cbuffer.data(pg_off + 24 * ccomps * dcomps);

            auto g_y_yyyy = cbuffer.data(pg_off + 25 * ccomps * dcomps);

            auto g_y_yyyz = cbuffer.data(pg_off + 26 * ccomps * dcomps);

            auto g_y_yyzz = cbuffer.data(pg_off + 27 * ccomps * dcomps);

            auto g_y_yzzz = cbuffer.data(pg_off + 28 * ccomps * dcomps);

            auto g_z_xxxx = cbuffer.data(pg_off + 30 * ccomps * dcomps);

            auto g_z_xxxy = cbuffer.data(pg_off + 31 * ccomps * dcomps);

            auto g_z_xxxz = cbuffer.data(pg_off + 32 * ccomps * dcomps);

            auto g_z_xxyy = cbuffer.data(pg_off + 33 * ccomps * dcomps);

            auto g_z_xxyz = cbuffer.data(pg_off + 34 * ccomps * dcomps);

            auto g_z_xxzz = cbuffer.data(pg_off + 35 * ccomps * dcomps);

            auto g_z_xyyy = cbuffer.data(pg_off + 36 * ccomps * dcomps);

            auto g_z_xyyz = cbuffer.data(pg_off + 37 * ccomps * dcomps);

            auto g_z_xyzz = cbuffer.data(pg_off + 38 * ccomps * dcomps);

            auto g_z_xzzz = cbuffer.data(pg_off + 39 * ccomps * dcomps);

            auto g_z_yyyy = cbuffer.data(pg_off + 40 * ccomps * dcomps);

            auto g_z_yyyz = cbuffer.data(pg_off + 41 * ccomps * dcomps);

            auto g_z_yyzz = cbuffer.data(pg_off + 42 * ccomps * dcomps);

            auto g_z_yzzz = cbuffer.data(pg_off + 43 * ccomps * dcomps);

            auto g_z_zzzz = cbuffer.data(pg_off + 44 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dfxx

            const auto df_off = idx_dfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_xx_xxx = cbuffer.data(df_off + 0 * ccomps * dcomps);

            auto g_xx_xxy = cbuffer.data(df_off + 1 * ccomps * dcomps);

            auto g_xx_xxz = cbuffer.data(df_off + 2 * ccomps * dcomps);

            auto g_xx_xyy = cbuffer.data(df_off + 3 * ccomps * dcomps);

            auto g_xx_xyz = cbuffer.data(df_off + 4 * ccomps * dcomps);

            auto g_xx_xzz = cbuffer.data(df_off + 5 * ccomps * dcomps);

            auto g_xx_yyy = cbuffer.data(df_off + 6 * ccomps * dcomps);

            auto g_xx_yyz = cbuffer.data(df_off + 7 * ccomps * dcomps);

            auto g_xx_yzz = cbuffer.data(df_off + 8 * ccomps * dcomps);

            auto g_xx_zzz = cbuffer.data(df_off + 9 * ccomps * dcomps);

#pragma omp simd aligned(g_x_xxx,      \
                             g_x_xxxx, \
                             g_x_xxxy, \
                             g_x_xxxz, \
                             g_x_xxy,  \
                             g_x_xxyy, \
                             g_x_xxyz, \
                             g_x_xxz,  \
                             g_x_xxzz, \
                             g_x_xyy,  \
                             g_x_xyyy, \
                             g_x_xyyz, \
                             g_x_xyz,  \
                             g_x_xyzz, \
                             g_x_xzz,  \
                             g_x_xzzz, \
                             g_x_yyy,  \
                             g_x_yyz,  \
                             g_x_yzz,  \
                             g_x_zzz,  \
                             g_xx_xxx, \
                             g_xx_xxy, \
                             g_xx_xxz, \
                             g_xx_xyy, \
                             g_xx_xyz, \
                             g_xx_xzz, \
                             g_xx_yyy, \
                             g_xx_yyz, \
                             g_xx_yzz, \
                             g_xx_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_xy_xxx = cbuffer.data(df_off + 10 * ccomps * dcomps);

            auto g_xy_xxy = cbuffer.data(df_off + 11 * ccomps * dcomps);

            auto g_xy_xxz = cbuffer.data(df_off + 12 * ccomps * dcomps);

            auto g_xy_xyy = cbuffer.data(df_off + 13 * ccomps * dcomps);

            auto g_xy_xyz = cbuffer.data(df_off + 14 * ccomps * dcomps);

            auto g_xy_xzz = cbuffer.data(df_off + 15 * ccomps * dcomps);

            auto g_xy_yyy = cbuffer.data(df_off + 16 * ccomps * dcomps);

            auto g_xy_yyz = cbuffer.data(df_off + 17 * ccomps * dcomps);

            auto g_xy_yzz = cbuffer.data(df_off + 18 * ccomps * dcomps);

            auto g_xy_zzz = cbuffer.data(df_off + 19 * ccomps * dcomps);

#pragma omp simd aligned(g_xy_xxx,     \
                             g_xy_xxy, \
                             g_xy_xxz, \
                             g_xy_xyy, \
                             g_xy_xyz, \
                             g_xy_xzz, \
                             g_xy_yyy, \
                             g_xy_yyz, \
                             g_xy_yzz, \
                             g_xy_zzz, \
                             g_y_xxx,  \
                             g_y_xxxx, \
                             g_y_xxxy, \
                             g_y_xxxz, \
                             g_y_xxy,  \
                             g_y_xxyy, \
                             g_y_xxyz, \
                             g_y_xxz,  \
                             g_y_xxzz, \
                             g_y_xyy,  \
                             g_y_xyyy, \
                             g_y_xyyz, \
                             g_y_xyz,  \
                             g_y_xyzz, \
                             g_y_xzz,  \
                             g_y_xzzz, \
                             g_y_yyy,  \
                             g_y_yyz,  \
                             g_y_yzz,  \
                             g_y_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_xz_xxx = cbuffer.data(df_off + 20 * ccomps * dcomps);

            auto g_xz_xxy = cbuffer.data(df_off + 21 * ccomps * dcomps);

            auto g_xz_xxz = cbuffer.data(df_off + 22 * ccomps * dcomps);

            auto g_xz_xyy = cbuffer.data(df_off + 23 * ccomps * dcomps);

            auto g_xz_xyz = cbuffer.data(df_off + 24 * ccomps * dcomps);

            auto g_xz_xzz = cbuffer.data(df_off + 25 * ccomps * dcomps);

            auto g_xz_yyy = cbuffer.data(df_off + 26 * ccomps * dcomps);

            auto g_xz_yyz = cbuffer.data(df_off + 27 * ccomps * dcomps);

            auto g_xz_yzz = cbuffer.data(df_off + 28 * ccomps * dcomps);

            auto g_xz_zzz = cbuffer.data(df_off + 29 * ccomps * dcomps);

#pragma omp simd aligned(g_xz_xxx,     \
                             g_xz_xxy, \
                             g_xz_xxz, \
                             g_xz_xyy, \
                             g_xz_xyz, \
                             g_xz_xzz, \
                             g_xz_yyy, \
                             g_xz_yyz, \
                             g_xz_yzz, \
                             g_xz_zzz, \
                             g_z_xxx,  \
                             g_z_xxxx, \
                             g_z_xxxy, \
                             g_z_xxxz, \
                             g_z_xxy,  \
                             g_z_xxyy, \
                             g_z_xxyz, \
                             g_z_xxz,  \
                             g_z_xxzz, \
                             g_z_xyy,  \
                             g_z_xyyy, \
                             g_z_xyyz, \
                             g_z_xyz,  \
                             g_z_xyzz, \
                             g_z_xzz,  \
                             g_z_xzzz, \
                             g_z_yyy,  \
                             g_z_yyz,  \
                             g_z_yzz,  \
                             g_z_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_yy_xxx = cbuffer.data(df_off + 30 * ccomps * dcomps);

            auto g_yy_xxy = cbuffer.data(df_off + 31 * ccomps * dcomps);

            auto g_yy_xxz = cbuffer.data(df_off + 32 * ccomps * dcomps);

            auto g_yy_xyy = cbuffer.data(df_off + 33 * ccomps * dcomps);

            auto g_yy_xyz = cbuffer.data(df_off + 34 * ccomps * dcomps);

            auto g_yy_xzz = cbuffer.data(df_off + 35 * ccomps * dcomps);

            auto g_yy_yyy = cbuffer.data(df_off + 36 * ccomps * dcomps);

            auto g_yy_yyz = cbuffer.data(df_off + 37 * ccomps * dcomps);

            auto g_yy_yzz = cbuffer.data(df_off + 38 * ccomps * dcomps);

            auto g_yy_zzz = cbuffer.data(df_off + 39 * ccomps * dcomps);

#pragma omp simd aligned(g_y_xxx,      \
                             g_y_xxxy, \
                             g_y_xxy,  \
                             g_y_xxyy, \
                             g_y_xxyz, \
                             g_y_xxz,  \
                             g_y_xyy,  \
                             g_y_xyyy, \
                             g_y_xyyz, \
                             g_y_xyz,  \
                             g_y_xyzz, \
                             g_y_xzz,  \
                             g_y_yyy,  \
                             g_y_yyyy, \
                             g_y_yyyz, \
                             g_y_yyz,  \
                             g_y_yyzz, \
                             g_y_yzz,  \
                             g_y_yzzz, \
                             g_y_zzz,  \
                             g_yy_xxx, \
                             g_yy_xxy, \
                             g_yy_xxz, \
                             g_yy_xyy, \
                             g_yy_xyz, \
                             g_yy_xzz, \
                             g_yy_yyy, \
                             g_yy_yyz, \
                             g_yy_yzz, \
                             g_yy_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_yz_xxx = cbuffer.data(df_off + 40 * ccomps * dcomps);

            auto g_yz_xxy = cbuffer.data(df_off + 41 * ccomps * dcomps);

            auto g_yz_xxz = cbuffer.data(df_off + 42 * ccomps * dcomps);

            auto g_yz_xyy = cbuffer.data(df_off + 43 * ccomps * dcomps);

            auto g_yz_xyz = cbuffer.data(df_off + 44 * ccomps * dcomps);

            auto g_yz_xzz = cbuffer.data(df_off + 45 * ccomps * dcomps);

            auto g_yz_yyy = cbuffer.data(df_off + 46 * ccomps * dcomps);

            auto g_yz_yyz = cbuffer.data(df_off + 47 * ccomps * dcomps);

            auto g_yz_yzz = cbuffer.data(df_off + 48 * ccomps * dcomps);

            auto g_yz_zzz = cbuffer.data(df_off + 49 * ccomps * dcomps);

#pragma omp simd aligned(g_yz_xxx,     \
                             g_yz_xxy, \
                             g_yz_xxz, \
                             g_yz_xyy, \
                             g_yz_xyz, \
                             g_yz_xzz, \
                             g_yz_yyy, \
                             g_yz_yyz, \
                             g_yz_yzz, \
                             g_yz_zzz, \
                             g_z_xxx,  \
                             g_z_xxxy, \
                             g_z_xxy,  \
                             g_z_xxyy, \
                             g_z_xxyz, \
                             g_z_xxz,  \
                             g_z_xyy,  \
                             g_z_xyyy, \
                             g_z_xyyz, \
                             g_z_xyz,  \
                             g_z_xyzz, \
                             g_z_xzz,  \
                             g_z_yyy,  \
                             g_z_yyyy, \
                             g_z_yyyz, \
                             g_z_yyz,  \
                             g_z_yyzz, \
                             g_z_yzz,  \
                             g_z_yzzz, \
                             g_z_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_zz_xxx = cbuffer.data(df_off + 50 * ccomps * dcomps);

            auto g_zz_xxy = cbuffer.data(df_off + 51 * ccomps * dcomps);

            auto g_zz_xxz = cbuffer.data(df_off + 52 * ccomps * dcomps);

            auto g_zz_xyy = cbuffer.data(df_off + 53 * ccomps * dcomps);

            auto g_zz_xyz = cbuffer.data(df_off + 54 * ccomps * dcomps);

            auto g_zz_xzz = cbuffer.data(df_off + 55 * ccomps * dcomps);

            auto g_zz_yyy = cbuffer.data(df_off + 56 * ccomps * dcomps);

            auto g_zz_yyz = cbuffer.data(df_off + 57 * ccomps * dcomps);

            auto g_zz_yzz = cbuffer.data(df_off + 58 * ccomps * dcomps);

            auto g_zz_zzz = cbuffer.data(df_off + 59 * ccomps * dcomps);

#pragma omp simd aligned(g_z_xxx,      \
                             g_z_xxxz, \
                             g_z_xxy,  \
                             g_z_xxyz, \
                             g_z_xxz,  \
                             g_z_xxzz, \
                             g_z_xyy,  \
                             g_z_xyyz, \
                             g_z_xyz,  \
                             g_z_xyzz, \
                             g_z_xzz,  \
                             g_z_xzzz, \
                             g_z_yyy,  \
                             g_z_yyyz, \
                             g_z_yyz,  \
                             g_z_yyzz, \
                             g_z_yzz,  \
                             g_z_yzzz, \
                             g_z_zzz,  \
                             g_z_zzzz, \
                             g_zz_xxx, \
                             g_zz_xxy, \
                             g_zz_xxz, \
                             g_zz_xyy, \
                             g_zz_xyz, \
                             g_zz_xzz, \
                             g_zz_yyy, \
                             g_zz_yyz, \
                             g_zz_yzz, \
                             g_zz_zzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

}  // namespace erirec