#include "ElectronRepulsionContrRecFDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_fdxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_fdxx,
                                     const size_t idx_ddxx,
                                     const size_t idx_dfxx,
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
            /// Set up components of auxilary buffer : DDSS

            const auto dd_off = idx_ddxx + i * dcomps + j;

            auto g_xx_xx = cbuffer.data(dd_off + 0 * ccomps * dcomps);

            auto g_xx_xy = cbuffer.data(dd_off + 1 * ccomps * dcomps);

            auto g_xx_xz = cbuffer.data(dd_off + 2 * ccomps * dcomps);

            auto g_xx_yy = cbuffer.data(dd_off + 3 * ccomps * dcomps);

            auto g_xx_yz = cbuffer.data(dd_off + 4 * ccomps * dcomps);

            auto g_xx_zz = cbuffer.data(dd_off + 5 * ccomps * dcomps);

            auto g_xy_xx = cbuffer.data(dd_off + 6 * ccomps * dcomps);

            auto g_xy_xy = cbuffer.data(dd_off + 7 * ccomps * dcomps);

            auto g_xy_xz = cbuffer.data(dd_off + 8 * ccomps * dcomps);

            auto g_xy_yy = cbuffer.data(dd_off + 9 * ccomps * dcomps);

            auto g_xy_yz = cbuffer.data(dd_off + 10 * ccomps * dcomps);

            auto g_xy_zz = cbuffer.data(dd_off + 11 * ccomps * dcomps);

            auto g_xz_xx = cbuffer.data(dd_off + 12 * ccomps * dcomps);

            auto g_xz_xy = cbuffer.data(dd_off + 13 * ccomps * dcomps);

            auto g_xz_xz = cbuffer.data(dd_off + 14 * ccomps * dcomps);

            auto g_xz_yy = cbuffer.data(dd_off + 15 * ccomps * dcomps);

            auto g_xz_yz = cbuffer.data(dd_off + 16 * ccomps * dcomps);

            auto g_xz_zz = cbuffer.data(dd_off + 17 * ccomps * dcomps);

            auto g_yy_xx = cbuffer.data(dd_off + 18 * ccomps * dcomps);

            auto g_yy_xy = cbuffer.data(dd_off + 19 * ccomps * dcomps);

            auto g_yy_xz = cbuffer.data(dd_off + 20 * ccomps * dcomps);

            auto g_yy_yy = cbuffer.data(dd_off + 21 * ccomps * dcomps);

            auto g_yy_yz = cbuffer.data(dd_off + 22 * ccomps * dcomps);

            auto g_yy_zz = cbuffer.data(dd_off + 23 * ccomps * dcomps);

            auto g_yz_xx = cbuffer.data(dd_off + 24 * ccomps * dcomps);

            auto g_yz_xy = cbuffer.data(dd_off + 25 * ccomps * dcomps);

            auto g_yz_xz = cbuffer.data(dd_off + 26 * ccomps * dcomps);

            auto g_yz_yy = cbuffer.data(dd_off + 27 * ccomps * dcomps);

            auto g_yz_yz = cbuffer.data(dd_off + 28 * ccomps * dcomps);

            auto g_yz_zz = cbuffer.data(dd_off + 29 * ccomps * dcomps);

            auto g_zz_xx = cbuffer.data(dd_off + 30 * ccomps * dcomps);

            auto g_zz_xy = cbuffer.data(dd_off + 31 * ccomps * dcomps);

            auto g_zz_xz = cbuffer.data(dd_off + 32 * ccomps * dcomps);

            auto g_zz_yy = cbuffer.data(dd_off + 33 * ccomps * dcomps);

            auto g_zz_yz = cbuffer.data(dd_off + 34 * ccomps * dcomps);

            auto g_zz_zz = cbuffer.data(dd_off + 35 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DFSS

            const auto df_off = idx_dfxx + i * dcomps + j;

            auto g_xx_xxx = cbuffer.data(df_off + 0 * ccomps * dcomps);

            auto g_xx_xxy = cbuffer.data(df_off + 1 * ccomps * dcomps);

            auto g_xx_xxz = cbuffer.data(df_off + 2 * ccomps * dcomps);

            auto g_xx_xyy = cbuffer.data(df_off + 3 * ccomps * dcomps);

            auto g_xx_xyz = cbuffer.data(df_off + 4 * ccomps * dcomps);

            auto g_xx_xzz = cbuffer.data(df_off + 5 * ccomps * dcomps);

            auto g_xy_xxx = cbuffer.data(df_off + 10 * ccomps * dcomps);

            auto g_xy_xxy = cbuffer.data(df_off + 11 * ccomps * dcomps);

            auto g_xy_xxz = cbuffer.data(df_off + 12 * ccomps * dcomps);

            auto g_xy_xyy = cbuffer.data(df_off + 13 * ccomps * dcomps);

            auto g_xy_xyz = cbuffer.data(df_off + 14 * ccomps * dcomps);

            auto g_xy_xzz = cbuffer.data(df_off + 15 * ccomps * dcomps);

            auto g_xz_xxx = cbuffer.data(df_off + 20 * ccomps * dcomps);

            auto g_xz_xxy = cbuffer.data(df_off + 21 * ccomps * dcomps);

            auto g_xz_xxz = cbuffer.data(df_off + 22 * ccomps * dcomps);

            auto g_xz_xyy = cbuffer.data(df_off + 23 * ccomps * dcomps);

            auto g_xz_xyz = cbuffer.data(df_off + 24 * ccomps * dcomps);

            auto g_xz_xzz = cbuffer.data(df_off + 25 * ccomps * dcomps);

            auto g_yy_xxx = cbuffer.data(df_off + 30 * ccomps * dcomps);

            auto g_yy_xxy = cbuffer.data(df_off + 31 * ccomps * dcomps);

            auto g_yy_xxz = cbuffer.data(df_off + 32 * ccomps * dcomps);

            auto g_yy_xyy = cbuffer.data(df_off + 33 * ccomps * dcomps);

            auto g_yy_xyz = cbuffer.data(df_off + 34 * ccomps * dcomps);

            auto g_yy_xzz = cbuffer.data(df_off + 35 * ccomps * dcomps);

            auto g_yy_yyy = cbuffer.data(df_off + 36 * ccomps * dcomps);

            auto g_yy_yyz = cbuffer.data(df_off + 37 * ccomps * dcomps);

            auto g_yy_yzz = cbuffer.data(df_off + 38 * ccomps * dcomps);

            auto g_yz_xxx = cbuffer.data(df_off + 40 * ccomps * dcomps);

            auto g_yz_xxy = cbuffer.data(df_off + 41 * ccomps * dcomps);

            auto g_yz_xxz = cbuffer.data(df_off + 42 * ccomps * dcomps);

            auto g_yz_xyy = cbuffer.data(df_off + 43 * ccomps * dcomps);

            auto g_yz_xyz = cbuffer.data(df_off + 44 * ccomps * dcomps);

            auto g_yz_xzz = cbuffer.data(df_off + 45 * ccomps * dcomps);

            auto g_yz_yyy = cbuffer.data(df_off + 46 * ccomps * dcomps);

            auto g_yz_yyz = cbuffer.data(df_off + 47 * ccomps * dcomps);

            auto g_yz_yzz = cbuffer.data(df_off + 48 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_fdxx

            const auto fd_off = idx_fdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_xxx_xx = cbuffer.data(fd_off + 0 * ccomps * dcomps);

            auto g_xxx_xy = cbuffer.data(fd_off + 1 * ccomps * dcomps);

            auto g_xxx_xz = cbuffer.data(fd_off + 2 * ccomps * dcomps);

            auto g_xxx_yy = cbuffer.data(fd_off + 3 * ccomps * dcomps);

            auto g_xxx_yz = cbuffer.data(fd_off + 4 * ccomps * dcomps);

            auto g_xxx_zz = cbuffer.data(fd_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_xx, g_xx_xxx, g_xx_xxy, g_xx_xxz, g_xx_xy, g_xx_xyy, g_xx_xyz, g_xx_xz, g_xx_xzz, g_xx_yy, g_xx_yz, g_xx_zz, g_xxx_xx, g_xxx_xy, g_xxx_xz, g_xxx_yy, g_xxx_yz, g_xxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxx_xx[k] = -g_xx_xx[k] * ab_x + g_xx_xxx[k];

                g_xxx_xy[k] = -g_xx_xy[k] * ab_x + g_xx_xxy[k];

                g_xxx_xz[k] = -g_xx_xz[k] * ab_x + g_xx_xxz[k];

                g_xxx_yy[k] = -g_xx_yy[k] * ab_x + g_xx_xyy[k];

                g_xxx_yz[k] = -g_xx_yz[k] * ab_x + g_xx_xyz[k];

                g_xxx_zz[k] = -g_xx_zz[k] * ab_x + g_xx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_xxy_xx = cbuffer.data(fd_off + 6 * ccomps * dcomps);

            auto g_xxy_xy = cbuffer.data(fd_off + 7 * ccomps * dcomps);

            auto g_xxy_xz = cbuffer.data(fd_off + 8 * ccomps * dcomps);

            auto g_xxy_yy = cbuffer.data(fd_off + 9 * ccomps * dcomps);

            auto g_xxy_yz = cbuffer.data(fd_off + 10 * ccomps * dcomps);

            auto g_xxy_zz = cbuffer.data(fd_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxy_xx, g_xxy_xy, g_xxy_xz, g_xxy_yy, g_xxy_yz, g_xxy_zz, g_xy_xx, g_xy_xxx, g_xy_xxy, g_xy_xxz, g_xy_xy, g_xy_xyy, g_xy_xyz, g_xy_xz, g_xy_xzz, g_xy_yy, g_xy_yz, g_xy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxy_xx[k] = -g_xy_xx[k] * ab_x + g_xy_xxx[k];

                g_xxy_xy[k] = -g_xy_xy[k] * ab_x + g_xy_xxy[k];

                g_xxy_xz[k] = -g_xy_xz[k] * ab_x + g_xy_xxz[k];

                g_xxy_yy[k] = -g_xy_yy[k] * ab_x + g_xy_xyy[k];

                g_xxy_yz[k] = -g_xy_yz[k] * ab_x + g_xy_xyz[k];

                g_xxy_zz[k] = -g_xy_zz[k] * ab_x + g_xy_xzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_xxz_xx = cbuffer.data(fd_off + 12 * ccomps * dcomps);

            auto g_xxz_xy = cbuffer.data(fd_off + 13 * ccomps * dcomps);

            auto g_xxz_xz = cbuffer.data(fd_off + 14 * ccomps * dcomps);

            auto g_xxz_yy = cbuffer.data(fd_off + 15 * ccomps * dcomps);

            auto g_xxz_yz = cbuffer.data(fd_off + 16 * ccomps * dcomps);

            auto g_xxz_zz = cbuffer.data(fd_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxz_xx, g_xxz_xy, g_xxz_xz, g_xxz_yy, g_xxz_yz, g_xxz_zz, g_xz_xx, g_xz_xxx, g_xz_xxy, g_xz_xxz, g_xz_xy, g_xz_xyy, g_xz_xyz, g_xz_xz, g_xz_xzz, g_xz_yy, g_xz_yz, g_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxz_xx[k] = -g_xz_xx[k] * ab_x + g_xz_xxx[k];

                g_xxz_xy[k] = -g_xz_xy[k] * ab_x + g_xz_xxy[k];

                g_xxz_xz[k] = -g_xz_xz[k] * ab_x + g_xz_xxz[k];

                g_xxz_yy[k] = -g_xz_yy[k] * ab_x + g_xz_xyy[k];

                g_xxz_yz[k] = -g_xz_yz[k] * ab_x + g_xz_xyz[k];

                g_xxz_zz[k] = -g_xz_zz[k] * ab_x + g_xz_xzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_xyy_xx = cbuffer.data(fd_off + 18 * ccomps * dcomps);

            auto g_xyy_xy = cbuffer.data(fd_off + 19 * ccomps * dcomps);

            auto g_xyy_xz = cbuffer.data(fd_off + 20 * ccomps * dcomps);

            auto g_xyy_yy = cbuffer.data(fd_off + 21 * ccomps * dcomps);

            auto g_xyy_yz = cbuffer.data(fd_off + 22 * ccomps * dcomps);

            auto g_xyy_zz = cbuffer.data(fd_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyy_xx, g_xyy_xy, g_xyy_xz, g_xyy_yy, g_xyy_yz, g_xyy_zz, g_yy_xx, g_yy_xxx, g_yy_xxy, g_yy_xxz, g_yy_xy, g_yy_xyy, g_yy_xyz, g_yy_xz, g_yy_xzz, g_yy_yy, g_yy_yz, g_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyy_xx[k] = -g_yy_xx[k] * ab_x + g_yy_xxx[k];

                g_xyy_xy[k] = -g_yy_xy[k] * ab_x + g_yy_xxy[k];

                g_xyy_xz[k] = -g_yy_xz[k] * ab_x + g_yy_xxz[k];

                g_xyy_yy[k] = -g_yy_yy[k] * ab_x + g_yy_xyy[k];

                g_xyy_yz[k] = -g_yy_yz[k] * ab_x + g_yy_xyz[k];

                g_xyy_zz[k] = -g_yy_zz[k] * ab_x + g_yy_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_xyz_xx = cbuffer.data(fd_off + 24 * ccomps * dcomps);

            auto g_xyz_xy = cbuffer.data(fd_off + 25 * ccomps * dcomps);

            auto g_xyz_xz = cbuffer.data(fd_off + 26 * ccomps * dcomps);

            auto g_xyz_yy = cbuffer.data(fd_off + 27 * ccomps * dcomps);

            auto g_xyz_yz = cbuffer.data(fd_off + 28 * ccomps * dcomps);

            auto g_xyz_zz = cbuffer.data(fd_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyz_xx, g_xyz_xy, g_xyz_xz, g_xyz_yy, g_xyz_yz, g_xyz_zz, g_yz_xx, g_yz_xxx, g_yz_xxy, g_yz_xxz, g_yz_xy, g_yz_xyy, g_yz_xyz, g_yz_xz, g_yz_xzz, g_yz_yy, g_yz_yz, g_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyz_xx[k] = -g_yz_xx[k] * ab_x + g_yz_xxx[k];

                g_xyz_xy[k] = -g_yz_xy[k] * ab_x + g_yz_xxy[k];

                g_xyz_xz[k] = -g_yz_xz[k] * ab_x + g_yz_xxz[k];

                g_xyz_yy[k] = -g_yz_yy[k] * ab_x + g_yz_xyy[k];

                g_xyz_yz[k] = -g_yz_yz[k] * ab_x + g_yz_xyz[k];

                g_xyz_zz[k] = -g_yz_zz[k] * ab_x + g_yz_xzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_xzz_xx = cbuffer.data(fd_off + 30 * ccomps * dcomps);

            auto g_xzz_xy = cbuffer.data(fd_off + 31 * ccomps * dcomps);

            auto g_xzz_xz = cbuffer.data(fd_off + 32 * ccomps * dcomps);

            auto g_xzz_yy = cbuffer.data(fd_off + 33 * ccomps * dcomps);

            auto g_xzz_yz = cbuffer.data(fd_off + 34 * ccomps * dcomps);

            auto g_xzz_zz = cbuffer.data(fd_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xzz_xx, g_xzz_xy, g_xzz_xz, g_xzz_yy, g_xzz_yz, g_xzz_zz, g_zz_xx, g_zz_xxx, g_zz_xxy, g_zz_xxz, g_zz_xy, g_zz_xyy, g_zz_xyz, g_zz_xz, g_zz_xzz, g_zz_yy, g_zz_yz, g_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xzz_xx[k] = -g_zz_xx[k] * ab_x + g_zz_xxx[k];

                g_xzz_xy[k] = -g_zz_xy[k] * ab_x + g_zz_xxy[k];

                g_xzz_xz[k] = -g_zz_xz[k] * ab_x + g_zz_xxz[k];

                g_xzz_yy[k] = -g_zz_yy[k] * ab_x + g_zz_xyy[k];

                g_xzz_yz[k] = -g_zz_yz[k] * ab_x + g_zz_xyz[k];

                g_xzz_zz[k] = -g_zz_zz[k] * ab_x + g_zz_xzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_yyy_xx = cbuffer.data(fd_off + 36 * ccomps * dcomps);

            auto g_yyy_xy = cbuffer.data(fd_off + 37 * ccomps * dcomps);

            auto g_yyy_xz = cbuffer.data(fd_off + 38 * ccomps * dcomps);

            auto g_yyy_yy = cbuffer.data(fd_off + 39 * ccomps * dcomps);

            auto g_yyy_yz = cbuffer.data(fd_off + 40 * ccomps * dcomps);

            auto g_yyy_zz = cbuffer.data(fd_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_xx, g_yy_xxy, g_yy_xy, g_yy_xyy, g_yy_xyz, g_yy_xz, g_yy_yy, g_yy_yyy, g_yy_yyz, g_yy_yz, g_yy_yzz, g_yy_zz, g_yyy_xx, g_yyy_xy, g_yyy_xz, g_yyy_yy, g_yyy_yz, g_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyy_xx[k] = -g_yy_xx[k] * ab_y + g_yy_xxy[k];

                g_yyy_xy[k] = -g_yy_xy[k] * ab_y + g_yy_xyy[k];

                g_yyy_xz[k] = -g_yy_xz[k] * ab_y + g_yy_xyz[k];

                g_yyy_yy[k] = -g_yy_yy[k] * ab_y + g_yy_yyy[k];

                g_yyy_yz[k] = -g_yy_yz[k] * ab_y + g_yy_yyz[k];

                g_yyy_zz[k] = -g_yy_zz[k] * ab_y + g_yy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_yyz_xx = cbuffer.data(fd_off + 42 * ccomps * dcomps);

            auto g_yyz_xy = cbuffer.data(fd_off + 43 * ccomps * dcomps);

            auto g_yyz_xz = cbuffer.data(fd_off + 44 * ccomps * dcomps);

            auto g_yyz_yy = cbuffer.data(fd_off + 45 * ccomps * dcomps);

            auto g_yyz_yz = cbuffer.data(fd_off + 46 * ccomps * dcomps);

            auto g_yyz_zz = cbuffer.data(fd_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyz_xx, g_yyz_xy, g_yyz_xz, g_yyz_yy, g_yyz_yz, g_yyz_zz, g_yz_xx, g_yz_xxy, g_yz_xy, g_yz_xyy, g_yz_xyz, g_yz_xz, g_yz_yy, g_yz_yyy, g_yz_yyz, g_yz_yz, g_yz_yzz, g_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyz_xx[k] = -g_yz_xx[k] * ab_y + g_yz_xxy[k];

                g_yyz_xy[k] = -g_yz_xy[k] * ab_y + g_yz_xyy[k];

                g_yyz_xz[k] = -g_yz_xz[k] * ab_y + g_yz_xyz[k];

                g_yyz_yy[k] = -g_yz_yy[k] * ab_y + g_yz_yyy[k];

                g_yyz_yz[k] = -g_yz_yz[k] * ab_y + g_yz_yyz[k];

                g_yyz_zz[k] = -g_yz_zz[k] * ab_y + g_yz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_yzz_xx = cbuffer.data(fd_off + 48 * ccomps * dcomps);

            auto g_yzz_xy = cbuffer.data(fd_off + 49 * ccomps * dcomps);

            auto g_yzz_xz = cbuffer.data(fd_off + 50 * ccomps * dcomps);

            auto g_yzz_yy = cbuffer.data(fd_off + 51 * ccomps * dcomps);

            auto g_yzz_yz = cbuffer.data(fd_off + 52 * ccomps * dcomps);

            auto g_yzz_zz = cbuffer.data(fd_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_yzz_xx, g_yzz_xy, g_yzz_xz, g_yzz_yy, g_yzz_yz, g_yzz_zz, g_zz_xx, g_zz_xxy, g_zz_xy, g_zz_xyy, g_zz_xyz, g_zz_xz, g_zz_yy, g_zz_yyy, g_zz_yyz, g_zz_yz, g_zz_yzz, g_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yzz_xx[k] = -g_zz_xx[k] * ab_y + g_zz_xxy[k];

                g_yzz_xy[k] = -g_zz_xy[k] * ab_y + g_zz_xyy[k];

                g_yzz_xz[k] = -g_zz_xz[k] * ab_y + g_zz_xyz[k];

                g_yzz_yy[k] = -g_zz_yy[k] * ab_y + g_zz_yyy[k];

                g_yzz_yz[k] = -g_zz_yz[k] * ab_y + g_zz_yyz[k];

                g_yzz_zz[k] = -g_zz_zz[k] * ab_y + g_zz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_zzz_xx = cbuffer.data(fd_off + 54 * ccomps * dcomps);

            auto g_zzz_xy = cbuffer.data(fd_off + 55 * ccomps * dcomps);

            auto g_zzz_xz = cbuffer.data(fd_off + 56 * ccomps * dcomps);

            auto g_zzz_yy = cbuffer.data(fd_off + 57 * ccomps * dcomps);

            auto g_zzz_yz = cbuffer.data(fd_off + 58 * ccomps * dcomps);

            auto g_zzz_zz = cbuffer.data(fd_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_xx, g_zz_xxz, g_zz_xy, g_zz_xyz, g_zz_xz, g_zz_xzz, g_zz_yy, g_zz_yyz, g_zz_yz, g_zz_yzz, g_zz_zz, g_zz_zzz, g_zzz_xx, g_zzz_xy, g_zzz_xz, g_zzz_yy, g_zzz_yz, g_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zzz_xx[k] = -g_zz_xx[k] * ab_z + g_zz_xxz[k];

                g_zzz_xy[k] = -g_zz_xy[k] * ab_z + g_zz_xyz[k];

                g_zzz_xz[k] = -g_zz_xz[k] * ab_z + g_zz_xzz[k];

                g_zzz_yy[k] = -g_zz_yy[k] * ab_z + g_zz_yyz[k];

                g_zzz_yz[k] = -g_zz_yz[k] * ab_z + g_zz_yzz[k];

                g_zzz_zz[k] = -g_zz_zz[k] * ab_z + g_zz_zzz[k];
            }
        }
    }
}

} // erirec namespace

