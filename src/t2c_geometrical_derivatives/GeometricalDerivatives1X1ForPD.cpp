#include "GeometricalDerivatives1X1ForPD.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_pd(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_pd,
                        const size_t              idx_op_sp,
                        const size_t              idx_op_sf,
                        const size_t              idx_op_dp,
                        const size_t              idx_op_df,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SP

        auto to_0_x = pbuffer.data(idx_op_sp + i * 3 + 0);

        auto to_0_y = pbuffer.data(idx_op_sp + i * 3 + 1);

        auto to_0_z = pbuffer.data(idx_op_sp + i * 3 + 2);

        // Set up components of auxiliary buffer : SF

        auto to_0_xxx = pbuffer.data(idx_op_sf + i * 10 + 0);

        auto to_0_xxy = pbuffer.data(idx_op_sf + i * 10 + 1);

        auto to_0_xxz = pbuffer.data(idx_op_sf + i * 10 + 2);

        auto to_0_xyy = pbuffer.data(idx_op_sf + i * 10 + 3);

        auto to_0_xyz = pbuffer.data(idx_op_sf + i * 10 + 4);

        auto to_0_xzz = pbuffer.data(idx_op_sf + i * 10 + 5);

        auto to_0_yyy = pbuffer.data(idx_op_sf + i * 10 + 6);

        auto to_0_yyz = pbuffer.data(idx_op_sf + i * 10 + 7);

        auto to_0_yzz = pbuffer.data(idx_op_sf + i * 10 + 8);

        auto to_0_zzz = pbuffer.data(idx_op_sf + i * 10 + 9);

        // Set up components of auxiliary buffer : DP

        auto to_xx_x = pbuffer.data(idx_op_dp + i * 18 + 0);

        auto to_xx_y = pbuffer.data(idx_op_dp + i * 18 + 1);

        auto to_xx_z = pbuffer.data(idx_op_dp + i * 18 + 2);

        auto to_xy_x = pbuffer.data(idx_op_dp + i * 18 + 3);

        auto to_xy_y = pbuffer.data(idx_op_dp + i * 18 + 4);

        auto to_xy_z = pbuffer.data(idx_op_dp + i * 18 + 5);

        auto to_xz_x = pbuffer.data(idx_op_dp + i * 18 + 6);

        auto to_xz_y = pbuffer.data(idx_op_dp + i * 18 + 7);

        auto to_xz_z = pbuffer.data(idx_op_dp + i * 18 + 8);

        auto to_yy_x = pbuffer.data(idx_op_dp + i * 18 + 9);

        auto to_yy_y = pbuffer.data(idx_op_dp + i * 18 + 10);

        auto to_yy_z = pbuffer.data(idx_op_dp + i * 18 + 11);

        auto to_yz_x = pbuffer.data(idx_op_dp + i * 18 + 12);

        auto to_yz_y = pbuffer.data(idx_op_dp + i * 18 + 13);

        auto to_yz_z = pbuffer.data(idx_op_dp + i * 18 + 14);

        auto to_zz_x = pbuffer.data(idx_op_dp + i * 18 + 15);

        auto to_zz_y = pbuffer.data(idx_op_dp + i * 18 + 16);

        auto to_zz_z = pbuffer.data(idx_op_dp + i * 18 + 17);

        // Set up components of auxiliary buffer : DF

        auto to_xx_xxx = pbuffer.data(idx_op_df + i * 60 + 0);

        auto to_xx_xxy = pbuffer.data(idx_op_df + i * 60 + 1);

        auto to_xx_xxz = pbuffer.data(idx_op_df + i * 60 + 2);

        auto to_xx_xyy = pbuffer.data(idx_op_df + i * 60 + 3);

        auto to_xx_xyz = pbuffer.data(idx_op_df + i * 60 + 4);

        auto to_xx_xzz = pbuffer.data(idx_op_df + i * 60 + 5);

        auto to_xx_yyy = pbuffer.data(idx_op_df + i * 60 + 6);

        auto to_xx_yyz = pbuffer.data(idx_op_df + i * 60 + 7);

        auto to_xx_yzz = pbuffer.data(idx_op_df + i * 60 + 8);

        auto to_xx_zzz = pbuffer.data(idx_op_df + i * 60 + 9);

        auto to_xy_xxx = pbuffer.data(idx_op_df + i * 60 + 10);

        auto to_xy_xxy = pbuffer.data(idx_op_df + i * 60 + 11);

        auto to_xy_xxz = pbuffer.data(idx_op_df + i * 60 + 12);

        auto to_xy_xyy = pbuffer.data(idx_op_df + i * 60 + 13);

        auto to_xy_xyz = pbuffer.data(idx_op_df + i * 60 + 14);

        auto to_xy_xzz = pbuffer.data(idx_op_df + i * 60 + 15);

        auto to_xy_yyy = pbuffer.data(idx_op_df + i * 60 + 16);

        auto to_xy_yyz = pbuffer.data(idx_op_df + i * 60 + 17);

        auto to_xy_yzz = pbuffer.data(idx_op_df + i * 60 + 18);

        auto to_xy_zzz = pbuffer.data(idx_op_df + i * 60 + 19);

        auto to_xz_xxx = pbuffer.data(idx_op_df + i * 60 + 20);

        auto to_xz_xxy = pbuffer.data(idx_op_df + i * 60 + 21);

        auto to_xz_xxz = pbuffer.data(idx_op_df + i * 60 + 22);

        auto to_xz_xyy = pbuffer.data(idx_op_df + i * 60 + 23);

        auto to_xz_xyz = pbuffer.data(idx_op_df + i * 60 + 24);

        auto to_xz_xzz = pbuffer.data(idx_op_df + i * 60 + 25);

        auto to_xz_yyy = pbuffer.data(idx_op_df + i * 60 + 26);

        auto to_xz_yyz = pbuffer.data(idx_op_df + i * 60 + 27);

        auto to_xz_yzz = pbuffer.data(idx_op_df + i * 60 + 28);

        auto to_xz_zzz = pbuffer.data(idx_op_df + i * 60 + 29);

        auto to_yy_xxx = pbuffer.data(idx_op_df + i * 60 + 30);

        auto to_yy_xxy = pbuffer.data(idx_op_df + i * 60 + 31);

        auto to_yy_xxz = pbuffer.data(idx_op_df + i * 60 + 32);

        auto to_yy_xyy = pbuffer.data(idx_op_df + i * 60 + 33);

        auto to_yy_xyz = pbuffer.data(idx_op_df + i * 60 + 34);

        auto to_yy_xzz = pbuffer.data(idx_op_df + i * 60 + 35);

        auto to_yy_yyy = pbuffer.data(idx_op_df + i * 60 + 36);

        auto to_yy_yyz = pbuffer.data(idx_op_df + i * 60 + 37);

        auto to_yy_yzz = pbuffer.data(idx_op_df + i * 60 + 38);

        auto to_yy_zzz = pbuffer.data(idx_op_df + i * 60 + 39);

        auto to_yz_xxx = pbuffer.data(idx_op_df + i * 60 + 40);

        auto to_yz_xxy = pbuffer.data(idx_op_df + i * 60 + 41);

        auto to_yz_xxz = pbuffer.data(idx_op_df + i * 60 + 42);

        auto to_yz_xyy = pbuffer.data(idx_op_df + i * 60 + 43);

        auto to_yz_xyz = pbuffer.data(idx_op_df + i * 60 + 44);

        auto to_yz_xzz = pbuffer.data(idx_op_df + i * 60 + 45);

        auto to_yz_yyy = pbuffer.data(idx_op_df + i * 60 + 46);

        auto to_yz_yyz = pbuffer.data(idx_op_df + i * 60 + 47);

        auto to_yz_yzz = pbuffer.data(idx_op_df + i * 60 + 48);

        auto to_yz_zzz = pbuffer.data(idx_op_df + i * 60 + 49);

        auto to_zz_xxx = pbuffer.data(idx_op_df + i * 60 + 50);

        auto to_zz_xxy = pbuffer.data(idx_op_df + i * 60 + 51);

        auto to_zz_xxz = pbuffer.data(idx_op_df + i * 60 + 52);

        auto to_zz_xyy = pbuffer.data(idx_op_df + i * 60 + 53);

        auto to_zz_xyz = pbuffer.data(idx_op_df + i * 60 + 54);

        auto to_zz_xzz = pbuffer.data(idx_op_df + i * 60 + 55);

        auto to_zz_yyy = pbuffer.data(idx_op_df + i * 60 + 56);

        auto to_zz_yyz = pbuffer.data(idx_op_df + i * 60 + 57);

        auto to_zz_yzz = pbuffer.data(idx_op_df + i * 60 + 58);

        auto to_zz_zzz = pbuffer.data(idx_op_df + i * 60 + 59);

        // Set up 0-6 components of targeted buffer : PD

        auto to_x_x_x_xx = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 0);

        auto to_x_x_x_xy = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 1);

        auto to_x_x_x_xz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 2);

        auto to_x_x_x_yy = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 3);

        auto to_x_x_x_yz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 4);

        auto to_x_x_x_zz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxx,    \
                             to_0_xxy,    \
                             to_0_xxz,    \
                             to_0_xyy,    \
                             to_0_xyz,    \
                             to_0_xzz,    \
                             to_0_y,      \
                             to_0_z,      \
                             to_x_x_x_xx, \
                             to_x_x_x_xy, \
                             to_x_x_x_xz, \
                             to_x_x_x_yy, \
                             to_x_x_x_yz, \
                             to_x_x_x_zz, \
                             to_xx_x,     \
                             to_xx_xxx,   \
                             to_xx_xxy,   \
                             to_xx_xxz,   \
                             to_xx_xyy,   \
                             to_xx_xyz,   \
                             to_xx_xzz,   \
                             to_xx_y,     \
                             to_xx_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_x_xx[k] = 2.0 * to_0_x[k] - 2.0 * to_0_xxx[k] * tke_0 - 4.0 * to_xx_x[k] * tbe_0 + 4.0 * to_xx_xxx[k] * tbe_0 * tke_0;

            to_x_x_x_xy[k] = to_0_y[k] - 2.0 * to_0_xxy[k] * tke_0 - 2.0 * to_xx_y[k] * tbe_0 + 4.0 * to_xx_xxy[k] * tbe_0 * tke_0;

            to_x_x_x_xz[k] = to_0_z[k] - 2.0 * to_0_xxz[k] * tke_0 - 2.0 * to_xx_z[k] * tbe_0 + 4.0 * to_xx_xxz[k] * tbe_0 * tke_0;

            to_x_x_x_yy[k] = -2.0 * to_0_xyy[k] * tke_0 + 4.0 * to_xx_xyy[k] * tbe_0 * tke_0;

            to_x_x_x_yz[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_xx_xyz[k] * tbe_0 * tke_0;

            to_x_x_x_zz[k] = -2.0 * to_0_xzz[k] * tke_0 + 4.0 * to_xx_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 6-12 components of targeted buffer : PD

        auto to_x_x_y_xx = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 6);

        auto to_x_x_y_xy = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 7);

        auto to_x_x_y_xz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 8);

        auto to_x_x_y_yy = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 9);

        auto to_x_x_y_yz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 10);

        auto to_x_x_y_zz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_x_x_y_xx,     \
                             to_x_x_y_xy, \
                             to_x_x_y_xz, \
                             to_x_x_y_yy, \
                             to_x_x_y_yz, \
                             to_x_x_y_zz, \
                             to_xy_x,     \
                             to_xy_xxx,   \
                             to_xy_xxy,   \
                             to_xy_xxz,   \
                             to_xy_xyy,   \
                             to_xy_xyz,   \
                             to_xy_xzz,   \
                             to_xy_y,     \
                             to_xy_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_y_xx[k] = -4.0 * to_xy_x[k] * tbe_0 + 4.0 * to_xy_xxx[k] * tbe_0 * tke_0;

            to_x_x_y_xy[k] = -2.0 * to_xy_y[k] * tbe_0 + 4.0 * to_xy_xxy[k] * tbe_0 * tke_0;

            to_x_x_y_xz[k] = -2.0 * to_xy_z[k] * tbe_0 + 4.0 * to_xy_xxz[k] * tbe_0 * tke_0;

            to_x_x_y_yy[k] = 4.0 * to_xy_xyy[k] * tbe_0 * tke_0;

            to_x_x_y_yz[k] = 4.0 * to_xy_xyz[k] * tbe_0 * tke_0;

            to_x_x_y_zz[k] = 4.0 * to_xy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 12-18 components of targeted buffer : PD

        auto to_x_x_z_xx = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 12);

        auto to_x_x_z_xy = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 13);

        auto to_x_x_z_xz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 14);

        auto to_x_x_z_yy = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 15);

        auto to_x_x_z_yz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 16);

        auto to_x_x_z_zz = pbuffer.data(idx_op_geom_101_pd + 0 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_x_x_z_xx,     \
                             to_x_x_z_xy, \
                             to_x_x_z_xz, \
                             to_x_x_z_yy, \
                             to_x_x_z_yz, \
                             to_x_x_z_zz, \
                             to_xz_x,     \
                             to_xz_xxx,   \
                             to_xz_xxy,   \
                             to_xz_xxz,   \
                             to_xz_xyy,   \
                             to_xz_xyz,   \
                             to_xz_xzz,   \
                             to_xz_y,     \
                             to_xz_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_z_xx[k] = -4.0 * to_xz_x[k] * tbe_0 + 4.0 * to_xz_xxx[k] * tbe_0 * tke_0;

            to_x_x_z_xy[k] = -2.0 * to_xz_y[k] * tbe_0 + 4.0 * to_xz_xxy[k] * tbe_0 * tke_0;

            to_x_x_z_xz[k] = -2.0 * to_xz_z[k] * tbe_0 + 4.0 * to_xz_xxz[k] * tbe_0 * tke_0;

            to_x_x_z_yy[k] = 4.0 * to_xz_xyy[k] * tbe_0 * tke_0;

            to_x_x_z_yz[k] = 4.0 * to_xz_xyz[k] * tbe_0 * tke_0;

            to_x_x_z_zz[k] = 4.0 * to_xz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 18-24 components of targeted buffer : PD

        auto to_x_y_x_xx = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 0);

        auto to_x_y_x_xy = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 1);

        auto to_x_y_x_xz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 2);

        auto to_x_y_x_yy = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 3);

        auto to_x_y_x_yz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 4);

        auto to_x_y_x_zz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxy,    \
                             to_0_xyy,    \
                             to_0_xyz,    \
                             to_0_y,      \
                             to_0_yyy,    \
                             to_0_yyz,    \
                             to_0_yzz,    \
                             to_0_z,      \
                             to_x_y_x_xx, \
                             to_x_y_x_xy, \
                             to_x_y_x_xz, \
                             to_x_y_x_yy, \
                             to_x_y_x_yz, \
                             to_x_y_x_zz, \
                             to_xx_x,     \
                             to_xx_xxy,   \
                             to_xx_xyy,   \
                             to_xx_xyz,   \
                             to_xx_y,     \
                             to_xx_yyy,   \
                             to_xx_yyz,   \
                             to_xx_yzz,   \
                             to_xx_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_x_xx[k] = -2.0 * to_0_xxy[k] * tke_0 + 4.0 * to_xx_xxy[k] * tbe_0 * tke_0;

            to_x_y_x_xy[k] = to_0_x[k] - 2.0 * to_0_xyy[k] * tke_0 - 2.0 * to_xx_x[k] * tbe_0 + 4.0 * to_xx_xyy[k] * tbe_0 * tke_0;

            to_x_y_x_xz[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_xx_xyz[k] * tbe_0 * tke_0;

            to_x_y_x_yy[k] = 2.0 * to_0_y[k] - 2.0 * to_0_yyy[k] * tke_0 - 4.0 * to_xx_y[k] * tbe_0 + 4.0 * to_xx_yyy[k] * tbe_0 * tke_0;

            to_x_y_x_yz[k] = to_0_z[k] - 2.0 * to_0_yyz[k] * tke_0 - 2.0 * to_xx_z[k] * tbe_0 + 4.0 * to_xx_yyz[k] * tbe_0 * tke_0;

            to_x_y_x_zz[k] = -2.0 * to_0_yzz[k] * tke_0 + 4.0 * to_xx_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 24-30 components of targeted buffer : PD

        auto to_x_y_y_xx = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 6);

        auto to_x_y_y_xy = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 7);

        auto to_x_y_y_xz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 8);

        auto to_x_y_y_yy = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 9);

        auto to_x_y_y_yz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 10);

        auto to_x_y_y_zz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_x_y_y_xx,     \
                             to_x_y_y_xy, \
                             to_x_y_y_xz, \
                             to_x_y_y_yy, \
                             to_x_y_y_yz, \
                             to_x_y_y_zz, \
                             to_xy_x,     \
                             to_xy_xxy,   \
                             to_xy_xyy,   \
                             to_xy_xyz,   \
                             to_xy_y,     \
                             to_xy_yyy,   \
                             to_xy_yyz,   \
                             to_xy_yzz,   \
                             to_xy_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_y_xx[k] = 4.0 * to_xy_xxy[k] * tbe_0 * tke_0;

            to_x_y_y_xy[k] = -2.0 * to_xy_x[k] * tbe_0 + 4.0 * to_xy_xyy[k] * tbe_0 * tke_0;

            to_x_y_y_xz[k] = 4.0 * to_xy_xyz[k] * tbe_0 * tke_0;

            to_x_y_y_yy[k] = -4.0 * to_xy_y[k] * tbe_0 + 4.0 * to_xy_yyy[k] * tbe_0 * tke_0;

            to_x_y_y_yz[k] = -2.0 * to_xy_z[k] * tbe_0 + 4.0 * to_xy_yyz[k] * tbe_0 * tke_0;

            to_x_y_y_zz[k] = 4.0 * to_xy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-36 components of targeted buffer : PD

        auto to_x_y_z_xx = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 12);

        auto to_x_y_z_xy = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 13);

        auto to_x_y_z_xz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 14);

        auto to_x_y_z_yy = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 15);

        auto to_x_y_z_yz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 16);

        auto to_x_y_z_zz = pbuffer.data(idx_op_geom_101_pd + 1 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_x_y_z_xx,     \
                             to_x_y_z_xy, \
                             to_x_y_z_xz, \
                             to_x_y_z_yy, \
                             to_x_y_z_yz, \
                             to_x_y_z_zz, \
                             to_xz_x,     \
                             to_xz_xxy,   \
                             to_xz_xyy,   \
                             to_xz_xyz,   \
                             to_xz_y,     \
                             to_xz_yyy,   \
                             to_xz_yyz,   \
                             to_xz_yzz,   \
                             to_xz_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_z_xx[k] = 4.0 * to_xz_xxy[k] * tbe_0 * tke_0;

            to_x_y_z_xy[k] = -2.0 * to_xz_x[k] * tbe_0 + 4.0 * to_xz_xyy[k] * tbe_0 * tke_0;

            to_x_y_z_xz[k] = 4.0 * to_xz_xyz[k] * tbe_0 * tke_0;

            to_x_y_z_yy[k] = -4.0 * to_xz_y[k] * tbe_0 + 4.0 * to_xz_yyy[k] * tbe_0 * tke_0;

            to_x_y_z_yz[k] = -2.0 * to_xz_z[k] * tbe_0 + 4.0 * to_xz_yyz[k] * tbe_0 * tke_0;

            to_x_y_z_zz[k] = 4.0 * to_xz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 36-42 components of targeted buffer : PD

        auto to_x_z_x_xx = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 0);

        auto to_x_z_x_xy = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 1);

        auto to_x_z_x_xz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 2);

        auto to_x_z_x_yy = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 3);

        auto to_x_z_x_yz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 4);

        auto to_x_z_x_zz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxz,    \
                             to_0_xyz,    \
                             to_0_xzz,    \
                             to_0_y,      \
                             to_0_yyz,    \
                             to_0_yzz,    \
                             to_0_z,      \
                             to_0_zzz,    \
                             to_x_z_x_xx, \
                             to_x_z_x_xy, \
                             to_x_z_x_xz, \
                             to_x_z_x_yy, \
                             to_x_z_x_yz, \
                             to_x_z_x_zz, \
                             to_xx_x,     \
                             to_xx_xxz,   \
                             to_xx_xyz,   \
                             to_xx_xzz,   \
                             to_xx_y,     \
                             to_xx_yyz,   \
                             to_xx_yzz,   \
                             to_xx_z,     \
                             to_xx_zzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_x_xx[k] = -2.0 * to_0_xxz[k] * tke_0 + 4.0 * to_xx_xxz[k] * tbe_0 * tke_0;

            to_x_z_x_xy[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_xx_xyz[k] * tbe_0 * tke_0;

            to_x_z_x_xz[k] = to_0_x[k] - 2.0 * to_0_xzz[k] * tke_0 - 2.0 * to_xx_x[k] * tbe_0 + 4.0 * to_xx_xzz[k] * tbe_0 * tke_0;

            to_x_z_x_yy[k] = -2.0 * to_0_yyz[k] * tke_0 + 4.0 * to_xx_yyz[k] * tbe_0 * tke_0;

            to_x_z_x_yz[k] = to_0_y[k] - 2.0 * to_0_yzz[k] * tke_0 - 2.0 * to_xx_y[k] * tbe_0 + 4.0 * to_xx_yzz[k] * tbe_0 * tke_0;

            to_x_z_x_zz[k] = 2.0 * to_0_z[k] - 2.0 * to_0_zzz[k] * tke_0 - 4.0 * to_xx_z[k] * tbe_0 + 4.0 * to_xx_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 42-48 components of targeted buffer : PD

        auto to_x_z_y_xx = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 6);

        auto to_x_z_y_xy = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 7);

        auto to_x_z_y_xz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 8);

        auto to_x_z_y_yy = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 9);

        auto to_x_z_y_yz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 10);

        auto to_x_z_y_zz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_x_z_y_xx,     \
                             to_x_z_y_xy, \
                             to_x_z_y_xz, \
                             to_x_z_y_yy, \
                             to_x_z_y_yz, \
                             to_x_z_y_zz, \
                             to_xy_x,     \
                             to_xy_xxz,   \
                             to_xy_xyz,   \
                             to_xy_xzz,   \
                             to_xy_y,     \
                             to_xy_yyz,   \
                             to_xy_yzz,   \
                             to_xy_z,     \
                             to_xy_zzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_y_xx[k] = 4.0 * to_xy_xxz[k] * tbe_0 * tke_0;

            to_x_z_y_xy[k] = 4.0 * to_xy_xyz[k] * tbe_0 * tke_0;

            to_x_z_y_xz[k] = -2.0 * to_xy_x[k] * tbe_0 + 4.0 * to_xy_xzz[k] * tbe_0 * tke_0;

            to_x_z_y_yy[k] = 4.0 * to_xy_yyz[k] * tbe_0 * tke_0;

            to_x_z_y_yz[k] = -2.0 * to_xy_y[k] * tbe_0 + 4.0 * to_xy_yzz[k] * tbe_0 * tke_0;

            to_x_z_y_zz[k] = -4.0 * to_xy_z[k] * tbe_0 + 4.0 * to_xy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 48-54 components of targeted buffer : PD

        auto to_x_z_z_xx = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 12);

        auto to_x_z_z_xy = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 13);

        auto to_x_z_z_xz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 14);

        auto to_x_z_z_yy = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 15);

        auto to_x_z_z_yz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 16);

        auto to_x_z_z_zz = pbuffer.data(idx_op_geom_101_pd + 2 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_x_z_z_xx,     \
                             to_x_z_z_xy, \
                             to_x_z_z_xz, \
                             to_x_z_z_yy, \
                             to_x_z_z_yz, \
                             to_x_z_z_zz, \
                             to_xz_x,     \
                             to_xz_xxz,   \
                             to_xz_xyz,   \
                             to_xz_xzz,   \
                             to_xz_y,     \
                             to_xz_yyz,   \
                             to_xz_yzz,   \
                             to_xz_z,     \
                             to_xz_zzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_z_xx[k] = 4.0 * to_xz_xxz[k] * tbe_0 * tke_0;

            to_x_z_z_xy[k] = 4.0 * to_xz_xyz[k] * tbe_0 * tke_0;

            to_x_z_z_xz[k] = -2.0 * to_xz_x[k] * tbe_0 + 4.0 * to_xz_xzz[k] * tbe_0 * tke_0;

            to_x_z_z_yy[k] = 4.0 * to_xz_yyz[k] * tbe_0 * tke_0;

            to_x_z_z_yz[k] = -2.0 * to_xz_y[k] * tbe_0 + 4.0 * to_xz_yzz[k] * tbe_0 * tke_0;

            to_x_z_z_zz[k] = -4.0 * to_xz_z[k] * tbe_0 + 4.0 * to_xz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 54-60 components of targeted buffer : PD

        auto to_y_x_x_xx = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 0);

        auto to_y_x_x_xy = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 1);

        auto to_y_x_x_xz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 2);

        auto to_y_x_x_yy = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 3);

        auto to_y_x_x_yz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 4);

        auto to_y_x_x_zz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_xy_x,         \
                             to_xy_xxx,   \
                             to_xy_xxy,   \
                             to_xy_xxz,   \
                             to_xy_xyy,   \
                             to_xy_xyz,   \
                             to_xy_xzz,   \
                             to_xy_y,     \
                             to_xy_z,     \
                             to_y_x_x_xx, \
                             to_y_x_x_xy, \
                             to_y_x_x_xz, \
                             to_y_x_x_yy, \
                             to_y_x_x_yz, \
                             to_y_x_x_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_x_xx[k] = -4.0 * to_xy_x[k] * tbe_0 + 4.0 * to_xy_xxx[k] * tbe_0 * tke_0;

            to_y_x_x_xy[k] = -2.0 * to_xy_y[k] * tbe_0 + 4.0 * to_xy_xxy[k] * tbe_0 * tke_0;

            to_y_x_x_xz[k] = -2.0 * to_xy_z[k] * tbe_0 + 4.0 * to_xy_xxz[k] * tbe_0 * tke_0;

            to_y_x_x_yy[k] = 4.0 * to_xy_xyy[k] * tbe_0 * tke_0;

            to_y_x_x_yz[k] = 4.0 * to_xy_xyz[k] * tbe_0 * tke_0;

            to_y_x_x_zz[k] = 4.0 * to_xy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-66 components of targeted buffer : PD

        auto to_y_x_y_xx = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 6);

        auto to_y_x_y_xy = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 7);

        auto to_y_x_y_xz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 8);

        auto to_y_x_y_yy = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 9);

        auto to_y_x_y_yz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 10);

        auto to_y_x_y_zz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxx,    \
                             to_0_xxy,    \
                             to_0_xxz,    \
                             to_0_xyy,    \
                             to_0_xyz,    \
                             to_0_xzz,    \
                             to_0_y,      \
                             to_0_z,      \
                             to_y_x_y_xx, \
                             to_y_x_y_xy, \
                             to_y_x_y_xz, \
                             to_y_x_y_yy, \
                             to_y_x_y_yz, \
                             to_y_x_y_zz, \
                             to_yy_x,     \
                             to_yy_xxx,   \
                             to_yy_xxy,   \
                             to_yy_xxz,   \
                             to_yy_xyy,   \
                             to_yy_xyz,   \
                             to_yy_xzz,   \
                             to_yy_y,     \
                             to_yy_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_y_xx[k] = 2.0 * to_0_x[k] - 2.0 * to_0_xxx[k] * tke_0 - 4.0 * to_yy_x[k] * tbe_0 + 4.0 * to_yy_xxx[k] * tbe_0 * tke_0;

            to_y_x_y_xy[k] = to_0_y[k] - 2.0 * to_0_xxy[k] * tke_0 - 2.0 * to_yy_y[k] * tbe_0 + 4.0 * to_yy_xxy[k] * tbe_0 * tke_0;

            to_y_x_y_xz[k] = to_0_z[k] - 2.0 * to_0_xxz[k] * tke_0 - 2.0 * to_yy_z[k] * tbe_0 + 4.0 * to_yy_xxz[k] * tbe_0 * tke_0;

            to_y_x_y_yy[k] = -2.0 * to_0_xyy[k] * tke_0 + 4.0 * to_yy_xyy[k] * tbe_0 * tke_0;

            to_y_x_y_yz[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_yy_xyz[k] * tbe_0 * tke_0;

            to_y_x_y_zz[k] = -2.0 * to_0_xzz[k] * tke_0 + 4.0 * to_yy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 66-72 components of targeted buffer : PD

        auto to_y_x_z_xx = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 12);

        auto to_y_x_z_xy = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 13);

        auto to_y_x_z_xz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 14);

        auto to_y_x_z_yy = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 15);

        auto to_y_x_z_yz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 16);

        auto to_y_x_z_zz = pbuffer.data(idx_op_geom_101_pd + 3 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_y_x_z_xx,     \
                             to_y_x_z_xy, \
                             to_y_x_z_xz, \
                             to_y_x_z_yy, \
                             to_y_x_z_yz, \
                             to_y_x_z_zz, \
                             to_yz_x,     \
                             to_yz_xxx,   \
                             to_yz_xxy,   \
                             to_yz_xxz,   \
                             to_yz_xyy,   \
                             to_yz_xyz,   \
                             to_yz_xzz,   \
                             to_yz_y,     \
                             to_yz_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_z_xx[k] = -4.0 * to_yz_x[k] * tbe_0 + 4.0 * to_yz_xxx[k] * tbe_0 * tke_0;

            to_y_x_z_xy[k] = -2.0 * to_yz_y[k] * tbe_0 + 4.0 * to_yz_xxy[k] * tbe_0 * tke_0;

            to_y_x_z_xz[k] = -2.0 * to_yz_z[k] * tbe_0 + 4.0 * to_yz_xxz[k] * tbe_0 * tke_0;

            to_y_x_z_yy[k] = 4.0 * to_yz_xyy[k] * tbe_0 * tke_0;

            to_y_x_z_yz[k] = 4.0 * to_yz_xyz[k] * tbe_0 * tke_0;

            to_y_x_z_zz[k] = 4.0 * to_yz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 72-78 components of targeted buffer : PD

        auto to_y_y_x_xx = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 0);

        auto to_y_y_x_xy = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 1);

        auto to_y_y_x_xz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 2);

        auto to_y_y_x_yy = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 3);

        auto to_y_y_x_yz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 4);

        auto to_y_y_x_zz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_xy_x,         \
                             to_xy_xxy,   \
                             to_xy_xyy,   \
                             to_xy_xyz,   \
                             to_xy_y,     \
                             to_xy_yyy,   \
                             to_xy_yyz,   \
                             to_xy_yzz,   \
                             to_xy_z,     \
                             to_y_y_x_xx, \
                             to_y_y_x_xy, \
                             to_y_y_x_xz, \
                             to_y_y_x_yy, \
                             to_y_y_x_yz, \
                             to_y_y_x_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_x_xx[k] = 4.0 * to_xy_xxy[k] * tbe_0 * tke_0;

            to_y_y_x_xy[k] = -2.0 * to_xy_x[k] * tbe_0 + 4.0 * to_xy_xyy[k] * tbe_0 * tke_0;

            to_y_y_x_xz[k] = 4.0 * to_xy_xyz[k] * tbe_0 * tke_0;

            to_y_y_x_yy[k] = -4.0 * to_xy_y[k] * tbe_0 + 4.0 * to_xy_yyy[k] * tbe_0 * tke_0;

            to_y_y_x_yz[k] = -2.0 * to_xy_z[k] * tbe_0 + 4.0 * to_xy_yyz[k] * tbe_0 * tke_0;

            to_y_y_x_zz[k] = 4.0 * to_xy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 78-84 components of targeted buffer : PD

        auto to_y_y_y_xx = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 6);

        auto to_y_y_y_xy = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 7);

        auto to_y_y_y_xz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 8);

        auto to_y_y_y_yy = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 9);

        auto to_y_y_y_yz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 10);

        auto to_y_y_y_zz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxy,    \
                             to_0_xyy,    \
                             to_0_xyz,    \
                             to_0_y,      \
                             to_0_yyy,    \
                             to_0_yyz,    \
                             to_0_yzz,    \
                             to_0_z,      \
                             to_y_y_y_xx, \
                             to_y_y_y_xy, \
                             to_y_y_y_xz, \
                             to_y_y_y_yy, \
                             to_y_y_y_yz, \
                             to_y_y_y_zz, \
                             to_yy_x,     \
                             to_yy_xxy,   \
                             to_yy_xyy,   \
                             to_yy_xyz,   \
                             to_yy_y,     \
                             to_yy_yyy,   \
                             to_yy_yyz,   \
                             to_yy_yzz,   \
                             to_yy_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_y_xx[k] = -2.0 * to_0_xxy[k] * tke_0 + 4.0 * to_yy_xxy[k] * tbe_0 * tke_0;

            to_y_y_y_xy[k] = to_0_x[k] - 2.0 * to_0_xyy[k] * tke_0 - 2.0 * to_yy_x[k] * tbe_0 + 4.0 * to_yy_xyy[k] * tbe_0 * tke_0;

            to_y_y_y_xz[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_yy_xyz[k] * tbe_0 * tke_0;

            to_y_y_y_yy[k] = 2.0 * to_0_y[k] - 2.0 * to_0_yyy[k] * tke_0 - 4.0 * to_yy_y[k] * tbe_0 + 4.0 * to_yy_yyy[k] * tbe_0 * tke_0;

            to_y_y_y_yz[k] = to_0_z[k] - 2.0 * to_0_yyz[k] * tke_0 - 2.0 * to_yy_z[k] * tbe_0 + 4.0 * to_yy_yyz[k] * tbe_0 * tke_0;

            to_y_y_y_zz[k] = -2.0 * to_0_yzz[k] * tke_0 + 4.0 * to_yy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 84-90 components of targeted buffer : PD

        auto to_y_y_z_xx = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 12);

        auto to_y_y_z_xy = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 13);

        auto to_y_y_z_xz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 14);

        auto to_y_y_z_yy = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 15);

        auto to_y_y_z_yz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 16);

        auto to_y_y_z_zz = pbuffer.data(idx_op_geom_101_pd + 4 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_y_y_z_xx,     \
                             to_y_y_z_xy, \
                             to_y_y_z_xz, \
                             to_y_y_z_yy, \
                             to_y_y_z_yz, \
                             to_y_y_z_zz, \
                             to_yz_x,     \
                             to_yz_xxy,   \
                             to_yz_xyy,   \
                             to_yz_xyz,   \
                             to_yz_y,     \
                             to_yz_yyy,   \
                             to_yz_yyz,   \
                             to_yz_yzz,   \
                             to_yz_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_z_xx[k] = 4.0 * to_yz_xxy[k] * tbe_0 * tke_0;

            to_y_y_z_xy[k] = -2.0 * to_yz_x[k] * tbe_0 + 4.0 * to_yz_xyy[k] * tbe_0 * tke_0;

            to_y_y_z_xz[k] = 4.0 * to_yz_xyz[k] * tbe_0 * tke_0;

            to_y_y_z_yy[k] = -4.0 * to_yz_y[k] * tbe_0 + 4.0 * to_yz_yyy[k] * tbe_0 * tke_0;

            to_y_y_z_yz[k] = -2.0 * to_yz_z[k] * tbe_0 + 4.0 * to_yz_yyz[k] * tbe_0 * tke_0;

            to_y_y_z_zz[k] = 4.0 * to_yz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-96 components of targeted buffer : PD

        auto to_y_z_x_xx = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 0);

        auto to_y_z_x_xy = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 1);

        auto to_y_z_x_xz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 2);

        auto to_y_z_x_yy = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 3);

        auto to_y_z_x_yz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 4);

        auto to_y_z_x_zz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_xy_x,         \
                             to_xy_xxz,   \
                             to_xy_xyz,   \
                             to_xy_xzz,   \
                             to_xy_y,     \
                             to_xy_yyz,   \
                             to_xy_yzz,   \
                             to_xy_z,     \
                             to_xy_zzz,   \
                             to_y_z_x_xx, \
                             to_y_z_x_xy, \
                             to_y_z_x_xz, \
                             to_y_z_x_yy, \
                             to_y_z_x_yz, \
                             to_y_z_x_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_x_xx[k] = 4.0 * to_xy_xxz[k] * tbe_0 * tke_0;

            to_y_z_x_xy[k] = 4.0 * to_xy_xyz[k] * tbe_0 * tke_0;

            to_y_z_x_xz[k] = -2.0 * to_xy_x[k] * tbe_0 + 4.0 * to_xy_xzz[k] * tbe_0 * tke_0;

            to_y_z_x_yy[k] = 4.0 * to_xy_yyz[k] * tbe_0 * tke_0;

            to_y_z_x_yz[k] = -2.0 * to_xy_y[k] * tbe_0 + 4.0 * to_xy_yzz[k] * tbe_0 * tke_0;

            to_y_z_x_zz[k] = -4.0 * to_xy_z[k] * tbe_0 + 4.0 * to_xy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 96-102 components of targeted buffer : PD

        auto to_y_z_y_xx = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 6);

        auto to_y_z_y_xy = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 7);

        auto to_y_z_y_xz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 8);

        auto to_y_z_y_yy = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 9);

        auto to_y_z_y_yz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 10);

        auto to_y_z_y_zz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxz,    \
                             to_0_xyz,    \
                             to_0_xzz,    \
                             to_0_y,      \
                             to_0_yyz,    \
                             to_0_yzz,    \
                             to_0_z,      \
                             to_0_zzz,    \
                             to_y_z_y_xx, \
                             to_y_z_y_xy, \
                             to_y_z_y_xz, \
                             to_y_z_y_yy, \
                             to_y_z_y_yz, \
                             to_y_z_y_zz, \
                             to_yy_x,     \
                             to_yy_xxz,   \
                             to_yy_xyz,   \
                             to_yy_xzz,   \
                             to_yy_y,     \
                             to_yy_yyz,   \
                             to_yy_yzz,   \
                             to_yy_z,     \
                             to_yy_zzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_y_xx[k] = -2.0 * to_0_xxz[k] * tke_0 + 4.0 * to_yy_xxz[k] * tbe_0 * tke_0;

            to_y_z_y_xy[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_yy_xyz[k] * tbe_0 * tke_0;

            to_y_z_y_xz[k] = to_0_x[k] - 2.0 * to_0_xzz[k] * tke_0 - 2.0 * to_yy_x[k] * tbe_0 + 4.0 * to_yy_xzz[k] * tbe_0 * tke_0;

            to_y_z_y_yy[k] = -2.0 * to_0_yyz[k] * tke_0 + 4.0 * to_yy_yyz[k] * tbe_0 * tke_0;

            to_y_z_y_yz[k] = to_0_y[k] - 2.0 * to_0_yzz[k] * tke_0 - 2.0 * to_yy_y[k] * tbe_0 + 4.0 * to_yy_yzz[k] * tbe_0 * tke_0;

            to_y_z_y_zz[k] = 2.0 * to_0_z[k] - 2.0 * to_0_zzz[k] * tke_0 - 4.0 * to_yy_z[k] * tbe_0 + 4.0 * to_yy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 102-108 components of targeted buffer : PD

        auto to_y_z_z_xx = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 12);

        auto to_y_z_z_xy = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 13);

        auto to_y_z_z_xz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 14);

        auto to_y_z_z_yy = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 15);

        auto to_y_z_z_yz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 16);

        auto to_y_z_z_zz = pbuffer.data(idx_op_geom_101_pd + 5 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_y_z_z_xx,     \
                             to_y_z_z_xy, \
                             to_y_z_z_xz, \
                             to_y_z_z_yy, \
                             to_y_z_z_yz, \
                             to_y_z_z_zz, \
                             to_yz_x,     \
                             to_yz_xxz,   \
                             to_yz_xyz,   \
                             to_yz_xzz,   \
                             to_yz_y,     \
                             to_yz_yyz,   \
                             to_yz_yzz,   \
                             to_yz_z,     \
                             to_yz_zzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_z_xx[k] = 4.0 * to_yz_xxz[k] * tbe_0 * tke_0;

            to_y_z_z_xy[k] = 4.0 * to_yz_xyz[k] * tbe_0 * tke_0;

            to_y_z_z_xz[k] = -2.0 * to_yz_x[k] * tbe_0 + 4.0 * to_yz_xzz[k] * tbe_0 * tke_0;

            to_y_z_z_yy[k] = 4.0 * to_yz_yyz[k] * tbe_0 * tke_0;

            to_y_z_z_yz[k] = -2.0 * to_yz_y[k] * tbe_0 + 4.0 * to_yz_yzz[k] * tbe_0 * tke_0;

            to_y_z_z_zz[k] = -4.0 * to_yz_z[k] * tbe_0 + 4.0 * to_yz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 108-114 components of targeted buffer : PD

        auto to_z_x_x_xx = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 0);

        auto to_z_x_x_xy = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 1);

        auto to_z_x_x_xz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 2);

        auto to_z_x_x_yy = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 3);

        auto to_z_x_x_yz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 4);

        auto to_z_x_x_zz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_xz_x,         \
                             to_xz_xxx,   \
                             to_xz_xxy,   \
                             to_xz_xxz,   \
                             to_xz_xyy,   \
                             to_xz_xyz,   \
                             to_xz_xzz,   \
                             to_xz_y,     \
                             to_xz_z,     \
                             to_z_x_x_xx, \
                             to_z_x_x_xy, \
                             to_z_x_x_xz, \
                             to_z_x_x_yy, \
                             to_z_x_x_yz, \
                             to_z_x_x_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_x_xx[k] = -4.0 * to_xz_x[k] * tbe_0 + 4.0 * to_xz_xxx[k] * tbe_0 * tke_0;

            to_z_x_x_xy[k] = -2.0 * to_xz_y[k] * tbe_0 + 4.0 * to_xz_xxy[k] * tbe_0 * tke_0;

            to_z_x_x_xz[k] = -2.0 * to_xz_z[k] * tbe_0 + 4.0 * to_xz_xxz[k] * tbe_0 * tke_0;

            to_z_x_x_yy[k] = 4.0 * to_xz_xyy[k] * tbe_0 * tke_0;

            to_z_x_x_yz[k] = 4.0 * to_xz_xyz[k] * tbe_0 * tke_0;

            to_z_x_x_zz[k] = 4.0 * to_xz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 114-120 components of targeted buffer : PD

        auto to_z_x_y_xx = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 6);

        auto to_z_x_y_xy = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 7);

        auto to_z_x_y_xz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 8);

        auto to_z_x_y_yy = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 9);

        auto to_z_x_y_yz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 10);

        auto to_z_x_y_zz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_yz_x,         \
                             to_yz_xxx,   \
                             to_yz_xxy,   \
                             to_yz_xxz,   \
                             to_yz_xyy,   \
                             to_yz_xyz,   \
                             to_yz_xzz,   \
                             to_yz_y,     \
                             to_yz_z,     \
                             to_z_x_y_xx, \
                             to_z_x_y_xy, \
                             to_z_x_y_xz, \
                             to_z_x_y_yy, \
                             to_z_x_y_yz, \
                             to_z_x_y_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_y_xx[k] = -4.0 * to_yz_x[k] * tbe_0 + 4.0 * to_yz_xxx[k] * tbe_0 * tke_0;

            to_z_x_y_xy[k] = -2.0 * to_yz_y[k] * tbe_0 + 4.0 * to_yz_xxy[k] * tbe_0 * tke_0;

            to_z_x_y_xz[k] = -2.0 * to_yz_z[k] * tbe_0 + 4.0 * to_yz_xxz[k] * tbe_0 * tke_0;

            to_z_x_y_yy[k] = 4.0 * to_yz_xyy[k] * tbe_0 * tke_0;

            to_z_x_y_yz[k] = 4.0 * to_yz_xyz[k] * tbe_0 * tke_0;

            to_z_x_y_zz[k] = 4.0 * to_yz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-126 components of targeted buffer : PD

        auto to_z_x_z_xx = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 12);

        auto to_z_x_z_xy = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 13);

        auto to_z_x_z_xz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 14);

        auto to_z_x_z_yy = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 15);

        auto to_z_x_z_yz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 16);

        auto to_z_x_z_zz = pbuffer.data(idx_op_geom_101_pd + 6 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxx,    \
                             to_0_xxy,    \
                             to_0_xxz,    \
                             to_0_xyy,    \
                             to_0_xyz,    \
                             to_0_xzz,    \
                             to_0_y,      \
                             to_0_z,      \
                             to_z_x_z_xx, \
                             to_z_x_z_xy, \
                             to_z_x_z_xz, \
                             to_z_x_z_yy, \
                             to_z_x_z_yz, \
                             to_z_x_z_zz, \
                             to_zz_x,     \
                             to_zz_xxx,   \
                             to_zz_xxy,   \
                             to_zz_xxz,   \
                             to_zz_xyy,   \
                             to_zz_xyz,   \
                             to_zz_xzz,   \
                             to_zz_y,     \
                             to_zz_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_z_xx[k] = 2.0 * to_0_x[k] - 2.0 * to_0_xxx[k] * tke_0 - 4.0 * to_zz_x[k] * tbe_0 + 4.0 * to_zz_xxx[k] * tbe_0 * tke_0;

            to_z_x_z_xy[k] = to_0_y[k] - 2.0 * to_0_xxy[k] * tke_0 - 2.0 * to_zz_y[k] * tbe_0 + 4.0 * to_zz_xxy[k] * tbe_0 * tke_0;

            to_z_x_z_xz[k] = to_0_z[k] - 2.0 * to_0_xxz[k] * tke_0 - 2.0 * to_zz_z[k] * tbe_0 + 4.0 * to_zz_xxz[k] * tbe_0 * tke_0;

            to_z_x_z_yy[k] = -2.0 * to_0_xyy[k] * tke_0 + 4.0 * to_zz_xyy[k] * tbe_0 * tke_0;

            to_z_x_z_yz[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_zz_xyz[k] * tbe_0 * tke_0;

            to_z_x_z_zz[k] = -2.0 * to_0_xzz[k] * tke_0 + 4.0 * to_zz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 126-132 components of targeted buffer : PD

        auto to_z_y_x_xx = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 0);

        auto to_z_y_x_xy = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 1);

        auto to_z_y_x_xz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 2);

        auto to_z_y_x_yy = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 3);

        auto to_z_y_x_yz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 4);

        auto to_z_y_x_zz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_xz_x,         \
                             to_xz_xxy,   \
                             to_xz_xyy,   \
                             to_xz_xyz,   \
                             to_xz_y,     \
                             to_xz_yyy,   \
                             to_xz_yyz,   \
                             to_xz_yzz,   \
                             to_xz_z,     \
                             to_z_y_x_xx, \
                             to_z_y_x_xy, \
                             to_z_y_x_xz, \
                             to_z_y_x_yy, \
                             to_z_y_x_yz, \
                             to_z_y_x_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_x_xx[k] = 4.0 * to_xz_xxy[k] * tbe_0 * tke_0;

            to_z_y_x_xy[k] = -2.0 * to_xz_x[k] * tbe_0 + 4.0 * to_xz_xyy[k] * tbe_0 * tke_0;

            to_z_y_x_xz[k] = 4.0 * to_xz_xyz[k] * tbe_0 * tke_0;

            to_z_y_x_yy[k] = -4.0 * to_xz_y[k] * tbe_0 + 4.0 * to_xz_yyy[k] * tbe_0 * tke_0;

            to_z_y_x_yz[k] = -2.0 * to_xz_z[k] * tbe_0 + 4.0 * to_xz_yyz[k] * tbe_0 * tke_0;

            to_z_y_x_zz[k] = 4.0 * to_xz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 132-138 components of targeted buffer : PD

        auto to_z_y_y_xx = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 6);

        auto to_z_y_y_xy = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 7);

        auto to_z_y_y_xz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 8);

        auto to_z_y_y_yy = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 9);

        auto to_z_y_y_yz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 10);

        auto to_z_y_y_zz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_yz_x,         \
                             to_yz_xxy,   \
                             to_yz_xyy,   \
                             to_yz_xyz,   \
                             to_yz_y,     \
                             to_yz_yyy,   \
                             to_yz_yyz,   \
                             to_yz_yzz,   \
                             to_yz_z,     \
                             to_z_y_y_xx, \
                             to_z_y_y_xy, \
                             to_z_y_y_xz, \
                             to_z_y_y_yy, \
                             to_z_y_y_yz, \
                             to_z_y_y_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_y_xx[k] = 4.0 * to_yz_xxy[k] * tbe_0 * tke_0;

            to_z_y_y_xy[k] = -2.0 * to_yz_x[k] * tbe_0 + 4.0 * to_yz_xyy[k] * tbe_0 * tke_0;

            to_z_y_y_xz[k] = 4.0 * to_yz_xyz[k] * tbe_0 * tke_0;

            to_z_y_y_yy[k] = -4.0 * to_yz_y[k] * tbe_0 + 4.0 * to_yz_yyy[k] * tbe_0 * tke_0;

            to_z_y_y_yz[k] = -2.0 * to_yz_z[k] * tbe_0 + 4.0 * to_yz_yyz[k] * tbe_0 * tke_0;

            to_z_y_y_zz[k] = 4.0 * to_yz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 138-144 components of targeted buffer : PD

        auto to_z_y_z_xx = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 12);

        auto to_z_y_z_xy = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 13);

        auto to_z_y_z_xz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 14);

        auto to_z_y_z_yy = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 15);

        auto to_z_y_z_yz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 16);

        auto to_z_y_z_zz = pbuffer.data(idx_op_geom_101_pd + 7 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxy,    \
                             to_0_xyy,    \
                             to_0_xyz,    \
                             to_0_y,      \
                             to_0_yyy,    \
                             to_0_yyz,    \
                             to_0_yzz,    \
                             to_0_z,      \
                             to_z_y_z_xx, \
                             to_z_y_z_xy, \
                             to_z_y_z_xz, \
                             to_z_y_z_yy, \
                             to_z_y_z_yz, \
                             to_z_y_z_zz, \
                             to_zz_x,     \
                             to_zz_xxy,   \
                             to_zz_xyy,   \
                             to_zz_xyz,   \
                             to_zz_y,     \
                             to_zz_yyy,   \
                             to_zz_yyz,   \
                             to_zz_yzz,   \
                             to_zz_z,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_z_xx[k] = -2.0 * to_0_xxy[k] * tke_0 + 4.0 * to_zz_xxy[k] * tbe_0 * tke_0;

            to_z_y_z_xy[k] = to_0_x[k] - 2.0 * to_0_xyy[k] * tke_0 - 2.0 * to_zz_x[k] * tbe_0 + 4.0 * to_zz_xyy[k] * tbe_0 * tke_0;

            to_z_y_z_xz[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_zz_xyz[k] * tbe_0 * tke_0;

            to_z_y_z_yy[k] = 2.0 * to_0_y[k] - 2.0 * to_0_yyy[k] * tke_0 - 4.0 * to_zz_y[k] * tbe_0 + 4.0 * to_zz_yyy[k] * tbe_0 * tke_0;

            to_z_y_z_yz[k] = to_0_z[k] - 2.0 * to_0_yyz[k] * tke_0 - 2.0 * to_zz_z[k] * tbe_0 + 4.0 * to_zz_yyz[k] * tbe_0 * tke_0;

            to_z_y_z_zz[k] = -2.0 * to_0_yzz[k] * tke_0 + 4.0 * to_zz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 144-150 components of targeted buffer : PD

        auto to_z_z_x_xx = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 0);

        auto to_z_z_x_xy = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 1);

        auto to_z_z_x_xz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 2);

        auto to_z_z_x_yy = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 3);

        auto to_z_z_x_yz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 4);

        auto to_z_z_x_zz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 5);

#pragma omp simd aligned(to_xz_x,         \
                             to_xz_xxz,   \
                             to_xz_xyz,   \
                             to_xz_xzz,   \
                             to_xz_y,     \
                             to_xz_yyz,   \
                             to_xz_yzz,   \
                             to_xz_z,     \
                             to_xz_zzz,   \
                             to_z_z_x_xx, \
                             to_z_z_x_xy, \
                             to_z_z_x_xz, \
                             to_z_z_x_yy, \
                             to_z_z_x_yz, \
                             to_z_z_x_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_x_xx[k] = 4.0 * to_xz_xxz[k] * tbe_0 * tke_0;

            to_z_z_x_xy[k] = 4.0 * to_xz_xyz[k] * tbe_0 * tke_0;

            to_z_z_x_xz[k] = -2.0 * to_xz_x[k] * tbe_0 + 4.0 * to_xz_xzz[k] * tbe_0 * tke_0;

            to_z_z_x_yy[k] = 4.0 * to_xz_yyz[k] * tbe_0 * tke_0;

            to_z_z_x_yz[k] = -2.0 * to_xz_y[k] * tbe_0 + 4.0 * to_xz_yzz[k] * tbe_0 * tke_0;

            to_z_z_x_zz[k] = -4.0 * to_xz_z[k] * tbe_0 + 4.0 * to_xz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-156 components of targeted buffer : PD

        auto to_z_z_y_xx = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 6);

        auto to_z_z_y_xy = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 7);

        auto to_z_z_y_xz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 8);

        auto to_z_z_y_yy = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 9);

        auto to_z_z_y_yz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 10);

        auto to_z_z_y_zz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 11);

#pragma omp simd aligned(to_yz_x,         \
                             to_yz_xxz,   \
                             to_yz_xyz,   \
                             to_yz_xzz,   \
                             to_yz_y,     \
                             to_yz_yyz,   \
                             to_yz_yzz,   \
                             to_yz_z,     \
                             to_yz_zzz,   \
                             to_z_z_y_xx, \
                             to_z_z_y_xy, \
                             to_z_z_y_xz, \
                             to_z_z_y_yy, \
                             to_z_z_y_yz, \
                             to_z_z_y_zz, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_y_xx[k] = 4.0 * to_yz_xxz[k] * tbe_0 * tke_0;

            to_z_z_y_xy[k] = 4.0 * to_yz_xyz[k] * tbe_0 * tke_0;

            to_z_z_y_xz[k] = -2.0 * to_yz_x[k] * tbe_0 + 4.0 * to_yz_xzz[k] * tbe_0 * tke_0;

            to_z_z_y_yy[k] = 4.0 * to_yz_yyz[k] * tbe_0 * tke_0;

            to_z_z_y_yz[k] = -2.0 * to_yz_y[k] * tbe_0 + 4.0 * to_yz_yzz[k] * tbe_0 * tke_0;

            to_z_z_y_zz[k] = -4.0 * to_yz_z[k] * tbe_0 + 4.0 * to_yz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 156-162 components of targeted buffer : PD

        auto to_z_z_z_xx = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 12);

        auto to_z_z_z_xy = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 13);

        auto to_z_z_z_xz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 14);

        auto to_z_z_z_yy = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 15);

        auto to_z_z_z_yz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 16);

        auto to_z_z_z_zz = pbuffer.data(idx_op_geom_101_pd + 8 * op_comps * 18 + i * 18 + 17);

#pragma omp simd aligned(to_0_x,          \
                             to_0_xxz,    \
                             to_0_xyz,    \
                             to_0_xzz,    \
                             to_0_y,      \
                             to_0_yyz,    \
                             to_0_yzz,    \
                             to_0_z,      \
                             to_0_zzz,    \
                             to_z_z_z_xx, \
                             to_z_z_z_xy, \
                             to_z_z_z_xz, \
                             to_z_z_z_yy, \
                             to_z_z_z_yz, \
                             to_z_z_z_zz, \
                             to_zz_x,     \
                             to_zz_xxz,   \
                             to_zz_xyz,   \
                             to_zz_xzz,   \
                             to_zz_y,     \
                             to_zz_yyz,   \
                             to_zz_yzz,   \
                             to_zz_z,     \
                             to_zz_zzz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_z_xx[k] = -2.0 * to_0_xxz[k] * tke_0 + 4.0 * to_zz_xxz[k] * tbe_0 * tke_0;

            to_z_z_z_xy[k] = -2.0 * to_0_xyz[k] * tke_0 + 4.0 * to_zz_xyz[k] * tbe_0 * tke_0;

            to_z_z_z_xz[k] = to_0_x[k] - 2.0 * to_0_xzz[k] * tke_0 - 2.0 * to_zz_x[k] * tbe_0 + 4.0 * to_zz_xzz[k] * tbe_0 * tke_0;

            to_z_z_z_yy[k] = -2.0 * to_0_yyz[k] * tke_0 + 4.0 * to_zz_yyz[k] * tbe_0 * tke_0;

            to_z_z_z_yz[k] = to_0_y[k] - 2.0 * to_0_yzz[k] * tke_0 - 2.0 * to_zz_y[k] * tbe_0 + 4.0 * to_zz_yzz[k] * tbe_0 * tke_0;

            to_z_z_z_zz[k] = 2.0 * to_0_z[k] - 2.0 * to_0_zzz[k] * tke_0 - 4.0 * to_zz_z[k] * tbe_0 + 4.0 * to_zz_zzz[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
