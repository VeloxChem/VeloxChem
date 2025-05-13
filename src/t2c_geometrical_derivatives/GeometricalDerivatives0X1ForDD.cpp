#include "GeometricalDerivatives0X1ForDD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_dd(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_dd,
                       const int idx_op_dp,
                       const int idx_op_df,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
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

        // Set up 0-6 components of targeted buffer : DD

        auto to_0_x_xx_xx = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 0);

        auto to_0_x_xx_xy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 1);

        auto to_0_x_xx_xz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 2);

        auto to_0_x_xx_yy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 3);

        auto to_0_x_xx_yz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 4);

        auto to_0_x_xx_zz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_0_x_xx_xx, to_0_x_xx_xy, to_0_x_xx_xz, to_0_x_xx_yy, to_0_x_xx_yz, to_0_x_xx_zz, to_xx_x, to_xx_xxx, to_xx_xxy, to_xx_xxz, to_xx_xyy, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xx_xx[k] = -2.0 * to_xx_x[k] + 2.0 * to_xx_xxx[k] * tke_0;

            to_0_x_xx_xy[k] = -to_xx_y[k] + 2.0 * to_xx_xxy[k] * tke_0;

            to_0_x_xx_xz[k] = -to_xx_z[k] + 2.0 * to_xx_xxz[k] * tke_0;

            to_0_x_xx_yy[k] = 2.0 * to_xx_xyy[k] * tke_0;

            to_0_x_xx_yz[k] = 2.0 * to_xx_xyz[k] * tke_0;

            to_0_x_xx_zz[k] = 2.0 * to_xx_xzz[k] * tke_0;
        }

        // Set up 6-12 components of targeted buffer : DD

        auto to_0_x_xy_xx = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 6);

        auto to_0_x_xy_xy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 7);

        auto to_0_x_xy_xz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 8);

        auto to_0_x_xy_yy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 9);

        auto to_0_x_xy_yz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 10);

        auto to_0_x_xy_zz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_0_x_xy_xx, to_0_x_xy_xy, to_0_x_xy_xz, to_0_x_xy_yy, to_0_x_xy_yz, to_0_x_xy_zz, to_xy_x, to_xy_xxx, to_xy_xxy, to_xy_xxz, to_xy_xyy, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xy_xx[k] = -2.0 * to_xy_x[k] + 2.0 * to_xy_xxx[k] * tke_0;

            to_0_x_xy_xy[k] = -to_xy_y[k] + 2.0 * to_xy_xxy[k] * tke_0;

            to_0_x_xy_xz[k] = -to_xy_z[k] + 2.0 * to_xy_xxz[k] * tke_0;

            to_0_x_xy_yy[k] = 2.0 * to_xy_xyy[k] * tke_0;

            to_0_x_xy_yz[k] = 2.0 * to_xy_xyz[k] * tke_0;

            to_0_x_xy_zz[k] = 2.0 * to_xy_xzz[k] * tke_0;
        }

        // Set up 12-18 components of targeted buffer : DD

        auto to_0_x_xz_xx = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 12);

        auto to_0_x_xz_xy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 13);

        auto to_0_x_xz_xz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 14);

        auto to_0_x_xz_yy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 15);

        auto to_0_x_xz_yz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 16);

        auto to_0_x_xz_zz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_0_x_xz_xx, to_0_x_xz_xy, to_0_x_xz_xz, to_0_x_xz_yy, to_0_x_xz_yz, to_0_x_xz_zz, to_xz_x, to_xz_xxx, to_xz_xxy, to_xz_xxz, to_xz_xyy, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xz_xx[k] = -2.0 * to_xz_x[k] + 2.0 * to_xz_xxx[k] * tke_0;

            to_0_x_xz_xy[k] = -to_xz_y[k] + 2.0 * to_xz_xxy[k] * tke_0;

            to_0_x_xz_xz[k] = -to_xz_z[k] + 2.0 * to_xz_xxz[k] * tke_0;

            to_0_x_xz_yy[k] = 2.0 * to_xz_xyy[k] * tke_0;

            to_0_x_xz_yz[k] = 2.0 * to_xz_xyz[k] * tke_0;

            to_0_x_xz_zz[k] = 2.0 * to_xz_xzz[k] * tke_0;
        }

        // Set up 18-24 components of targeted buffer : DD

        auto to_0_x_yy_xx = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 18);

        auto to_0_x_yy_xy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 19);

        auto to_0_x_yy_xz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 20);

        auto to_0_x_yy_yy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 21);

        auto to_0_x_yy_yz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 22);

        auto to_0_x_yy_zz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_0_x_yy_xx, to_0_x_yy_xy, to_0_x_yy_xz, to_0_x_yy_yy, to_0_x_yy_yz, to_0_x_yy_zz, to_yy_x, to_yy_xxx, to_yy_xxy, to_yy_xxz, to_yy_xyy, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yy_xx[k] = -2.0 * to_yy_x[k] + 2.0 * to_yy_xxx[k] * tke_0;

            to_0_x_yy_xy[k] = -to_yy_y[k] + 2.0 * to_yy_xxy[k] * tke_0;

            to_0_x_yy_xz[k] = -to_yy_z[k] + 2.0 * to_yy_xxz[k] * tke_0;

            to_0_x_yy_yy[k] = 2.0 * to_yy_xyy[k] * tke_0;

            to_0_x_yy_yz[k] = 2.0 * to_yy_xyz[k] * tke_0;

            to_0_x_yy_zz[k] = 2.0 * to_yy_xzz[k] * tke_0;
        }

        // Set up 24-30 components of targeted buffer : DD

        auto to_0_x_yz_xx = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 24);

        auto to_0_x_yz_xy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 25);

        auto to_0_x_yz_xz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 26);

        auto to_0_x_yz_yy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 27);

        auto to_0_x_yz_yz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 28);

        auto to_0_x_yz_zz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_0_x_yz_xx, to_0_x_yz_xy, to_0_x_yz_xz, to_0_x_yz_yy, to_0_x_yz_yz, to_0_x_yz_zz, to_yz_x, to_yz_xxx, to_yz_xxy, to_yz_xxz, to_yz_xyy, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yz_xx[k] = -2.0 * to_yz_x[k] + 2.0 * to_yz_xxx[k] * tke_0;

            to_0_x_yz_xy[k] = -to_yz_y[k] + 2.0 * to_yz_xxy[k] * tke_0;

            to_0_x_yz_xz[k] = -to_yz_z[k] + 2.0 * to_yz_xxz[k] * tke_0;

            to_0_x_yz_yy[k] = 2.0 * to_yz_xyy[k] * tke_0;

            to_0_x_yz_yz[k] = 2.0 * to_yz_xyz[k] * tke_0;

            to_0_x_yz_zz[k] = 2.0 * to_yz_xzz[k] * tke_0;
        }

        // Set up 30-36 components of targeted buffer : DD

        auto to_0_x_zz_xx = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 30);

        auto to_0_x_zz_xy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 31);

        auto to_0_x_zz_xz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 32);

        auto to_0_x_zz_yy = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 33);

        auto to_0_x_zz_yz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 34);

        auto to_0_x_zz_zz = pbuffer.data(idx_op_geom_001_dd + 0 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_0_x_zz_xx, to_0_x_zz_xy, to_0_x_zz_xz, to_0_x_zz_yy, to_0_x_zz_yz, to_0_x_zz_zz, to_zz_x, to_zz_xxx, to_zz_xxy, to_zz_xxz, to_zz_xyy, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zz_xx[k] = -2.0 * to_zz_x[k] + 2.0 * to_zz_xxx[k] * tke_0;

            to_0_x_zz_xy[k] = -to_zz_y[k] + 2.0 * to_zz_xxy[k] * tke_0;

            to_0_x_zz_xz[k] = -to_zz_z[k] + 2.0 * to_zz_xxz[k] * tke_0;

            to_0_x_zz_yy[k] = 2.0 * to_zz_xyy[k] * tke_0;

            to_0_x_zz_yz[k] = 2.0 * to_zz_xyz[k] * tke_0;

            to_0_x_zz_zz[k] = 2.0 * to_zz_xzz[k] * tke_0;
        }

        // Set up 36-42 components of targeted buffer : DD

        auto to_0_y_xx_xx = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 0);

        auto to_0_y_xx_xy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 1);

        auto to_0_y_xx_xz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 2);

        auto to_0_y_xx_yy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 3);

        auto to_0_y_xx_yz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 4);

        auto to_0_y_xx_zz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_0_y_xx_xx, to_0_y_xx_xy, to_0_y_xx_xz, to_0_y_xx_yy, to_0_y_xx_yz, to_0_y_xx_zz, to_xx_x, to_xx_xxy, to_xx_xyy, to_xx_xyz, to_xx_y, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xx_xx[k] = 2.0 * to_xx_xxy[k] * tke_0;

            to_0_y_xx_xy[k] = -to_xx_x[k] + 2.0 * to_xx_xyy[k] * tke_0;

            to_0_y_xx_xz[k] = 2.0 * to_xx_xyz[k] * tke_0;

            to_0_y_xx_yy[k] = -2.0 * to_xx_y[k] + 2.0 * to_xx_yyy[k] * tke_0;

            to_0_y_xx_yz[k] = -to_xx_z[k] + 2.0 * to_xx_yyz[k] * tke_0;

            to_0_y_xx_zz[k] = 2.0 * to_xx_yzz[k] * tke_0;
        }

        // Set up 42-48 components of targeted buffer : DD

        auto to_0_y_xy_xx = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 6);

        auto to_0_y_xy_xy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 7);

        auto to_0_y_xy_xz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 8);

        auto to_0_y_xy_yy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 9);

        auto to_0_y_xy_yz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 10);

        auto to_0_y_xy_zz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_0_y_xy_xx, to_0_y_xy_xy, to_0_y_xy_xz, to_0_y_xy_yy, to_0_y_xy_yz, to_0_y_xy_zz, to_xy_x, to_xy_xxy, to_xy_xyy, to_xy_xyz, to_xy_y, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xy_xx[k] = 2.0 * to_xy_xxy[k] * tke_0;

            to_0_y_xy_xy[k] = -to_xy_x[k] + 2.0 * to_xy_xyy[k] * tke_0;

            to_0_y_xy_xz[k] = 2.0 * to_xy_xyz[k] * tke_0;

            to_0_y_xy_yy[k] = -2.0 * to_xy_y[k] + 2.0 * to_xy_yyy[k] * tke_0;

            to_0_y_xy_yz[k] = -to_xy_z[k] + 2.0 * to_xy_yyz[k] * tke_0;

            to_0_y_xy_zz[k] = 2.0 * to_xy_yzz[k] * tke_0;
        }

        // Set up 48-54 components of targeted buffer : DD

        auto to_0_y_xz_xx = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 12);

        auto to_0_y_xz_xy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 13);

        auto to_0_y_xz_xz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 14);

        auto to_0_y_xz_yy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 15);

        auto to_0_y_xz_yz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 16);

        auto to_0_y_xz_zz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_0_y_xz_xx, to_0_y_xz_xy, to_0_y_xz_xz, to_0_y_xz_yy, to_0_y_xz_yz, to_0_y_xz_zz, to_xz_x, to_xz_xxy, to_xz_xyy, to_xz_xyz, to_xz_y, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xz_xx[k] = 2.0 * to_xz_xxy[k] * tke_0;

            to_0_y_xz_xy[k] = -to_xz_x[k] + 2.0 * to_xz_xyy[k] * tke_0;

            to_0_y_xz_xz[k] = 2.0 * to_xz_xyz[k] * tke_0;

            to_0_y_xz_yy[k] = -2.0 * to_xz_y[k] + 2.0 * to_xz_yyy[k] * tke_0;

            to_0_y_xz_yz[k] = -to_xz_z[k] + 2.0 * to_xz_yyz[k] * tke_0;

            to_0_y_xz_zz[k] = 2.0 * to_xz_yzz[k] * tke_0;
        }

        // Set up 54-60 components of targeted buffer : DD

        auto to_0_y_yy_xx = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 18);

        auto to_0_y_yy_xy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 19);

        auto to_0_y_yy_xz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 20);

        auto to_0_y_yy_yy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 21);

        auto to_0_y_yy_yz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 22);

        auto to_0_y_yy_zz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_0_y_yy_xx, to_0_y_yy_xy, to_0_y_yy_xz, to_0_y_yy_yy, to_0_y_yy_yz, to_0_y_yy_zz, to_yy_x, to_yy_xxy, to_yy_xyy, to_yy_xyz, to_yy_y, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yy_xx[k] = 2.0 * to_yy_xxy[k] * tke_0;

            to_0_y_yy_xy[k] = -to_yy_x[k] + 2.0 * to_yy_xyy[k] * tke_0;

            to_0_y_yy_xz[k] = 2.0 * to_yy_xyz[k] * tke_0;

            to_0_y_yy_yy[k] = -2.0 * to_yy_y[k] + 2.0 * to_yy_yyy[k] * tke_0;

            to_0_y_yy_yz[k] = -to_yy_z[k] + 2.0 * to_yy_yyz[k] * tke_0;

            to_0_y_yy_zz[k] = 2.0 * to_yy_yzz[k] * tke_0;
        }

        // Set up 60-66 components of targeted buffer : DD

        auto to_0_y_yz_xx = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 24);

        auto to_0_y_yz_xy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 25);

        auto to_0_y_yz_xz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 26);

        auto to_0_y_yz_yy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 27);

        auto to_0_y_yz_yz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 28);

        auto to_0_y_yz_zz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_0_y_yz_xx, to_0_y_yz_xy, to_0_y_yz_xz, to_0_y_yz_yy, to_0_y_yz_yz, to_0_y_yz_zz, to_yz_x, to_yz_xxy, to_yz_xyy, to_yz_xyz, to_yz_y, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yz_xx[k] = 2.0 * to_yz_xxy[k] * tke_0;

            to_0_y_yz_xy[k] = -to_yz_x[k] + 2.0 * to_yz_xyy[k] * tke_0;

            to_0_y_yz_xz[k] = 2.0 * to_yz_xyz[k] * tke_0;

            to_0_y_yz_yy[k] = -2.0 * to_yz_y[k] + 2.0 * to_yz_yyy[k] * tke_0;

            to_0_y_yz_yz[k] = -to_yz_z[k] + 2.0 * to_yz_yyz[k] * tke_0;

            to_0_y_yz_zz[k] = 2.0 * to_yz_yzz[k] * tke_0;
        }

        // Set up 66-72 components of targeted buffer : DD

        auto to_0_y_zz_xx = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 30);

        auto to_0_y_zz_xy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 31);

        auto to_0_y_zz_xz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 32);

        auto to_0_y_zz_yy = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 33);

        auto to_0_y_zz_yz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 34);

        auto to_0_y_zz_zz = pbuffer.data(idx_op_geom_001_dd + 1 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_0_y_zz_xx, to_0_y_zz_xy, to_0_y_zz_xz, to_0_y_zz_yy, to_0_y_zz_yz, to_0_y_zz_zz, to_zz_x, to_zz_xxy, to_zz_xyy, to_zz_xyz, to_zz_y, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zz_xx[k] = 2.0 * to_zz_xxy[k] * tke_0;

            to_0_y_zz_xy[k] = -to_zz_x[k] + 2.0 * to_zz_xyy[k] * tke_0;

            to_0_y_zz_xz[k] = 2.0 * to_zz_xyz[k] * tke_0;

            to_0_y_zz_yy[k] = -2.0 * to_zz_y[k] + 2.0 * to_zz_yyy[k] * tke_0;

            to_0_y_zz_yz[k] = -to_zz_z[k] + 2.0 * to_zz_yyz[k] * tke_0;

            to_0_y_zz_zz[k] = 2.0 * to_zz_yzz[k] * tke_0;
        }

        // Set up 72-78 components of targeted buffer : DD

        auto to_0_z_xx_xx = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 0);

        auto to_0_z_xx_xy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 1);

        auto to_0_z_xx_xz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 2);

        auto to_0_z_xx_yy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 3);

        auto to_0_z_xx_yz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 4);

        auto to_0_z_xx_zz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 5);

        #pragma omp simd aligned(to_0_z_xx_xx, to_0_z_xx_xy, to_0_z_xx_xz, to_0_z_xx_yy, to_0_z_xx_yz, to_0_z_xx_zz, to_xx_x, to_xx_xxz, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_yyz, to_xx_yzz, to_xx_z, to_xx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xx_xx[k] = 2.0 * to_xx_xxz[k] * tke_0;

            to_0_z_xx_xy[k] = 2.0 * to_xx_xyz[k] * tke_0;

            to_0_z_xx_xz[k] = -to_xx_x[k] + 2.0 * to_xx_xzz[k] * tke_0;

            to_0_z_xx_yy[k] = 2.0 * to_xx_yyz[k] * tke_0;

            to_0_z_xx_yz[k] = -to_xx_y[k] + 2.0 * to_xx_yzz[k] * tke_0;

            to_0_z_xx_zz[k] = -2.0 * to_xx_z[k] + 2.0 * to_xx_zzz[k] * tke_0;
        }

        // Set up 78-84 components of targeted buffer : DD

        auto to_0_z_xy_xx = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 6);

        auto to_0_z_xy_xy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 7);

        auto to_0_z_xy_xz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 8);

        auto to_0_z_xy_yy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 9);

        auto to_0_z_xy_yz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 10);

        auto to_0_z_xy_zz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 11);

        #pragma omp simd aligned(to_0_z_xy_xx, to_0_z_xy_xy, to_0_z_xy_xz, to_0_z_xy_yy, to_0_z_xy_yz, to_0_z_xy_zz, to_xy_x, to_xy_xxz, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_yyz, to_xy_yzz, to_xy_z, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xy_xx[k] = 2.0 * to_xy_xxz[k] * tke_0;

            to_0_z_xy_xy[k] = 2.0 * to_xy_xyz[k] * tke_0;

            to_0_z_xy_xz[k] = -to_xy_x[k] + 2.0 * to_xy_xzz[k] * tke_0;

            to_0_z_xy_yy[k] = 2.0 * to_xy_yyz[k] * tke_0;

            to_0_z_xy_yz[k] = -to_xy_y[k] + 2.0 * to_xy_yzz[k] * tke_0;

            to_0_z_xy_zz[k] = -2.0 * to_xy_z[k] + 2.0 * to_xy_zzz[k] * tke_0;
        }

        // Set up 84-90 components of targeted buffer : DD

        auto to_0_z_xz_xx = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 12);

        auto to_0_z_xz_xy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 13);

        auto to_0_z_xz_xz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 14);

        auto to_0_z_xz_yy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 15);

        auto to_0_z_xz_yz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 16);

        auto to_0_z_xz_zz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 17);

        #pragma omp simd aligned(to_0_z_xz_xx, to_0_z_xz_xy, to_0_z_xz_xz, to_0_z_xz_yy, to_0_z_xz_yz, to_0_z_xz_zz, to_xz_x, to_xz_xxz, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_yyz, to_xz_yzz, to_xz_z, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xz_xx[k] = 2.0 * to_xz_xxz[k] * tke_0;

            to_0_z_xz_xy[k] = 2.0 * to_xz_xyz[k] * tke_0;

            to_0_z_xz_xz[k] = -to_xz_x[k] + 2.0 * to_xz_xzz[k] * tke_0;

            to_0_z_xz_yy[k] = 2.0 * to_xz_yyz[k] * tke_0;

            to_0_z_xz_yz[k] = -to_xz_y[k] + 2.0 * to_xz_yzz[k] * tke_0;

            to_0_z_xz_zz[k] = -2.0 * to_xz_z[k] + 2.0 * to_xz_zzz[k] * tke_0;
        }

        // Set up 90-96 components of targeted buffer : DD

        auto to_0_z_yy_xx = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 18);

        auto to_0_z_yy_xy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 19);

        auto to_0_z_yy_xz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 20);

        auto to_0_z_yy_yy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 21);

        auto to_0_z_yy_yz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 22);

        auto to_0_z_yy_zz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 23);

        #pragma omp simd aligned(to_0_z_yy_xx, to_0_z_yy_xy, to_0_z_yy_xz, to_0_z_yy_yy, to_0_z_yy_yz, to_0_z_yy_zz, to_yy_x, to_yy_xxz, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_yyz, to_yy_yzz, to_yy_z, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yy_xx[k] = 2.0 * to_yy_xxz[k] * tke_0;

            to_0_z_yy_xy[k] = 2.0 * to_yy_xyz[k] * tke_0;

            to_0_z_yy_xz[k] = -to_yy_x[k] + 2.0 * to_yy_xzz[k] * tke_0;

            to_0_z_yy_yy[k] = 2.0 * to_yy_yyz[k] * tke_0;

            to_0_z_yy_yz[k] = -to_yy_y[k] + 2.0 * to_yy_yzz[k] * tke_0;

            to_0_z_yy_zz[k] = -2.0 * to_yy_z[k] + 2.0 * to_yy_zzz[k] * tke_0;
        }

        // Set up 96-102 components of targeted buffer : DD

        auto to_0_z_yz_xx = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 24);

        auto to_0_z_yz_xy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 25);

        auto to_0_z_yz_xz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 26);

        auto to_0_z_yz_yy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 27);

        auto to_0_z_yz_yz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 28);

        auto to_0_z_yz_zz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 29);

        #pragma omp simd aligned(to_0_z_yz_xx, to_0_z_yz_xy, to_0_z_yz_xz, to_0_z_yz_yy, to_0_z_yz_yz, to_0_z_yz_zz, to_yz_x, to_yz_xxz, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_yyz, to_yz_yzz, to_yz_z, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yz_xx[k] = 2.0 * to_yz_xxz[k] * tke_0;

            to_0_z_yz_xy[k] = 2.0 * to_yz_xyz[k] * tke_0;

            to_0_z_yz_xz[k] = -to_yz_x[k] + 2.0 * to_yz_xzz[k] * tke_0;

            to_0_z_yz_yy[k] = 2.0 * to_yz_yyz[k] * tke_0;

            to_0_z_yz_yz[k] = -to_yz_y[k] + 2.0 * to_yz_yzz[k] * tke_0;

            to_0_z_yz_zz[k] = -2.0 * to_yz_z[k] + 2.0 * to_yz_zzz[k] * tke_0;
        }

        // Set up 102-108 components of targeted buffer : DD

        auto to_0_z_zz_xx = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 30);

        auto to_0_z_zz_xy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 31);

        auto to_0_z_zz_xz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 32);

        auto to_0_z_zz_yy = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 33);

        auto to_0_z_zz_yz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 34);

        auto to_0_z_zz_zz = pbuffer.data(idx_op_geom_001_dd + 2 * op_comps * 36 + i * 36 + 35);

        #pragma omp simd aligned(to_0_z_zz_xx, to_0_z_zz_xy, to_0_z_zz_xz, to_0_z_zz_yy, to_0_z_zz_yz, to_0_z_zz_zz, to_zz_x, to_zz_xxz, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_yyz, to_zz_yzz, to_zz_z, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zz_xx[k] = 2.0 * to_zz_xxz[k] * tke_0;

            to_0_z_zz_xy[k] = 2.0 * to_zz_xyz[k] * tke_0;

            to_0_z_zz_xz[k] = -to_zz_x[k] + 2.0 * to_zz_xzz[k] * tke_0;

            to_0_z_zz_yy[k] = 2.0 * to_zz_yyz[k] * tke_0;

            to_0_z_zz_yz[k] = -to_zz_y[k] + 2.0 * to_zz_yzz[k] * tke_0;

            to_0_z_zz_zz[k] = -2.0 * to_zz_z[k] + 2.0 * to_zz_zzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

