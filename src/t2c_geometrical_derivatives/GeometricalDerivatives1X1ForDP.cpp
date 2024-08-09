#include "GeometricalDerivatives1X1ForDP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_dp(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_dp,
                        const size_t idx_op_ps,
                        const size_t idx_op_pd,
                        const size_t idx_op_fs,
                        const size_t idx_op_fd,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PS

        auto to_x_0 = pbuffer.data(idx_op_ps + i * 3 + 0);

        auto to_y_0 = pbuffer.data(idx_op_ps + i * 3 + 1);

        auto to_z_0 = pbuffer.data(idx_op_ps + i * 3 + 2);

        // Set up components of auxiliary buffer : PD

        auto to_x_xx = pbuffer.data(idx_op_pd + i * 18 + 0);

        auto to_x_xy = pbuffer.data(idx_op_pd + i * 18 + 1);

        auto to_x_xz = pbuffer.data(idx_op_pd + i * 18 + 2);

        auto to_x_yy = pbuffer.data(idx_op_pd + i * 18 + 3);

        auto to_x_yz = pbuffer.data(idx_op_pd + i * 18 + 4);

        auto to_x_zz = pbuffer.data(idx_op_pd + i * 18 + 5);

        auto to_y_xx = pbuffer.data(idx_op_pd + i * 18 + 6);

        auto to_y_xy = pbuffer.data(idx_op_pd + i * 18 + 7);

        auto to_y_xz = pbuffer.data(idx_op_pd + i * 18 + 8);

        auto to_y_yy = pbuffer.data(idx_op_pd + i * 18 + 9);

        auto to_y_yz = pbuffer.data(idx_op_pd + i * 18 + 10);

        auto to_y_zz = pbuffer.data(idx_op_pd + i * 18 + 11);

        auto to_z_xx = pbuffer.data(idx_op_pd + i * 18 + 12);

        auto to_z_xy = pbuffer.data(idx_op_pd + i * 18 + 13);

        auto to_z_xz = pbuffer.data(idx_op_pd + i * 18 + 14);

        auto to_z_yy = pbuffer.data(idx_op_pd + i * 18 + 15);

        auto to_z_yz = pbuffer.data(idx_op_pd + i * 18 + 16);

        auto to_z_zz = pbuffer.data(idx_op_pd + i * 18 + 17);

        // Set up components of auxiliary buffer : FS

        auto to_xxx_0 = pbuffer.data(idx_op_fs + i * 10 + 0);

        auto to_xxy_0 = pbuffer.data(idx_op_fs + i * 10 + 1);

        auto to_xxz_0 = pbuffer.data(idx_op_fs + i * 10 + 2);

        auto to_xyy_0 = pbuffer.data(idx_op_fs + i * 10 + 3);

        auto to_xyz_0 = pbuffer.data(idx_op_fs + i * 10 + 4);

        auto to_xzz_0 = pbuffer.data(idx_op_fs + i * 10 + 5);

        auto to_yyy_0 = pbuffer.data(idx_op_fs + i * 10 + 6);

        auto to_yyz_0 = pbuffer.data(idx_op_fs + i * 10 + 7);

        auto to_yzz_0 = pbuffer.data(idx_op_fs + i * 10 + 8);

        auto to_zzz_0 = pbuffer.data(idx_op_fs + i * 10 + 9);

        // Set up components of auxiliary buffer : FD

        auto to_xxx_xx = pbuffer.data(idx_op_fd + i * 60 + 0);

        auto to_xxx_xy = pbuffer.data(idx_op_fd + i * 60 + 1);

        auto to_xxx_xz = pbuffer.data(idx_op_fd + i * 60 + 2);

        auto to_xxx_yy = pbuffer.data(idx_op_fd + i * 60 + 3);

        auto to_xxx_yz = pbuffer.data(idx_op_fd + i * 60 + 4);

        auto to_xxx_zz = pbuffer.data(idx_op_fd + i * 60 + 5);

        auto to_xxy_xx = pbuffer.data(idx_op_fd + i * 60 + 6);

        auto to_xxy_xy = pbuffer.data(idx_op_fd + i * 60 + 7);

        auto to_xxy_xz = pbuffer.data(idx_op_fd + i * 60 + 8);

        auto to_xxy_yy = pbuffer.data(idx_op_fd + i * 60 + 9);

        auto to_xxy_yz = pbuffer.data(idx_op_fd + i * 60 + 10);

        auto to_xxy_zz = pbuffer.data(idx_op_fd + i * 60 + 11);

        auto to_xxz_xx = pbuffer.data(idx_op_fd + i * 60 + 12);

        auto to_xxz_xy = pbuffer.data(idx_op_fd + i * 60 + 13);

        auto to_xxz_xz = pbuffer.data(idx_op_fd + i * 60 + 14);

        auto to_xxz_yy = pbuffer.data(idx_op_fd + i * 60 + 15);

        auto to_xxz_yz = pbuffer.data(idx_op_fd + i * 60 + 16);

        auto to_xxz_zz = pbuffer.data(idx_op_fd + i * 60 + 17);

        auto to_xyy_xx = pbuffer.data(idx_op_fd + i * 60 + 18);

        auto to_xyy_xy = pbuffer.data(idx_op_fd + i * 60 + 19);

        auto to_xyy_xz = pbuffer.data(idx_op_fd + i * 60 + 20);

        auto to_xyy_yy = pbuffer.data(idx_op_fd + i * 60 + 21);

        auto to_xyy_yz = pbuffer.data(idx_op_fd + i * 60 + 22);

        auto to_xyy_zz = pbuffer.data(idx_op_fd + i * 60 + 23);

        auto to_xyz_xx = pbuffer.data(idx_op_fd + i * 60 + 24);

        auto to_xyz_xy = pbuffer.data(idx_op_fd + i * 60 + 25);

        auto to_xyz_xz = pbuffer.data(idx_op_fd + i * 60 + 26);

        auto to_xyz_yy = pbuffer.data(idx_op_fd + i * 60 + 27);

        auto to_xyz_yz = pbuffer.data(idx_op_fd + i * 60 + 28);

        auto to_xyz_zz = pbuffer.data(idx_op_fd + i * 60 + 29);

        auto to_xzz_xx = pbuffer.data(idx_op_fd + i * 60 + 30);

        auto to_xzz_xy = pbuffer.data(idx_op_fd + i * 60 + 31);

        auto to_xzz_xz = pbuffer.data(idx_op_fd + i * 60 + 32);

        auto to_xzz_yy = pbuffer.data(idx_op_fd + i * 60 + 33);

        auto to_xzz_yz = pbuffer.data(idx_op_fd + i * 60 + 34);

        auto to_xzz_zz = pbuffer.data(idx_op_fd + i * 60 + 35);

        auto to_yyy_xx = pbuffer.data(idx_op_fd + i * 60 + 36);

        auto to_yyy_xy = pbuffer.data(idx_op_fd + i * 60 + 37);

        auto to_yyy_xz = pbuffer.data(idx_op_fd + i * 60 + 38);

        auto to_yyy_yy = pbuffer.data(idx_op_fd + i * 60 + 39);

        auto to_yyy_yz = pbuffer.data(idx_op_fd + i * 60 + 40);

        auto to_yyy_zz = pbuffer.data(idx_op_fd + i * 60 + 41);

        auto to_yyz_xx = pbuffer.data(idx_op_fd + i * 60 + 42);

        auto to_yyz_xy = pbuffer.data(idx_op_fd + i * 60 + 43);

        auto to_yyz_xz = pbuffer.data(idx_op_fd + i * 60 + 44);

        auto to_yyz_yy = pbuffer.data(idx_op_fd + i * 60 + 45);

        auto to_yyz_yz = pbuffer.data(idx_op_fd + i * 60 + 46);

        auto to_yyz_zz = pbuffer.data(idx_op_fd + i * 60 + 47);

        auto to_yzz_xx = pbuffer.data(idx_op_fd + i * 60 + 48);

        auto to_yzz_xy = pbuffer.data(idx_op_fd + i * 60 + 49);

        auto to_yzz_xz = pbuffer.data(idx_op_fd + i * 60 + 50);

        auto to_yzz_yy = pbuffer.data(idx_op_fd + i * 60 + 51);

        auto to_yzz_yz = pbuffer.data(idx_op_fd + i * 60 + 52);

        auto to_yzz_zz = pbuffer.data(idx_op_fd + i * 60 + 53);

        auto to_zzz_xx = pbuffer.data(idx_op_fd + i * 60 + 54);

        auto to_zzz_xy = pbuffer.data(idx_op_fd + i * 60 + 55);

        auto to_zzz_xz = pbuffer.data(idx_op_fd + i * 60 + 56);

        auto to_zzz_yy = pbuffer.data(idx_op_fd + i * 60 + 57);

        auto to_zzz_yz = pbuffer.data(idx_op_fd + i * 60 + 58);

        auto to_zzz_zz = pbuffer.data(idx_op_fd + i * 60 + 59);

        // Set up 0-3 components of targeted buffer : DP

        auto to_x_x_xx_x = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 0);

        auto to_x_x_xx_y = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 1);

        auto to_x_x_xx_z = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_x_0, to_x_x_xx_x, to_x_x_xx_y, to_x_x_xx_z, to_x_xx, to_x_xy, to_x_xz, to_xxx_0, to_xxx_xx, to_xxx_xy, to_xxx_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xx_x[k] = 2.0 * to_x_0[k] - 4.0 * to_x_xx[k] * tke_0 - 2.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxx_xx[k] * tbe_0 * tke_0;

            to_x_x_xx_y[k] = -4.0 * to_x_xy[k] * tke_0 + 4.0 * to_xxx_xy[k] * tbe_0 * tke_0;

            to_x_x_xx_z[k] = -4.0 * to_x_xz[k] * tke_0 + 4.0 * to_xxx_xz[k] * tbe_0 * tke_0;
        }

        // Set up 3-6 components of targeted buffer : DP

        auto to_x_x_xy_x = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 3);

        auto to_x_x_xy_y = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 4);

        auto to_x_x_xy_z = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_x_x_xy_x, to_x_x_xy_y, to_x_x_xy_z, to_xxy_0, to_xxy_xx, to_xxy_xy, to_xxy_xz, to_y_0, to_y_xx, to_y_xy, to_y_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xy_x[k] = to_y_0[k] - 2.0 * to_y_xx[k] * tke_0 - 2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxy_xx[k] * tbe_0 * tke_0;

            to_x_x_xy_y[k] = -2.0 * to_y_xy[k] * tke_0 + 4.0 * to_xxy_xy[k] * tbe_0 * tke_0;

            to_x_x_xy_z[k] = -2.0 * to_y_xz[k] * tke_0 + 4.0 * to_xxy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 6-9 components of targeted buffer : DP

        auto to_x_x_xz_x = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 6);

        auto to_x_x_xz_y = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 7);

        auto to_x_x_xz_z = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_x_x_xz_x, to_x_x_xz_y, to_x_x_xz_z, to_xxz_0, to_xxz_xx, to_xxz_xy, to_xxz_xz, to_z_0, to_z_xx, to_z_xy, to_z_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xz_x[k] = to_z_0[k] - 2.0 * to_z_xx[k] * tke_0 - 2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxz_xx[k] * tbe_0 * tke_0;

            to_x_x_xz_y[k] = -2.0 * to_z_xy[k] * tke_0 + 4.0 * to_xxz_xy[k] * tbe_0 * tke_0;

            to_x_x_xz_z[k] = -2.0 * to_z_xz[k] * tke_0 + 4.0 * to_xxz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 9-12 components of targeted buffer : DP

        auto to_x_x_yy_x = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 9);

        auto to_x_x_yy_y = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 10);

        auto to_x_x_yy_z = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_x_x_yy_x, to_x_x_yy_y, to_x_x_yy_z, to_xyy_0, to_xyy_xx, to_xyy_xy, to_xyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yy_x[k] = -2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyy_xx[k] * tbe_0 * tke_0;

            to_x_x_yy_y[k] = 4.0 * to_xyy_xy[k] * tbe_0 * tke_0;

            to_x_x_yy_z[k] = 4.0 * to_xyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 12-15 components of targeted buffer : DP

        auto to_x_x_yz_x = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 12);

        auto to_x_x_yz_y = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 13);

        auto to_x_x_yz_z = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_x_x_yz_x, to_x_x_yz_y, to_x_x_yz_z, to_xyz_0, to_xyz_xx, to_xyz_xy, to_xyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yz_x[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_xx[k] * tbe_0 * tke_0;

            to_x_x_yz_y[k] = 4.0 * to_xyz_xy[k] * tbe_0 * tke_0;

            to_x_x_yz_z[k] = 4.0 * to_xyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 15-18 components of targeted buffer : DP

        auto to_x_x_zz_x = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 15);

        auto to_x_x_zz_y = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 16);

        auto to_x_x_zz_z = pbuffer.data(idx_op_geom_101_dp + 0 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_x_x_zz_x, to_x_x_zz_y, to_x_x_zz_z, to_xzz_0, to_xzz_xx, to_xzz_xy, to_xzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zz_x[k] = -2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzz_xx[k] * tbe_0 * tke_0;

            to_x_x_zz_y[k] = 4.0 * to_xzz_xy[k] * tbe_0 * tke_0;

            to_x_x_zz_z[k] = 4.0 * to_xzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 18-21 components of targeted buffer : DP

        auto to_x_y_xx_x = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 0);

        auto to_x_y_xx_y = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 1);

        auto to_x_y_xx_z = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_x_0, to_x_xy, to_x_y_xx_x, to_x_y_xx_y, to_x_y_xx_z, to_x_yy, to_x_yz, to_xxx_0, to_xxx_xy, to_xxx_yy, to_xxx_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xx_x[k] = -4.0 * to_x_xy[k] * tke_0 + 4.0 * to_xxx_xy[k] * tbe_0 * tke_0;

            to_x_y_xx_y[k] = 2.0 * to_x_0[k] - 4.0 * to_x_yy[k] * tke_0 - 2.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxx_yy[k] * tbe_0 * tke_0;

            to_x_y_xx_z[k] = -4.0 * to_x_yz[k] * tke_0 + 4.0 * to_xxx_yz[k] * tbe_0 * tke_0;
        }

        // Set up 21-24 components of targeted buffer : DP

        auto to_x_y_xy_x = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 3);

        auto to_x_y_xy_y = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 4);

        auto to_x_y_xy_z = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_x_y_xy_x, to_x_y_xy_y, to_x_y_xy_z, to_xxy_0, to_xxy_xy, to_xxy_yy, to_xxy_yz, to_y_0, to_y_xy, to_y_yy, to_y_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xy_x[k] = -2.0 * to_y_xy[k] * tke_0 + 4.0 * to_xxy_xy[k] * tbe_0 * tke_0;

            to_x_y_xy_y[k] = to_y_0[k] - 2.0 * to_y_yy[k] * tke_0 - 2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxy_yy[k] * tbe_0 * tke_0;

            to_x_y_xy_z[k] = -2.0 * to_y_yz[k] * tke_0 + 4.0 * to_xxy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 24-27 components of targeted buffer : DP

        auto to_x_y_xz_x = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 6);

        auto to_x_y_xz_y = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 7);

        auto to_x_y_xz_z = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_x_y_xz_x, to_x_y_xz_y, to_x_y_xz_z, to_xxz_0, to_xxz_xy, to_xxz_yy, to_xxz_yz, to_z_0, to_z_xy, to_z_yy, to_z_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xz_x[k] = -2.0 * to_z_xy[k] * tke_0 + 4.0 * to_xxz_xy[k] * tbe_0 * tke_0;

            to_x_y_xz_y[k] = to_z_0[k] - 2.0 * to_z_yy[k] * tke_0 - 2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxz_yy[k] * tbe_0 * tke_0;

            to_x_y_xz_z[k] = -2.0 * to_z_yz[k] * tke_0 + 4.0 * to_xxz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 27-30 components of targeted buffer : DP

        auto to_x_y_yy_x = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 9);

        auto to_x_y_yy_y = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 10);

        auto to_x_y_yy_z = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_x_y_yy_x, to_x_y_yy_y, to_x_y_yy_z, to_xyy_0, to_xyy_xy, to_xyy_yy, to_xyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yy_x[k] = 4.0 * to_xyy_xy[k] * tbe_0 * tke_0;

            to_x_y_yy_y[k] = -2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyy_yy[k] * tbe_0 * tke_0;

            to_x_y_yy_z[k] = 4.0 * to_xyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 30-33 components of targeted buffer : DP

        auto to_x_y_yz_x = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 12);

        auto to_x_y_yz_y = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 13);

        auto to_x_y_yz_z = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_x_y_yz_x, to_x_y_yz_y, to_x_y_yz_z, to_xyz_0, to_xyz_xy, to_xyz_yy, to_xyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yz_x[k] = 4.0 * to_xyz_xy[k] * tbe_0 * tke_0;

            to_x_y_yz_y[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_yy[k] * tbe_0 * tke_0;

            to_x_y_yz_z[k] = 4.0 * to_xyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 33-36 components of targeted buffer : DP

        auto to_x_y_zz_x = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 15);

        auto to_x_y_zz_y = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 16);

        auto to_x_y_zz_z = pbuffer.data(idx_op_geom_101_dp + 1 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_x_y_zz_x, to_x_y_zz_y, to_x_y_zz_z, to_xzz_0, to_xzz_xy, to_xzz_yy, to_xzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zz_x[k] = 4.0 * to_xzz_xy[k] * tbe_0 * tke_0;

            to_x_y_zz_y[k] = -2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzz_yy[k] * tbe_0 * tke_0;

            to_x_y_zz_z[k] = 4.0 * to_xzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 36-39 components of targeted buffer : DP

        auto to_x_z_xx_x = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 0);

        auto to_x_z_xx_y = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 1);

        auto to_x_z_xx_z = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_x_0, to_x_xz, to_x_yz, to_x_z_xx_x, to_x_z_xx_y, to_x_z_xx_z, to_x_zz, to_xxx_0, to_xxx_xz, to_xxx_yz, to_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xx_x[k] = -4.0 * to_x_xz[k] * tke_0 + 4.0 * to_xxx_xz[k] * tbe_0 * tke_0;

            to_x_z_xx_y[k] = -4.0 * to_x_yz[k] * tke_0 + 4.0 * to_xxx_yz[k] * tbe_0 * tke_0;

            to_x_z_xx_z[k] = 2.0 * to_x_0[k] - 4.0 * to_x_zz[k] * tke_0 - 2.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxx_zz[k] * tbe_0 * tke_0;
        }

        // Set up 39-42 components of targeted buffer : DP

        auto to_x_z_xy_x = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 3);

        auto to_x_z_xy_y = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 4);

        auto to_x_z_xy_z = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_x_z_xy_x, to_x_z_xy_y, to_x_z_xy_z, to_xxy_0, to_xxy_xz, to_xxy_yz, to_xxy_zz, to_y_0, to_y_xz, to_y_yz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xy_x[k] = -2.0 * to_y_xz[k] * tke_0 + 4.0 * to_xxy_xz[k] * tbe_0 * tke_0;

            to_x_z_xy_y[k] = -2.0 * to_y_yz[k] * tke_0 + 4.0 * to_xxy_yz[k] * tbe_0 * tke_0;

            to_x_z_xy_z[k] = to_y_0[k] - 2.0 * to_y_zz[k] * tke_0 - 2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 42-45 components of targeted buffer : DP

        auto to_x_z_xz_x = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 6);

        auto to_x_z_xz_y = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 7);

        auto to_x_z_xz_z = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_x_z_xz_x, to_x_z_xz_y, to_x_z_xz_z, to_xxz_0, to_xxz_xz, to_xxz_yz, to_xxz_zz, to_z_0, to_z_xz, to_z_yz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xz_x[k] = -2.0 * to_z_xz[k] * tke_0 + 4.0 * to_xxz_xz[k] * tbe_0 * tke_0;

            to_x_z_xz_y[k] = -2.0 * to_z_yz[k] * tke_0 + 4.0 * to_xxz_yz[k] * tbe_0 * tke_0;

            to_x_z_xz_z[k] = to_z_0[k] - 2.0 * to_z_zz[k] * tke_0 - 2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 45-48 components of targeted buffer : DP

        auto to_x_z_yy_x = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 9);

        auto to_x_z_yy_y = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 10);

        auto to_x_z_yy_z = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_x_z_yy_x, to_x_z_yy_y, to_x_z_yy_z, to_xyy_0, to_xyy_xz, to_xyy_yz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yy_x[k] = 4.0 * to_xyy_xz[k] * tbe_0 * tke_0;

            to_x_z_yy_y[k] = 4.0 * to_xyy_yz[k] * tbe_0 * tke_0;

            to_x_z_yy_z[k] = -2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 48-51 components of targeted buffer : DP

        auto to_x_z_yz_x = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 12);

        auto to_x_z_yz_y = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 13);

        auto to_x_z_yz_z = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_x_z_yz_x, to_x_z_yz_y, to_x_z_yz_z, to_xyz_0, to_xyz_xz, to_xyz_yz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yz_x[k] = 4.0 * to_xyz_xz[k] * tbe_0 * tke_0;

            to_x_z_yz_y[k] = 4.0 * to_xyz_yz[k] * tbe_0 * tke_0;

            to_x_z_yz_z[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 51-54 components of targeted buffer : DP

        auto to_x_z_zz_x = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 15);

        auto to_x_z_zz_y = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 16);

        auto to_x_z_zz_z = pbuffer.data(idx_op_geom_101_dp + 2 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_x_z_zz_x, to_x_z_zz_y, to_x_z_zz_z, to_xzz_0, to_xzz_xz, to_xzz_yz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zz_x[k] = 4.0 * to_xzz_xz[k] * tbe_0 * tke_0;

            to_x_z_zz_y[k] = 4.0 * to_xzz_yz[k] * tbe_0 * tke_0;

            to_x_z_zz_z[k] = -2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 54-57 components of targeted buffer : DP

        auto to_y_x_xx_x = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 0);

        auto to_y_x_xx_y = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 1);

        auto to_y_x_xx_z = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_xxy_0, to_xxy_xx, to_xxy_xy, to_xxy_xz, to_y_x_xx_x, to_y_x_xx_y, to_y_x_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xx_x[k] = -2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxy_xx[k] * tbe_0 * tke_0;

            to_y_x_xx_y[k] = 4.0 * to_xxy_xy[k] * tbe_0 * tke_0;

            to_y_x_xx_z[k] = 4.0 * to_xxy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 57-60 components of targeted buffer : DP

        auto to_y_x_xy_x = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 3);

        auto to_y_x_xy_y = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 4);

        auto to_y_x_xy_z = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_x_0, to_x_xx, to_x_xy, to_x_xz, to_xyy_0, to_xyy_xx, to_xyy_xy, to_xyy_xz, to_y_x_xy_x, to_y_x_xy_y, to_y_x_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xy_x[k] = to_x_0[k] - 2.0 * to_x_xx[k] * tke_0 - 2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyy_xx[k] * tbe_0 * tke_0;

            to_y_x_xy_y[k] = -2.0 * to_x_xy[k] * tke_0 + 4.0 * to_xyy_xy[k] * tbe_0 * tke_0;

            to_y_x_xy_z[k] = -2.0 * to_x_xz[k] * tke_0 + 4.0 * to_xyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 60-63 components of targeted buffer : DP

        auto to_y_x_xz_x = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 6);

        auto to_y_x_xz_y = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 7);

        auto to_y_x_xz_z = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_xyz_0, to_xyz_xx, to_xyz_xy, to_xyz_xz, to_y_x_xz_x, to_y_x_xz_y, to_y_x_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xz_x[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_xx[k] * tbe_0 * tke_0;

            to_y_x_xz_y[k] = 4.0 * to_xyz_xy[k] * tbe_0 * tke_0;

            to_y_x_xz_z[k] = 4.0 * to_xyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 63-66 components of targeted buffer : DP

        auto to_y_x_yy_x = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 9);

        auto to_y_x_yy_y = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 10);

        auto to_y_x_yy_z = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_y_0, to_y_x_yy_x, to_y_x_yy_y, to_y_x_yy_z, to_y_xx, to_y_xy, to_y_xz, to_yyy_0, to_yyy_xx, to_yyy_xy, to_yyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yy_x[k] = 2.0 * to_y_0[k] - 4.0 * to_y_xx[k] * tke_0 - 2.0 * to_yyy_0[k] * tbe_0 + 4.0 * to_yyy_xx[k] * tbe_0 * tke_0;

            to_y_x_yy_y[k] = -4.0 * to_y_xy[k] * tke_0 + 4.0 * to_yyy_xy[k] * tbe_0 * tke_0;

            to_y_x_yy_z[k] = -4.0 * to_y_xz[k] * tke_0 + 4.0 * to_yyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 66-69 components of targeted buffer : DP

        auto to_y_x_yz_x = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 12);

        auto to_y_x_yz_y = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 13);

        auto to_y_x_yz_z = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_y_x_yz_x, to_y_x_yz_y, to_y_x_yz_z, to_yyz_0, to_yyz_xx, to_yyz_xy, to_yyz_xz, to_z_0, to_z_xx, to_z_xy, to_z_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yz_x[k] = to_z_0[k] - 2.0 * to_z_xx[k] * tke_0 - 2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyz_xx[k] * tbe_0 * tke_0;

            to_y_x_yz_y[k] = -2.0 * to_z_xy[k] * tke_0 + 4.0 * to_yyz_xy[k] * tbe_0 * tke_0;

            to_y_x_yz_z[k] = -2.0 * to_z_xz[k] * tke_0 + 4.0 * to_yyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 69-72 components of targeted buffer : DP

        auto to_y_x_zz_x = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 15);

        auto to_y_x_zz_y = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 16);

        auto to_y_x_zz_z = pbuffer.data(idx_op_geom_101_dp + 3 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_y_x_zz_x, to_y_x_zz_y, to_y_x_zz_z, to_yzz_0, to_yzz_xx, to_yzz_xy, to_yzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zz_x[k] = -2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzz_xx[k] * tbe_0 * tke_0;

            to_y_x_zz_y[k] = 4.0 * to_yzz_xy[k] * tbe_0 * tke_0;

            to_y_x_zz_z[k] = 4.0 * to_yzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 72-75 components of targeted buffer : DP

        auto to_y_y_xx_x = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 0);

        auto to_y_y_xx_y = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 1);

        auto to_y_y_xx_z = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_xxy_0, to_xxy_xy, to_xxy_yy, to_xxy_yz, to_y_y_xx_x, to_y_y_xx_y, to_y_y_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xx_x[k] = 4.0 * to_xxy_xy[k] * tbe_0 * tke_0;

            to_y_y_xx_y[k] = -2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxy_yy[k] * tbe_0 * tke_0;

            to_y_y_xx_z[k] = 4.0 * to_xxy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 75-78 components of targeted buffer : DP

        auto to_y_y_xy_x = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 3);

        auto to_y_y_xy_y = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 4);

        auto to_y_y_xy_z = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_x_0, to_x_xy, to_x_yy, to_x_yz, to_xyy_0, to_xyy_xy, to_xyy_yy, to_xyy_yz, to_y_y_xy_x, to_y_y_xy_y, to_y_y_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xy_x[k] = -2.0 * to_x_xy[k] * tke_0 + 4.0 * to_xyy_xy[k] * tbe_0 * tke_0;

            to_y_y_xy_y[k] = to_x_0[k] - 2.0 * to_x_yy[k] * tke_0 - 2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyy_yy[k] * tbe_0 * tke_0;

            to_y_y_xy_z[k] = -2.0 * to_x_yz[k] * tke_0 + 4.0 * to_xyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 78-81 components of targeted buffer : DP

        auto to_y_y_xz_x = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 6);

        auto to_y_y_xz_y = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 7);

        auto to_y_y_xz_z = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_xyz_0, to_xyz_xy, to_xyz_yy, to_xyz_yz, to_y_y_xz_x, to_y_y_xz_y, to_y_y_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xz_x[k] = 4.0 * to_xyz_xy[k] * tbe_0 * tke_0;

            to_y_y_xz_y[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_yy[k] * tbe_0 * tke_0;

            to_y_y_xz_z[k] = 4.0 * to_xyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 81-84 components of targeted buffer : DP

        auto to_y_y_yy_x = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 9);

        auto to_y_y_yy_y = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 10);

        auto to_y_y_yy_z = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_y_0, to_y_xy, to_y_y_yy_x, to_y_y_yy_y, to_y_y_yy_z, to_y_yy, to_y_yz, to_yyy_0, to_yyy_xy, to_yyy_yy, to_yyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yy_x[k] = -4.0 * to_y_xy[k] * tke_0 + 4.0 * to_yyy_xy[k] * tbe_0 * tke_0;

            to_y_y_yy_y[k] = 2.0 * to_y_0[k] - 4.0 * to_y_yy[k] * tke_0 - 2.0 * to_yyy_0[k] * tbe_0 + 4.0 * to_yyy_yy[k] * tbe_0 * tke_0;

            to_y_y_yy_z[k] = -4.0 * to_y_yz[k] * tke_0 + 4.0 * to_yyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 84-87 components of targeted buffer : DP

        auto to_y_y_yz_x = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 12);

        auto to_y_y_yz_y = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 13);

        auto to_y_y_yz_z = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_y_y_yz_x, to_y_y_yz_y, to_y_y_yz_z, to_yyz_0, to_yyz_xy, to_yyz_yy, to_yyz_yz, to_z_0, to_z_xy, to_z_yy, to_z_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yz_x[k] = -2.0 * to_z_xy[k] * tke_0 + 4.0 * to_yyz_xy[k] * tbe_0 * tke_0;

            to_y_y_yz_y[k] = to_z_0[k] - 2.0 * to_z_yy[k] * tke_0 - 2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyz_yy[k] * tbe_0 * tke_0;

            to_y_y_yz_z[k] = -2.0 * to_z_yz[k] * tke_0 + 4.0 * to_yyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 87-90 components of targeted buffer : DP

        auto to_y_y_zz_x = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 15);

        auto to_y_y_zz_y = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 16);

        auto to_y_y_zz_z = pbuffer.data(idx_op_geom_101_dp + 4 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_y_y_zz_x, to_y_y_zz_y, to_y_y_zz_z, to_yzz_0, to_yzz_xy, to_yzz_yy, to_yzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zz_x[k] = 4.0 * to_yzz_xy[k] * tbe_0 * tke_0;

            to_y_y_zz_y[k] = -2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzz_yy[k] * tbe_0 * tke_0;

            to_y_y_zz_z[k] = 4.0 * to_yzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 90-93 components of targeted buffer : DP

        auto to_y_z_xx_x = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 0);

        auto to_y_z_xx_y = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 1);

        auto to_y_z_xx_z = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_xxy_0, to_xxy_xz, to_xxy_yz, to_xxy_zz, to_y_z_xx_x, to_y_z_xx_y, to_y_z_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xx_x[k] = 4.0 * to_xxy_xz[k] * tbe_0 * tke_0;

            to_y_z_xx_y[k] = 4.0 * to_xxy_yz[k] * tbe_0 * tke_0;

            to_y_z_xx_z[k] = -2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 93-96 components of targeted buffer : DP

        auto to_y_z_xy_x = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 3);

        auto to_y_z_xy_y = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 4);

        auto to_y_z_xy_z = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_x_0, to_x_xz, to_x_yz, to_x_zz, to_xyy_0, to_xyy_xz, to_xyy_yz, to_xyy_zz, to_y_z_xy_x, to_y_z_xy_y, to_y_z_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xy_x[k] = -2.0 * to_x_xz[k] * tke_0 + 4.0 * to_xyy_xz[k] * tbe_0 * tke_0;

            to_y_z_xy_y[k] = -2.0 * to_x_yz[k] * tke_0 + 4.0 * to_xyy_yz[k] * tbe_0 * tke_0;

            to_y_z_xy_z[k] = to_x_0[k] - 2.0 * to_x_zz[k] * tke_0 - 2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 96-99 components of targeted buffer : DP

        auto to_y_z_xz_x = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 6);

        auto to_y_z_xz_y = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 7);

        auto to_y_z_xz_z = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_xyz_0, to_xyz_xz, to_xyz_yz, to_xyz_zz, to_y_z_xz_x, to_y_z_xz_y, to_y_z_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xz_x[k] = 4.0 * to_xyz_xz[k] * tbe_0 * tke_0;

            to_y_z_xz_y[k] = 4.0 * to_xyz_yz[k] * tbe_0 * tke_0;

            to_y_z_xz_z[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 99-102 components of targeted buffer : DP

        auto to_y_z_yy_x = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 9);

        auto to_y_z_yy_y = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 10);

        auto to_y_z_yy_z = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_y_0, to_y_xz, to_y_yz, to_y_z_yy_x, to_y_z_yy_y, to_y_z_yy_z, to_y_zz, to_yyy_0, to_yyy_xz, to_yyy_yz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yy_x[k] = -4.0 * to_y_xz[k] * tke_0 + 4.0 * to_yyy_xz[k] * tbe_0 * tke_0;

            to_y_z_yy_y[k] = -4.0 * to_y_yz[k] * tke_0 + 4.0 * to_yyy_yz[k] * tbe_0 * tke_0;

            to_y_z_yy_z[k] = 2.0 * to_y_0[k] - 4.0 * to_y_zz[k] * tke_0 - 2.0 * to_yyy_0[k] * tbe_0 + 4.0 * to_yyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 102-105 components of targeted buffer : DP

        auto to_y_z_yz_x = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 12);

        auto to_y_z_yz_y = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 13);

        auto to_y_z_yz_z = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_y_z_yz_x, to_y_z_yz_y, to_y_z_yz_z, to_yyz_0, to_yyz_xz, to_yyz_yz, to_yyz_zz, to_z_0, to_z_xz, to_z_yz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yz_x[k] = -2.0 * to_z_xz[k] * tke_0 + 4.0 * to_yyz_xz[k] * tbe_0 * tke_0;

            to_y_z_yz_y[k] = -2.0 * to_z_yz[k] * tke_0 + 4.0 * to_yyz_yz[k] * tbe_0 * tke_0;

            to_y_z_yz_z[k] = to_z_0[k] - 2.0 * to_z_zz[k] * tke_0 - 2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 105-108 components of targeted buffer : DP

        auto to_y_z_zz_x = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 15);

        auto to_y_z_zz_y = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 16);

        auto to_y_z_zz_z = pbuffer.data(idx_op_geom_101_dp + 5 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_y_z_zz_x, to_y_z_zz_y, to_y_z_zz_z, to_yzz_0, to_yzz_xz, to_yzz_yz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zz_x[k] = 4.0 * to_yzz_xz[k] * tbe_0 * tke_0;

            to_y_z_zz_y[k] = 4.0 * to_yzz_yz[k] * tbe_0 * tke_0;

            to_y_z_zz_z[k] = -2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 108-111 components of targeted buffer : DP

        auto to_z_x_xx_x = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 0);

        auto to_z_x_xx_y = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 1);

        auto to_z_x_xx_z = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_xxz_0, to_xxz_xx, to_xxz_xy, to_xxz_xz, to_z_x_xx_x, to_z_x_xx_y, to_z_x_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xx_x[k] = -2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxz_xx[k] * tbe_0 * tke_0;

            to_z_x_xx_y[k] = 4.0 * to_xxz_xy[k] * tbe_0 * tke_0;

            to_z_x_xx_z[k] = 4.0 * to_xxz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 111-114 components of targeted buffer : DP

        auto to_z_x_xy_x = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 3);

        auto to_z_x_xy_y = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 4);

        auto to_z_x_xy_z = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_xyz_0, to_xyz_xx, to_xyz_xy, to_xyz_xz, to_z_x_xy_x, to_z_x_xy_y, to_z_x_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xy_x[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_xx[k] * tbe_0 * tke_0;

            to_z_x_xy_y[k] = 4.0 * to_xyz_xy[k] * tbe_0 * tke_0;

            to_z_x_xy_z[k] = 4.0 * to_xyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 114-117 components of targeted buffer : DP

        auto to_z_x_xz_x = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 6);

        auto to_z_x_xz_y = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 7);

        auto to_z_x_xz_z = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_x_0, to_x_xx, to_x_xy, to_x_xz, to_xzz_0, to_xzz_xx, to_xzz_xy, to_xzz_xz, to_z_x_xz_x, to_z_x_xz_y, to_z_x_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xz_x[k] = to_x_0[k] - 2.0 * to_x_xx[k] * tke_0 - 2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xz_y[k] = -2.0 * to_x_xy[k] * tke_0 + 4.0 * to_xzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xz_z[k] = -2.0 * to_x_xz[k] * tke_0 + 4.0 * to_xzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 117-120 components of targeted buffer : DP

        auto to_z_x_yy_x = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 9);

        auto to_z_x_yy_y = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 10);

        auto to_z_x_yy_z = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_yyz_0, to_yyz_xx, to_yyz_xy, to_yyz_xz, to_z_x_yy_x, to_z_x_yy_y, to_z_x_yy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yy_x[k] = -2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyz_xx[k] * tbe_0 * tke_0;

            to_z_x_yy_y[k] = 4.0 * to_yyz_xy[k] * tbe_0 * tke_0;

            to_z_x_yy_z[k] = 4.0 * to_yyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 120-123 components of targeted buffer : DP

        auto to_z_x_yz_x = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 12);

        auto to_z_x_yz_y = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 13);

        auto to_z_x_yz_z = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_y_0, to_y_xx, to_y_xy, to_y_xz, to_yzz_0, to_yzz_xx, to_yzz_xy, to_yzz_xz, to_z_x_yz_x, to_z_x_yz_y, to_z_x_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yz_x[k] = to_y_0[k] - 2.0 * to_y_xx[k] * tke_0 - 2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzz_xx[k] * tbe_0 * tke_0;

            to_z_x_yz_y[k] = -2.0 * to_y_xy[k] * tke_0 + 4.0 * to_yzz_xy[k] * tbe_0 * tke_0;

            to_z_x_yz_z[k] = -2.0 * to_y_xz[k] * tke_0 + 4.0 * to_yzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 123-126 components of targeted buffer : DP

        auto to_z_x_zz_x = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 15);

        auto to_z_x_zz_y = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 16);

        auto to_z_x_zz_z = pbuffer.data(idx_op_geom_101_dp + 6 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_z_0, to_z_x_zz_x, to_z_x_zz_y, to_z_x_zz_z, to_z_xx, to_z_xy, to_z_xz, to_zzz_0, to_zzz_xx, to_zzz_xy, to_zzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zz_x[k] = 2.0 * to_z_0[k] - 4.0 * to_z_xx[k] * tke_0 - 2.0 * to_zzz_0[k] * tbe_0 + 4.0 * to_zzz_xx[k] * tbe_0 * tke_0;

            to_z_x_zz_y[k] = -4.0 * to_z_xy[k] * tke_0 + 4.0 * to_zzz_xy[k] * tbe_0 * tke_0;

            to_z_x_zz_z[k] = -4.0 * to_z_xz[k] * tke_0 + 4.0 * to_zzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 126-129 components of targeted buffer : DP

        auto to_z_y_xx_x = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 0);

        auto to_z_y_xx_y = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 1);

        auto to_z_y_xx_z = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_xxz_0, to_xxz_xy, to_xxz_yy, to_xxz_yz, to_z_y_xx_x, to_z_y_xx_y, to_z_y_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xx_x[k] = 4.0 * to_xxz_xy[k] * tbe_0 * tke_0;

            to_z_y_xx_y[k] = -2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxz_yy[k] * tbe_0 * tke_0;

            to_z_y_xx_z[k] = 4.0 * to_xxz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 129-132 components of targeted buffer : DP

        auto to_z_y_xy_x = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 3);

        auto to_z_y_xy_y = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 4);

        auto to_z_y_xy_z = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_xyz_0, to_xyz_xy, to_xyz_yy, to_xyz_yz, to_z_y_xy_x, to_z_y_xy_y, to_z_y_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xy_x[k] = 4.0 * to_xyz_xy[k] * tbe_0 * tke_0;

            to_z_y_xy_y[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_yy[k] * tbe_0 * tke_0;

            to_z_y_xy_z[k] = 4.0 * to_xyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 132-135 components of targeted buffer : DP

        auto to_z_y_xz_x = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 6);

        auto to_z_y_xz_y = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 7);

        auto to_z_y_xz_z = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_x_0, to_x_xy, to_x_yy, to_x_yz, to_xzz_0, to_xzz_xy, to_xzz_yy, to_xzz_yz, to_z_y_xz_x, to_z_y_xz_y, to_z_y_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xz_x[k] = -2.0 * to_x_xy[k] * tke_0 + 4.0 * to_xzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xz_y[k] = to_x_0[k] - 2.0 * to_x_yy[k] * tke_0 - 2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xz_z[k] = -2.0 * to_x_yz[k] * tke_0 + 4.0 * to_xzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 135-138 components of targeted buffer : DP

        auto to_z_y_yy_x = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 9);

        auto to_z_y_yy_y = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 10);

        auto to_z_y_yy_z = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_yyz_0, to_yyz_xy, to_yyz_yy, to_yyz_yz, to_z_y_yy_x, to_z_y_yy_y, to_z_y_yy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yy_x[k] = 4.0 * to_yyz_xy[k] * tbe_0 * tke_0;

            to_z_y_yy_y[k] = -2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyz_yy[k] * tbe_0 * tke_0;

            to_z_y_yy_z[k] = 4.0 * to_yyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 138-141 components of targeted buffer : DP

        auto to_z_y_yz_x = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 12);

        auto to_z_y_yz_y = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 13);

        auto to_z_y_yz_z = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_y_0, to_y_xy, to_y_yy, to_y_yz, to_yzz_0, to_yzz_xy, to_yzz_yy, to_yzz_yz, to_z_y_yz_x, to_z_y_yz_y, to_z_y_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yz_x[k] = -2.0 * to_y_xy[k] * tke_0 + 4.0 * to_yzz_xy[k] * tbe_0 * tke_0;

            to_z_y_yz_y[k] = to_y_0[k] - 2.0 * to_y_yy[k] * tke_0 - 2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzz_yy[k] * tbe_0 * tke_0;

            to_z_y_yz_z[k] = -2.0 * to_y_yz[k] * tke_0 + 4.0 * to_yzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 141-144 components of targeted buffer : DP

        auto to_z_y_zz_x = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 15);

        auto to_z_y_zz_y = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 16);

        auto to_z_y_zz_z = pbuffer.data(idx_op_geom_101_dp + 7 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_z_0, to_z_xy, to_z_y_zz_x, to_z_y_zz_y, to_z_y_zz_z, to_z_yy, to_z_yz, to_zzz_0, to_zzz_xy, to_zzz_yy, to_zzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zz_x[k] = -4.0 * to_z_xy[k] * tke_0 + 4.0 * to_zzz_xy[k] * tbe_0 * tke_0;

            to_z_y_zz_y[k] = 2.0 * to_z_0[k] - 4.0 * to_z_yy[k] * tke_0 - 2.0 * to_zzz_0[k] * tbe_0 + 4.0 * to_zzz_yy[k] * tbe_0 * tke_0;

            to_z_y_zz_z[k] = -4.0 * to_z_yz[k] * tke_0 + 4.0 * to_zzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 144-147 components of targeted buffer : DP

        auto to_z_z_xx_x = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 0);

        auto to_z_z_xx_y = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 1);

        auto to_z_z_xx_z = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_xxz_0, to_xxz_xz, to_xxz_yz, to_xxz_zz, to_z_z_xx_x, to_z_z_xx_y, to_z_z_xx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xx_x[k] = 4.0 * to_xxz_xz[k] * tbe_0 * tke_0;

            to_z_z_xx_y[k] = 4.0 * to_xxz_yz[k] * tbe_0 * tke_0;

            to_z_z_xx_z[k] = -2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 147-150 components of targeted buffer : DP

        auto to_z_z_xy_x = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 3);

        auto to_z_z_xy_y = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 4);

        auto to_z_z_xy_z = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_xyz_0, to_xyz_xz, to_xyz_yz, to_xyz_zz, to_z_z_xy_x, to_z_z_xy_y, to_z_z_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xy_x[k] = 4.0 * to_xyz_xz[k] * tbe_0 * tke_0;

            to_z_z_xy_y[k] = 4.0 * to_xyz_yz[k] * tbe_0 * tke_0;

            to_z_z_xy_z[k] = -2.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 150-153 components of targeted buffer : DP

        auto to_z_z_xz_x = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 6);

        auto to_z_z_xz_y = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 7);

        auto to_z_z_xz_z = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_x_0, to_x_xz, to_x_yz, to_x_zz, to_xzz_0, to_xzz_xz, to_xzz_yz, to_xzz_zz, to_z_z_xz_x, to_z_z_xz_y, to_z_z_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xz_x[k] = -2.0 * to_x_xz[k] * tke_0 + 4.0 * to_xzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xz_y[k] = -2.0 * to_x_yz[k] * tke_0 + 4.0 * to_xzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xz_z[k] = to_x_0[k] - 2.0 * to_x_zz[k] * tke_0 - 2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 153-156 components of targeted buffer : DP

        auto to_z_z_yy_x = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 9);

        auto to_z_z_yy_y = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 10);

        auto to_z_z_yy_z = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_yyz_0, to_yyz_xz, to_yyz_yz, to_yyz_zz, to_z_z_yy_x, to_z_z_yy_y, to_z_z_yy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yy_x[k] = 4.0 * to_yyz_xz[k] * tbe_0 * tke_0;

            to_z_z_yy_y[k] = 4.0 * to_yyz_yz[k] * tbe_0 * tke_0;

            to_z_z_yy_z[k] = -2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 156-159 components of targeted buffer : DP

        auto to_z_z_yz_x = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 12);

        auto to_z_z_yz_y = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 13);

        auto to_z_z_yz_z = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_y_0, to_y_xz, to_y_yz, to_y_zz, to_yzz_0, to_yzz_xz, to_yzz_yz, to_yzz_zz, to_z_z_yz_x, to_z_z_yz_y, to_z_z_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yz_x[k] = -2.0 * to_y_xz[k] * tke_0 + 4.0 * to_yzz_xz[k] * tbe_0 * tke_0;

            to_z_z_yz_y[k] = -2.0 * to_y_yz[k] * tke_0 + 4.0 * to_yzz_yz[k] * tbe_0 * tke_0;

            to_z_z_yz_z[k] = to_y_0[k] - 2.0 * to_y_zz[k] * tke_0 - 2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 159-162 components of targeted buffer : DP

        auto to_z_z_zz_x = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 15);

        auto to_z_z_zz_y = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 16);

        auto to_z_z_zz_z = pbuffer.data(idx_op_geom_101_dp + 8 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_z_0, to_z_xz, to_z_yz, to_z_z_zz_x, to_z_z_zz_y, to_z_z_zz_z, to_z_zz, to_zzz_0, to_zzz_xz, to_zzz_yz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zz_x[k] = -4.0 * to_z_xz[k] * tke_0 + 4.0 * to_zzz_xz[k] * tbe_0 * tke_0;

            to_z_z_zz_y[k] = -4.0 * to_z_yz[k] * tke_0 + 4.0 * to_zzz_yz[k] * tbe_0 * tke_0;

            to_z_z_zz_z[k] = 2.0 * to_z_0[k] - 4.0 * to_z_zz[k] * tke_0 - 2.0 * to_zzz_0[k] * tbe_0 + 4.0 * to_zzz_zz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

