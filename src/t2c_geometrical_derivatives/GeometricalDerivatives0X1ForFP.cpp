#include "GeometricalDerivatives0X1ForFP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_fp(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_fp,
                       const int idx_op_fs,
                       const int idx_op_fd,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
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

        // Set up 0-3 components of targeted buffer : FP

        auto to_0_x_xxx_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 0);

        auto to_0_x_xxx_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 1);

        auto to_0_x_xxx_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_0_x_xxx_x, to_0_x_xxx_y, to_0_x_xxx_z, to_xxx_0, to_xxx_xx, to_xxx_xy, to_xxx_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxx_x[k] = -to_xxx_0[k] + 2.0 * to_xxx_xx[k] * tke_0;

            to_0_x_xxx_y[k] = 2.0 * to_xxx_xy[k] * tke_0;

            to_0_x_xxx_z[k] = 2.0 * to_xxx_xz[k] * tke_0;
        }

        // Set up 3-6 components of targeted buffer : FP

        auto to_0_x_xxy_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 3);

        auto to_0_x_xxy_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 4);

        auto to_0_x_xxy_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_0_x_xxy_x, to_0_x_xxy_y, to_0_x_xxy_z, to_xxy_0, to_xxy_xx, to_xxy_xy, to_xxy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxy_x[k] = -to_xxy_0[k] + 2.0 * to_xxy_xx[k] * tke_0;

            to_0_x_xxy_y[k] = 2.0 * to_xxy_xy[k] * tke_0;

            to_0_x_xxy_z[k] = 2.0 * to_xxy_xz[k] * tke_0;
        }

        // Set up 6-9 components of targeted buffer : FP

        auto to_0_x_xxz_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 6);

        auto to_0_x_xxz_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 7);

        auto to_0_x_xxz_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_0_x_xxz_x, to_0_x_xxz_y, to_0_x_xxz_z, to_xxz_0, to_xxz_xx, to_xxz_xy, to_xxz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxz_x[k] = -to_xxz_0[k] + 2.0 * to_xxz_xx[k] * tke_0;

            to_0_x_xxz_y[k] = 2.0 * to_xxz_xy[k] * tke_0;

            to_0_x_xxz_z[k] = 2.0 * to_xxz_xz[k] * tke_0;
        }

        // Set up 9-12 components of targeted buffer : FP

        auto to_0_x_xyy_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 9);

        auto to_0_x_xyy_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 10);

        auto to_0_x_xyy_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_0_x_xyy_x, to_0_x_xyy_y, to_0_x_xyy_z, to_xyy_0, to_xyy_xx, to_xyy_xy, to_xyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyy_x[k] = -to_xyy_0[k] + 2.0 * to_xyy_xx[k] * tke_0;

            to_0_x_xyy_y[k] = 2.0 * to_xyy_xy[k] * tke_0;

            to_0_x_xyy_z[k] = 2.0 * to_xyy_xz[k] * tke_0;
        }

        // Set up 12-15 components of targeted buffer : FP

        auto to_0_x_xyz_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 12);

        auto to_0_x_xyz_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 13);

        auto to_0_x_xyz_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_0_x_xyz_x, to_0_x_xyz_y, to_0_x_xyz_z, to_xyz_0, to_xyz_xx, to_xyz_xy, to_xyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyz_x[k] = -to_xyz_0[k] + 2.0 * to_xyz_xx[k] * tke_0;

            to_0_x_xyz_y[k] = 2.0 * to_xyz_xy[k] * tke_0;

            to_0_x_xyz_z[k] = 2.0 * to_xyz_xz[k] * tke_0;
        }

        // Set up 15-18 components of targeted buffer : FP

        auto to_0_x_xzz_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 15);

        auto to_0_x_xzz_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 16);

        auto to_0_x_xzz_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_0_x_xzz_x, to_0_x_xzz_y, to_0_x_xzz_z, to_xzz_0, to_xzz_xx, to_xzz_xy, to_xzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzz_x[k] = -to_xzz_0[k] + 2.0 * to_xzz_xx[k] * tke_0;

            to_0_x_xzz_y[k] = 2.0 * to_xzz_xy[k] * tke_0;

            to_0_x_xzz_z[k] = 2.0 * to_xzz_xz[k] * tke_0;
        }

        // Set up 18-21 components of targeted buffer : FP

        auto to_0_x_yyy_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 18);

        auto to_0_x_yyy_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 19);

        auto to_0_x_yyy_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_0_x_yyy_x, to_0_x_yyy_y, to_0_x_yyy_z, to_yyy_0, to_yyy_xx, to_yyy_xy, to_yyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyy_x[k] = -to_yyy_0[k] + 2.0 * to_yyy_xx[k] * tke_0;

            to_0_x_yyy_y[k] = 2.0 * to_yyy_xy[k] * tke_0;

            to_0_x_yyy_z[k] = 2.0 * to_yyy_xz[k] * tke_0;
        }

        // Set up 21-24 components of targeted buffer : FP

        auto to_0_x_yyz_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 21);

        auto to_0_x_yyz_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 22);

        auto to_0_x_yyz_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_0_x_yyz_x, to_0_x_yyz_y, to_0_x_yyz_z, to_yyz_0, to_yyz_xx, to_yyz_xy, to_yyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyz_x[k] = -to_yyz_0[k] + 2.0 * to_yyz_xx[k] * tke_0;

            to_0_x_yyz_y[k] = 2.0 * to_yyz_xy[k] * tke_0;

            to_0_x_yyz_z[k] = 2.0 * to_yyz_xz[k] * tke_0;
        }

        // Set up 24-27 components of targeted buffer : FP

        auto to_0_x_yzz_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 24);

        auto to_0_x_yzz_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 25);

        auto to_0_x_yzz_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_0_x_yzz_x, to_0_x_yzz_y, to_0_x_yzz_z, to_yzz_0, to_yzz_xx, to_yzz_xy, to_yzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzz_x[k] = -to_yzz_0[k] + 2.0 * to_yzz_xx[k] * tke_0;

            to_0_x_yzz_y[k] = 2.0 * to_yzz_xy[k] * tke_0;

            to_0_x_yzz_z[k] = 2.0 * to_yzz_xz[k] * tke_0;
        }

        // Set up 27-30 components of targeted buffer : FP

        auto to_0_x_zzz_x = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 27);

        auto to_0_x_zzz_y = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 28);

        auto to_0_x_zzz_z = pbuffer.data(idx_op_geom_001_fp + 0 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_x_zzz_x, to_0_x_zzz_y, to_0_x_zzz_z, to_zzz_0, to_zzz_xx, to_zzz_xy, to_zzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzz_x[k] = -to_zzz_0[k] + 2.0 * to_zzz_xx[k] * tke_0;

            to_0_x_zzz_y[k] = 2.0 * to_zzz_xy[k] * tke_0;

            to_0_x_zzz_z[k] = 2.0 * to_zzz_xz[k] * tke_0;
        }

        // Set up 30-33 components of targeted buffer : FP

        auto to_0_y_xxx_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 0);

        auto to_0_y_xxx_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 1);

        auto to_0_y_xxx_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_0_y_xxx_x, to_0_y_xxx_y, to_0_y_xxx_z, to_xxx_0, to_xxx_xy, to_xxx_yy, to_xxx_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxx_x[k] = 2.0 * to_xxx_xy[k] * tke_0;

            to_0_y_xxx_y[k] = -to_xxx_0[k] + 2.0 * to_xxx_yy[k] * tke_0;

            to_0_y_xxx_z[k] = 2.0 * to_xxx_yz[k] * tke_0;
        }

        // Set up 33-36 components of targeted buffer : FP

        auto to_0_y_xxy_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 3);

        auto to_0_y_xxy_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 4);

        auto to_0_y_xxy_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_0_y_xxy_x, to_0_y_xxy_y, to_0_y_xxy_z, to_xxy_0, to_xxy_xy, to_xxy_yy, to_xxy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxy_x[k] = 2.0 * to_xxy_xy[k] * tke_0;

            to_0_y_xxy_y[k] = -to_xxy_0[k] + 2.0 * to_xxy_yy[k] * tke_0;

            to_0_y_xxy_z[k] = 2.0 * to_xxy_yz[k] * tke_0;
        }

        // Set up 36-39 components of targeted buffer : FP

        auto to_0_y_xxz_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 6);

        auto to_0_y_xxz_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 7);

        auto to_0_y_xxz_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_0_y_xxz_x, to_0_y_xxz_y, to_0_y_xxz_z, to_xxz_0, to_xxz_xy, to_xxz_yy, to_xxz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxz_x[k] = 2.0 * to_xxz_xy[k] * tke_0;

            to_0_y_xxz_y[k] = -to_xxz_0[k] + 2.0 * to_xxz_yy[k] * tke_0;

            to_0_y_xxz_z[k] = 2.0 * to_xxz_yz[k] * tke_0;
        }

        // Set up 39-42 components of targeted buffer : FP

        auto to_0_y_xyy_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 9);

        auto to_0_y_xyy_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 10);

        auto to_0_y_xyy_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_0_y_xyy_x, to_0_y_xyy_y, to_0_y_xyy_z, to_xyy_0, to_xyy_xy, to_xyy_yy, to_xyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyy_x[k] = 2.0 * to_xyy_xy[k] * tke_0;

            to_0_y_xyy_y[k] = -to_xyy_0[k] + 2.0 * to_xyy_yy[k] * tke_0;

            to_0_y_xyy_z[k] = 2.0 * to_xyy_yz[k] * tke_0;
        }

        // Set up 42-45 components of targeted buffer : FP

        auto to_0_y_xyz_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 12);

        auto to_0_y_xyz_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 13);

        auto to_0_y_xyz_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_0_y_xyz_x, to_0_y_xyz_y, to_0_y_xyz_z, to_xyz_0, to_xyz_xy, to_xyz_yy, to_xyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyz_x[k] = 2.0 * to_xyz_xy[k] * tke_0;

            to_0_y_xyz_y[k] = -to_xyz_0[k] + 2.0 * to_xyz_yy[k] * tke_0;

            to_0_y_xyz_z[k] = 2.0 * to_xyz_yz[k] * tke_0;
        }

        // Set up 45-48 components of targeted buffer : FP

        auto to_0_y_xzz_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 15);

        auto to_0_y_xzz_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 16);

        auto to_0_y_xzz_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_0_y_xzz_x, to_0_y_xzz_y, to_0_y_xzz_z, to_xzz_0, to_xzz_xy, to_xzz_yy, to_xzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzz_x[k] = 2.0 * to_xzz_xy[k] * tke_0;

            to_0_y_xzz_y[k] = -to_xzz_0[k] + 2.0 * to_xzz_yy[k] * tke_0;

            to_0_y_xzz_z[k] = 2.0 * to_xzz_yz[k] * tke_0;
        }

        // Set up 48-51 components of targeted buffer : FP

        auto to_0_y_yyy_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 18);

        auto to_0_y_yyy_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 19);

        auto to_0_y_yyy_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_0_y_yyy_x, to_0_y_yyy_y, to_0_y_yyy_z, to_yyy_0, to_yyy_xy, to_yyy_yy, to_yyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyy_x[k] = 2.0 * to_yyy_xy[k] * tke_0;

            to_0_y_yyy_y[k] = -to_yyy_0[k] + 2.0 * to_yyy_yy[k] * tke_0;

            to_0_y_yyy_z[k] = 2.0 * to_yyy_yz[k] * tke_0;
        }

        // Set up 51-54 components of targeted buffer : FP

        auto to_0_y_yyz_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 21);

        auto to_0_y_yyz_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 22);

        auto to_0_y_yyz_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_0_y_yyz_x, to_0_y_yyz_y, to_0_y_yyz_z, to_yyz_0, to_yyz_xy, to_yyz_yy, to_yyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyz_x[k] = 2.0 * to_yyz_xy[k] * tke_0;

            to_0_y_yyz_y[k] = -to_yyz_0[k] + 2.0 * to_yyz_yy[k] * tke_0;

            to_0_y_yyz_z[k] = 2.0 * to_yyz_yz[k] * tke_0;
        }

        // Set up 54-57 components of targeted buffer : FP

        auto to_0_y_yzz_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 24);

        auto to_0_y_yzz_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 25);

        auto to_0_y_yzz_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_0_y_yzz_x, to_0_y_yzz_y, to_0_y_yzz_z, to_yzz_0, to_yzz_xy, to_yzz_yy, to_yzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzz_x[k] = 2.0 * to_yzz_xy[k] * tke_0;

            to_0_y_yzz_y[k] = -to_yzz_0[k] + 2.0 * to_yzz_yy[k] * tke_0;

            to_0_y_yzz_z[k] = 2.0 * to_yzz_yz[k] * tke_0;
        }

        // Set up 57-60 components of targeted buffer : FP

        auto to_0_y_zzz_x = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 27);

        auto to_0_y_zzz_y = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 28);

        auto to_0_y_zzz_z = pbuffer.data(idx_op_geom_001_fp + 1 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_y_zzz_x, to_0_y_zzz_y, to_0_y_zzz_z, to_zzz_0, to_zzz_xy, to_zzz_yy, to_zzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzz_x[k] = 2.0 * to_zzz_xy[k] * tke_0;

            to_0_y_zzz_y[k] = -to_zzz_0[k] + 2.0 * to_zzz_yy[k] * tke_0;

            to_0_y_zzz_z[k] = 2.0 * to_zzz_yz[k] * tke_0;
        }

        // Set up 60-63 components of targeted buffer : FP

        auto to_0_z_xxx_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 0);

        auto to_0_z_xxx_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 1);

        auto to_0_z_xxx_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_0_z_xxx_x, to_0_z_xxx_y, to_0_z_xxx_z, to_xxx_0, to_xxx_xz, to_xxx_yz, to_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxx_x[k] = 2.0 * to_xxx_xz[k] * tke_0;

            to_0_z_xxx_y[k] = 2.0 * to_xxx_yz[k] * tke_0;

            to_0_z_xxx_z[k] = -to_xxx_0[k] + 2.0 * to_xxx_zz[k] * tke_0;
        }

        // Set up 63-66 components of targeted buffer : FP

        auto to_0_z_xxy_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 3);

        auto to_0_z_xxy_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 4);

        auto to_0_z_xxy_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_0_z_xxy_x, to_0_z_xxy_y, to_0_z_xxy_z, to_xxy_0, to_xxy_xz, to_xxy_yz, to_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxy_x[k] = 2.0 * to_xxy_xz[k] * tke_0;

            to_0_z_xxy_y[k] = 2.0 * to_xxy_yz[k] * tke_0;

            to_0_z_xxy_z[k] = -to_xxy_0[k] + 2.0 * to_xxy_zz[k] * tke_0;
        }

        // Set up 66-69 components of targeted buffer : FP

        auto to_0_z_xxz_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 6);

        auto to_0_z_xxz_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 7);

        auto to_0_z_xxz_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_0_z_xxz_x, to_0_z_xxz_y, to_0_z_xxz_z, to_xxz_0, to_xxz_xz, to_xxz_yz, to_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxz_x[k] = 2.0 * to_xxz_xz[k] * tke_0;

            to_0_z_xxz_y[k] = 2.0 * to_xxz_yz[k] * tke_0;

            to_0_z_xxz_z[k] = -to_xxz_0[k] + 2.0 * to_xxz_zz[k] * tke_0;
        }

        // Set up 69-72 components of targeted buffer : FP

        auto to_0_z_xyy_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 9);

        auto to_0_z_xyy_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 10);

        auto to_0_z_xyy_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_0_z_xyy_x, to_0_z_xyy_y, to_0_z_xyy_z, to_xyy_0, to_xyy_xz, to_xyy_yz, to_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyy_x[k] = 2.0 * to_xyy_xz[k] * tke_0;

            to_0_z_xyy_y[k] = 2.0 * to_xyy_yz[k] * tke_0;

            to_0_z_xyy_z[k] = -to_xyy_0[k] + 2.0 * to_xyy_zz[k] * tke_0;
        }

        // Set up 72-75 components of targeted buffer : FP

        auto to_0_z_xyz_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 12);

        auto to_0_z_xyz_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 13);

        auto to_0_z_xyz_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_0_z_xyz_x, to_0_z_xyz_y, to_0_z_xyz_z, to_xyz_0, to_xyz_xz, to_xyz_yz, to_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyz_x[k] = 2.0 * to_xyz_xz[k] * tke_0;

            to_0_z_xyz_y[k] = 2.0 * to_xyz_yz[k] * tke_0;

            to_0_z_xyz_z[k] = -to_xyz_0[k] + 2.0 * to_xyz_zz[k] * tke_0;
        }

        // Set up 75-78 components of targeted buffer : FP

        auto to_0_z_xzz_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 15);

        auto to_0_z_xzz_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 16);

        auto to_0_z_xzz_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_0_z_xzz_x, to_0_z_xzz_y, to_0_z_xzz_z, to_xzz_0, to_xzz_xz, to_xzz_yz, to_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzz_x[k] = 2.0 * to_xzz_xz[k] * tke_0;

            to_0_z_xzz_y[k] = 2.0 * to_xzz_yz[k] * tke_0;

            to_0_z_xzz_z[k] = -to_xzz_0[k] + 2.0 * to_xzz_zz[k] * tke_0;
        }

        // Set up 78-81 components of targeted buffer : FP

        auto to_0_z_yyy_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 18);

        auto to_0_z_yyy_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 19);

        auto to_0_z_yyy_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_0_z_yyy_x, to_0_z_yyy_y, to_0_z_yyy_z, to_yyy_0, to_yyy_xz, to_yyy_yz, to_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyy_x[k] = 2.0 * to_yyy_xz[k] * tke_0;

            to_0_z_yyy_y[k] = 2.0 * to_yyy_yz[k] * tke_0;

            to_0_z_yyy_z[k] = -to_yyy_0[k] + 2.0 * to_yyy_zz[k] * tke_0;
        }

        // Set up 81-84 components of targeted buffer : FP

        auto to_0_z_yyz_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 21);

        auto to_0_z_yyz_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 22);

        auto to_0_z_yyz_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_0_z_yyz_x, to_0_z_yyz_y, to_0_z_yyz_z, to_yyz_0, to_yyz_xz, to_yyz_yz, to_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyz_x[k] = 2.0 * to_yyz_xz[k] * tke_0;

            to_0_z_yyz_y[k] = 2.0 * to_yyz_yz[k] * tke_0;

            to_0_z_yyz_z[k] = -to_yyz_0[k] + 2.0 * to_yyz_zz[k] * tke_0;
        }

        // Set up 84-87 components of targeted buffer : FP

        auto to_0_z_yzz_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 24);

        auto to_0_z_yzz_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 25);

        auto to_0_z_yzz_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_0_z_yzz_x, to_0_z_yzz_y, to_0_z_yzz_z, to_yzz_0, to_yzz_xz, to_yzz_yz, to_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzz_x[k] = 2.0 * to_yzz_xz[k] * tke_0;

            to_0_z_yzz_y[k] = 2.0 * to_yzz_yz[k] * tke_0;

            to_0_z_yzz_z[k] = -to_yzz_0[k] + 2.0 * to_yzz_zz[k] * tke_0;
        }

        // Set up 87-90 components of targeted buffer : FP

        auto to_0_z_zzz_x = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 27);

        auto to_0_z_zzz_y = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 28);

        auto to_0_z_zzz_z = pbuffer.data(idx_op_geom_001_fp + 2 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_z_zzz_x, to_0_z_zzz_y, to_0_z_zzz_z, to_zzz_0, to_zzz_xz, to_zzz_yz, to_zzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzz_x[k] = 2.0 * to_zzz_xz[k] * tke_0;

            to_0_z_zzz_y[k] = 2.0 * to_zzz_yz[k] * tke_0;

            to_0_z_zzz_z[k] = -to_zzz_0[k] + 2.0 * to_zzz_zz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

