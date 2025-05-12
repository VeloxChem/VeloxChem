#include "GeometricalDerivatives0X1ForFD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_fd(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_fd,
                       const int idx_op_fp,
                       const int idx_op_ff,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : FP

        auto to_xxx_x = pbuffer.data(idx_op_fp + i * 30 + 0);

        auto to_xxx_y = pbuffer.data(idx_op_fp + i * 30 + 1);

        auto to_xxx_z = pbuffer.data(idx_op_fp + i * 30 + 2);

        auto to_xxy_x = pbuffer.data(idx_op_fp + i * 30 + 3);

        auto to_xxy_y = pbuffer.data(idx_op_fp + i * 30 + 4);

        auto to_xxy_z = pbuffer.data(idx_op_fp + i * 30 + 5);

        auto to_xxz_x = pbuffer.data(idx_op_fp + i * 30 + 6);

        auto to_xxz_y = pbuffer.data(idx_op_fp + i * 30 + 7);

        auto to_xxz_z = pbuffer.data(idx_op_fp + i * 30 + 8);

        auto to_xyy_x = pbuffer.data(idx_op_fp + i * 30 + 9);

        auto to_xyy_y = pbuffer.data(idx_op_fp + i * 30 + 10);

        auto to_xyy_z = pbuffer.data(idx_op_fp + i * 30 + 11);

        auto to_xyz_x = pbuffer.data(idx_op_fp + i * 30 + 12);

        auto to_xyz_y = pbuffer.data(idx_op_fp + i * 30 + 13);

        auto to_xyz_z = pbuffer.data(idx_op_fp + i * 30 + 14);

        auto to_xzz_x = pbuffer.data(idx_op_fp + i * 30 + 15);

        auto to_xzz_y = pbuffer.data(idx_op_fp + i * 30 + 16);

        auto to_xzz_z = pbuffer.data(idx_op_fp + i * 30 + 17);

        auto to_yyy_x = pbuffer.data(idx_op_fp + i * 30 + 18);

        auto to_yyy_y = pbuffer.data(idx_op_fp + i * 30 + 19);

        auto to_yyy_z = pbuffer.data(idx_op_fp + i * 30 + 20);

        auto to_yyz_x = pbuffer.data(idx_op_fp + i * 30 + 21);

        auto to_yyz_y = pbuffer.data(idx_op_fp + i * 30 + 22);

        auto to_yyz_z = pbuffer.data(idx_op_fp + i * 30 + 23);

        auto to_yzz_x = pbuffer.data(idx_op_fp + i * 30 + 24);

        auto to_yzz_y = pbuffer.data(idx_op_fp + i * 30 + 25);

        auto to_yzz_z = pbuffer.data(idx_op_fp + i * 30 + 26);

        auto to_zzz_x = pbuffer.data(idx_op_fp + i * 30 + 27);

        auto to_zzz_y = pbuffer.data(idx_op_fp + i * 30 + 28);

        auto to_zzz_z = pbuffer.data(idx_op_fp + i * 30 + 29);

        // Set up components of auxiliary buffer : FF

        auto to_xxx_xxx = pbuffer.data(idx_op_ff + i * 100 + 0);

        auto to_xxx_xxy = pbuffer.data(idx_op_ff + i * 100 + 1);

        auto to_xxx_xxz = pbuffer.data(idx_op_ff + i * 100 + 2);

        auto to_xxx_xyy = pbuffer.data(idx_op_ff + i * 100 + 3);

        auto to_xxx_xyz = pbuffer.data(idx_op_ff + i * 100 + 4);

        auto to_xxx_xzz = pbuffer.data(idx_op_ff + i * 100 + 5);

        auto to_xxx_yyy = pbuffer.data(idx_op_ff + i * 100 + 6);

        auto to_xxx_yyz = pbuffer.data(idx_op_ff + i * 100 + 7);

        auto to_xxx_yzz = pbuffer.data(idx_op_ff + i * 100 + 8);

        auto to_xxx_zzz = pbuffer.data(idx_op_ff + i * 100 + 9);

        auto to_xxy_xxx = pbuffer.data(idx_op_ff + i * 100 + 10);

        auto to_xxy_xxy = pbuffer.data(idx_op_ff + i * 100 + 11);

        auto to_xxy_xxz = pbuffer.data(idx_op_ff + i * 100 + 12);

        auto to_xxy_xyy = pbuffer.data(idx_op_ff + i * 100 + 13);

        auto to_xxy_xyz = pbuffer.data(idx_op_ff + i * 100 + 14);

        auto to_xxy_xzz = pbuffer.data(idx_op_ff + i * 100 + 15);

        auto to_xxy_yyy = pbuffer.data(idx_op_ff + i * 100 + 16);

        auto to_xxy_yyz = pbuffer.data(idx_op_ff + i * 100 + 17);

        auto to_xxy_yzz = pbuffer.data(idx_op_ff + i * 100 + 18);

        auto to_xxy_zzz = pbuffer.data(idx_op_ff + i * 100 + 19);

        auto to_xxz_xxx = pbuffer.data(idx_op_ff + i * 100 + 20);

        auto to_xxz_xxy = pbuffer.data(idx_op_ff + i * 100 + 21);

        auto to_xxz_xxz = pbuffer.data(idx_op_ff + i * 100 + 22);

        auto to_xxz_xyy = pbuffer.data(idx_op_ff + i * 100 + 23);

        auto to_xxz_xyz = pbuffer.data(idx_op_ff + i * 100 + 24);

        auto to_xxz_xzz = pbuffer.data(idx_op_ff + i * 100 + 25);

        auto to_xxz_yyy = pbuffer.data(idx_op_ff + i * 100 + 26);

        auto to_xxz_yyz = pbuffer.data(idx_op_ff + i * 100 + 27);

        auto to_xxz_yzz = pbuffer.data(idx_op_ff + i * 100 + 28);

        auto to_xxz_zzz = pbuffer.data(idx_op_ff + i * 100 + 29);

        auto to_xyy_xxx = pbuffer.data(idx_op_ff + i * 100 + 30);

        auto to_xyy_xxy = pbuffer.data(idx_op_ff + i * 100 + 31);

        auto to_xyy_xxz = pbuffer.data(idx_op_ff + i * 100 + 32);

        auto to_xyy_xyy = pbuffer.data(idx_op_ff + i * 100 + 33);

        auto to_xyy_xyz = pbuffer.data(idx_op_ff + i * 100 + 34);

        auto to_xyy_xzz = pbuffer.data(idx_op_ff + i * 100 + 35);

        auto to_xyy_yyy = pbuffer.data(idx_op_ff + i * 100 + 36);

        auto to_xyy_yyz = pbuffer.data(idx_op_ff + i * 100 + 37);

        auto to_xyy_yzz = pbuffer.data(idx_op_ff + i * 100 + 38);

        auto to_xyy_zzz = pbuffer.data(idx_op_ff + i * 100 + 39);

        auto to_xyz_xxx = pbuffer.data(idx_op_ff + i * 100 + 40);

        auto to_xyz_xxy = pbuffer.data(idx_op_ff + i * 100 + 41);

        auto to_xyz_xxz = pbuffer.data(idx_op_ff + i * 100 + 42);

        auto to_xyz_xyy = pbuffer.data(idx_op_ff + i * 100 + 43);

        auto to_xyz_xyz = pbuffer.data(idx_op_ff + i * 100 + 44);

        auto to_xyz_xzz = pbuffer.data(idx_op_ff + i * 100 + 45);

        auto to_xyz_yyy = pbuffer.data(idx_op_ff + i * 100 + 46);

        auto to_xyz_yyz = pbuffer.data(idx_op_ff + i * 100 + 47);

        auto to_xyz_yzz = pbuffer.data(idx_op_ff + i * 100 + 48);

        auto to_xyz_zzz = pbuffer.data(idx_op_ff + i * 100 + 49);

        auto to_xzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 50);

        auto to_xzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 51);

        auto to_xzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 52);

        auto to_xzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 53);

        auto to_xzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 54);

        auto to_xzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 55);

        auto to_xzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 56);

        auto to_xzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 57);

        auto to_xzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 58);

        auto to_xzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 59);

        auto to_yyy_xxx = pbuffer.data(idx_op_ff + i * 100 + 60);

        auto to_yyy_xxy = pbuffer.data(idx_op_ff + i * 100 + 61);

        auto to_yyy_xxz = pbuffer.data(idx_op_ff + i * 100 + 62);

        auto to_yyy_xyy = pbuffer.data(idx_op_ff + i * 100 + 63);

        auto to_yyy_xyz = pbuffer.data(idx_op_ff + i * 100 + 64);

        auto to_yyy_xzz = pbuffer.data(idx_op_ff + i * 100 + 65);

        auto to_yyy_yyy = pbuffer.data(idx_op_ff + i * 100 + 66);

        auto to_yyy_yyz = pbuffer.data(idx_op_ff + i * 100 + 67);

        auto to_yyy_yzz = pbuffer.data(idx_op_ff + i * 100 + 68);

        auto to_yyy_zzz = pbuffer.data(idx_op_ff + i * 100 + 69);

        auto to_yyz_xxx = pbuffer.data(idx_op_ff + i * 100 + 70);

        auto to_yyz_xxy = pbuffer.data(idx_op_ff + i * 100 + 71);

        auto to_yyz_xxz = pbuffer.data(idx_op_ff + i * 100 + 72);

        auto to_yyz_xyy = pbuffer.data(idx_op_ff + i * 100 + 73);

        auto to_yyz_xyz = pbuffer.data(idx_op_ff + i * 100 + 74);

        auto to_yyz_xzz = pbuffer.data(idx_op_ff + i * 100 + 75);

        auto to_yyz_yyy = pbuffer.data(idx_op_ff + i * 100 + 76);

        auto to_yyz_yyz = pbuffer.data(idx_op_ff + i * 100 + 77);

        auto to_yyz_yzz = pbuffer.data(idx_op_ff + i * 100 + 78);

        auto to_yyz_zzz = pbuffer.data(idx_op_ff + i * 100 + 79);

        auto to_yzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 80);

        auto to_yzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 81);

        auto to_yzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 82);

        auto to_yzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 83);

        auto to_yzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 84);

        auto to_yzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 85);

        auto to_yzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 86);

        auto to_yzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 87);

        auto to_yzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 88);

        auto to_yzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 89);

        auto to_zzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 90);

        auto to_zzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 91);

        auto to_zzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 92);

        auto to_zzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 93);

        auto to_zzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 94);

        auto to_zzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 95);

        auto to_zzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 96);

        auto to_zzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 97);

        auto to_zzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 98);

        auto to_zzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 99);

        // Set up 0-6 components of targeted buffer : FD

        auto to_0_x_xxx_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 0);

        auto to_0_x_xxx_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 1);

        auto to_0_x_xxx_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 2);

        auto to_0_x_xxx_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 3);

        auto to_0_x_xxx_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 4);

        auto to_0_x_xxx_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_0_x_xxx_xx, to_0_x_xxx_xy, to_0_x_xxx_xz, to_0_x_xxx_yy, to_0_x_xxx_yz, to_0_x_xxx_zz, to_xxx_x, to_xxx_xxx, to_xxx_xxy, to_xxx_xxz, to_xxx_xyy, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxx_xx[k] = -2.0 * to_xxx_x[k] + 2.0 * to_xxx_xxx[k] * tke_0;

            to_0_x_xxx_xy[k] = -to_xxx_y[k] + 2.0 * to_xxx_xxy[k] * tke_0;

            to_0_x_xxx_xz[k] = -to_xxx_z[k] + 2.0 * to_xxx_xxz[k] * tke_0;

            to_0_x_xxx_yy[k] = 2.0 * to_xxx_xyy[k] * tke_0;

            to_0_x_xxx_yz[k] = 2.0 * to_xxx_xyz[k] * tke_0;

            to_0_x_xxx_zz[k] = 2.0 * to_xxx_xzz[k] * tke_0;
        }

        // Set up 6-12 components of targeted buffer : FD

        auto to_0_x_xxy_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 6);

        auto to_0_x_xxy_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 7);

        auto to_0_x_xxy_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 8);

        auto to_0_x_xxy_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 9);

        auto to_0_x_xxy_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 10);

        auto to_0_x_xxy_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_0_x_xxy_xx, to_0_x_xxy_xy, to_0_x_xxy_xz, to_0_x_xxy_yy, to_0_x_xxy_yz, to_0_x_xxy_zz, to_xxy_x, to_xxy_xxx, to_xxy_xxy, to_xxy_xxz, to_xxy_xyy, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxy_xx[k] = -2.0 * to_xxy_x[k] + 2.0 * to_xxy_xxx[k] * tke_0;

            to_0_x_xxy_xy[k] = -to_xxy_y[k] + 2.0 * to_xxy_xxy[k] * tke_0;

            to_0_x_xxy_xz[k] = -to_xxy_z[k] + 2.0 * to_xxy_xxz[k] * tke_0;

            to_0_x_xxy_yy[k] = 2.0 * to_xxy_xyy[k] * tke_0;

            to_0_x_xxy_yz[k] = 2.0 * to_xxy_xyz[k] * tke_0;

            to_0_x_xxy_zz[k] = 2.0 * to_xxy_xzz[k] * tke_0;
        }

        // Set up 12-18 components of targeted buffer : FD

        auto to_0_x_xxz_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 12);

        auto to_0_x_xxz_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 13);

        auto to_0_x_xxz_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 14);

        auto to_0_x_xxz_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 15);

        auto to_0_x_xxz_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 16);

        auto to_0_x_xxz_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_0_x_xxz_xx, to_0_x_xxz_xy, to_0_x_xxz_xz, to_0_x_xxz_yy, to_0_x_xxz_yz, to_0_x_xxz_zz, to_xxz_x, to_xxz_xxx, to_xxz_xxy, to_xxz_xxz, to_xxz_xyy, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxz_xx[k] = -2.0 * to_xxz_x[k] + 2.0 * to_xxz_xxx[k] * tke_0;

            to_0_x_xxz_xy[k] = -to_xxz_y[k] + 2.0 * to_xxz_xxy[k] * tke_0;

            to_0_x_xxz_xz[k] = -to_xxz_z[k] + 2.0 * to_xxz_xxz[k] * tke_0;

            to_0_x_xxz_yy[k] = 2.0 * to_xxz_xyy[k] * tke_0;

            to_0_x_xxz_yz[k] = 2.0 * to_xxz_xyz[k] * tke_0;

            to_0_x_xxz_zz[k] = 2.0 * to_xxz_xzz[k] * tke_0;
        }

        // Set up 18-24 components of targeted buffer : FD

        auto to_0_x_xyy_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 18);

        auto to_0_x_xyy_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 19);

        auto to_0_x_xyy_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 20);

        auto to_0_x_xyy_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 21);

        auto to_0_x_xyy_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 22);

        auto to_0_x_xyy_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_0_x_xyy_xx, to_0_x_xyy_xy, to_0_x_xyy_xz, to_0_x_xyy_yy, to_0_x_xyy_yz, to_0_x_xyy_zz, to_xyy_x, to_xyy_xxx, to_xyy_xxy, to_xyy_xxz, to_xyy_xyy, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyy_xx[k] = -2.0 * to_xyy_x[k] + 2.0 * to_xyy_xxx[k] * tke_0;

            to_0_x_xyy_xy[k] = -to_xyy_y[k] + 2.0 * to_xyy_xxy[k] * tke_0;

            to_0_x_xyy_xz[k] = -to_xyy_z[k] + 2.0 * to_xyy_xxz[k] * tke_0;

            to_0_x_xyy_yy[k] = 2.0 * to_xyy_xyy[k] * tke_0;

            to_0_x_xyy_yz[k] = 2.0 * to_xyy_xyz[k] * tke_0;

            to_0_x_xyy_zz[k] = 2.0 * to_xyy_xzz[k] * tke_0;
        }

        // Set up 24-30 components of targeted buffer : FD

        auto to_0_x_xyz_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 24);

        auto to_0_x_xyz_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 25);

        auto to_0_x_xyz_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 26);

        auto to_0_x_xyz_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 27);

        auto to_0_x_xyz_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 28);

        auto to_0_x_xyz_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_0_x_xyz_xx, to_0_x_xyz_xy, to_0_x_xyz_xz, to_0_x_xyz_yy, to_0_x_xyz_yz, to_0_x_xyz_zz, to_xyz_x, to_xyz_xxx, to_xyz_xxy, to_xyz_xxz, to_xyz_xyy, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyz_xx[k] = -2.0 * to_xyz_x[k] + 2.0 * to_xyz_xxx[k] * tke_0;

            to_0_x_xyz_xy[k] = -to_xyz_y[k] + 2.0 * to_xyz_xxy[k] * tke_0;

            to_0_x_xyz_xz[k] = -to_xyz_z[k] + 2.0 * to_xyz_xxz[k] * tke_0;

            to_0_x_xyz_yy[k] = 2.0 * to_xyz_xyy[k] * tke_0;

            to_0_x_xyz_yz[k] = 2.0 * to_xyz_xyz[k] * tke_0;

            to_0_x_xyz_zz[k] = 2.0 * to_xyz_xzz[k] * tke_0;
        }

        // Set up 30-36 components of targeted buffer : FD

        auto to_0_x_xzz_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 30);

        auto to_0_x_xzz_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 31);

        auto to_0_x_xzz_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 32);

        auto to_0_x_xzz_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 33);

        auto to_0_x_xzz_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 34);

        auto to_0_x_xzz_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_0_x_xzz_xx, to_0_x_xzz_xy, to_0_x_xzz_xz, to_0_x_xzz_yy, to_0_x_xzz_yz, to_0_x_xzz_zz, to_xzz_x, to_xzz_xxx, to_xzz_xxy, to_xzz_xxz, to_xzz_xyy, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzz_xx[k] = -2.0 * to_xzz_x[k] + 2.0 * to_xzz_xxx[k] * tke_0;

            to_0_x_xzz_xy[k] = -to_xzz_y[k] + 2.0 * to_xzz_xxy[k] * tke_0;

            to_0_x_xzz_xz[k] = -to_xzz_z[k] + 2.0 * to_xzz_xxz[k] * tke_0;

            to_0_x_xzz_yy[k] = 2.0 * to_xzz_xyy[k] * tke_0;

            to_0_x_xzz_yz[k] = 2.0 * to_xzz_xyz[k] * tke_0;

            to_0_x_xzz_zz[k] = 2.0 * to_xzz_xzz[k] * tke_0;
        }

        // Set up 36-42 components of targeted buffer : FD

        auto to_0_x_yyy_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 36);

        auto to_0_x_yyy_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 37);

        auto to_0_x_yyy_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 38);

        auto to_0_x_yyy_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 39);

        auto to_0_x_yyy_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 40);

        auto to_0_x_yyy_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_0_x_yyy_xx, to_0_x_yyy_xy, to_0_x_yyy_xz, to_0_x_yyy_yy, to_0_x_yyy_yz, to_0_x_yyy_zz, to_yyy_x, to_yyy_xxx, to_yyy_xxy, to_yyy_xxz, to_yyy_xyy, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyy_xx[k] = -2.0 * to_yyy_x[k] + 2.0 * to_yyy_xxx[k] * tke_0;

            to_0_x_yyy_xy[k] = -to_yyy_y[k] + 2.0 * to_yyy_xxy[k] * tke_0;

            to_0_x_yyy_xz[k] = -to_yyy_z[k] + 2.0 * to_yyy_xxz[k] * tke_0;

            to_0_x_yyy_yy[k] = 2.0 * to_yyy_xyy[k] * tke_0;

            to_0_x_yyy_yz[k] = 2.0 * to_yyy_xyz[k] * tke_0;

            to_0_x_yyy_zz[k] = 2.0 * to_yyy_xzz[k] * tke_0;
        }

        // Set up 42-48 components of targeted buffer : FD

        auto to_0_x_yyz_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 42);

        auto to_0_x_yyz_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 43);

        auto to_0_x_yyz_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 44);

        auto to_0_x_yyz_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 45);

        auto to_0_x_yyz_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 46);

        auto to_0_x_yyz_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_0_x_yyz_xx, to_0_x_yyz_xy, to_0_x_yyz_xz, to_0_x_yyz_yy, to_0_x_yyz_yz, to_0_x_yyz_zz, to_yyz_x, to_yyz_xxx, to_yyz_xxy, to_yyz_xxz, to_yyz_xyy, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyz_xx[k] = -2.0 * to_yyz_x[k] + 2.0 * to_yyz_xxx[k] * tke_0;

            to_0_x_yyz_xy[k] = -to_yyz_y[k] + 2.0 * to_yyz_xxy[k] * tke_0;

            to_0_x_yyz_xz[k] = -to_yyz_z[k] + 2.0 * to_yyz_xxz[k] * tke_0;

            to_0_x_yyz_yy[k] = 2.0 * to_yyz_xyy[k] * tke_0;

            to_0_x_yyz_yz[k] = 2.0 * to_yyz_xyz[k] * tke_0;

            to_0_x_yyz_zz[k] = 2.0 * to_yyz_xzz[k] * tke_0;
        }

        // Set up 48-54 components of targeted buffer : FD

        auto to_0_x_yzz_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 48);

        auto to_0_x_yzz_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 49);

        auto to_0_x_yzz_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 50);

        auto to_0_x_yzz_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 51);

        auto to_0_x_yzz_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 52);

        auto to_0_x_yzz_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_0_x_yzz_xx, to_0_x_yzz_xy, to_0_x_yzz_xz, to_0_x_yzz_yy, to_0_x_yzz_yz, to_0_x_yzz_zz, to_yzz_x, to_yzz_xxx, to_yzz_xxy, to_yzz_xxz, to_yzz_xyy, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzz_xx[k] = -2.0 * to_yzz_x[k] + 2.0 * to_yzz_xxx[k] * tke_0;

            to_0_x_yzz_xy[k] = -to_yzz_y[k] + 2.0 * to_yzz_xxy[k] * tke_0;

            to_0_x_yzz_xz[k] = -to_yzz_z[k] + 2.0 * to_yzz_xxz[k] * tke_0;

            to_0_x_yzz_yy[k] = 2.0 * to_yzz_xyy[k] * tke_0;

            to_0_x_yzz_yz[k] = 2.0 * to_yzz_xyz[k] * tke_0;

            to_0_x_yzz_zz[k] = 2.0 * to_yzz_xzz[k] * tke_0;
        }

        // Set up 54-60 components of targeted buffer : FD

        auto to_0_x_zzz_xx = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 54);

        auto to_0_x_zzz_xy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 55);

        auto to_0_x_zzz_xz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 56);

        auto to_0_x_zzz_yy = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 57);

        auto to_0_x_zzz_yz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 58);

        auto to_0_x_zzz_zz = pbuffer.data(idx_op_geom_001_fd + 0 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_0_x_zzz_xx, to_0_x_zzz_xy, to_0_x_zzz_xz, to_0_x_zzz_yy, to_0_x_zzz_yz, to_0_x_zzz_zz, to_zzz_x, to_zzz_xxx, to_zzz_xxy, to_zzz_xxz, to_zzz_xyy, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzz_xx[k] = -2.0 * to_zzz_x[k] + 2.0 * to_zzz_xxx[k] * tke_0;

            to_0_x_zzz_xy[k] = -to_zzz_y[k] + 2.0 * to_zzz_xxy[k] * tke_0;

            to_0_x_zzz_xz[k] = -to_zzz_z[k] + 2.0 * to_zzz_xxz[k] * tke_0;

            to_0_x_zzz_yy[k] = 2.0 * to_zzz_xyy[k] * tke_0;

            to_0_x_zzz_yz[k] = 2.0 * to_zzz_xyz[k] * tke_0;

            to_0_x_zzz_zz[k] = 2.0 * to_zzz_xzz[k] * tke_0;
        }

        // Set up 60-66 components of targeted buffer : FD

        auto to_0_y_xxx_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 0);

        auto to_0_y_xxx_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 1);

        auto to_0_y_xxx_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 2);

        auto to_0_y_xxx_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 3);

        auto to_0_y_xxx_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 4);

        auto to_0_y_xxx_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_0_y_xxx_xx, to_0_y_xxx_xy, to_0_y_xxx_xz, to_0_y_xxx_yy, to_0_y_xxx_yz, to_0_y_xxx_zz, to_xxx_x, to_xxx_xxy, to_xxx_xyy, to_xxx_xyz, to_xxx_y, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxx_xx[k] = 2.0 * to_xxx_xxy[k] * tke_0;

            to_0_y_xxx_xy[k] = -to_xxx_x[k] + 2.0 * to_xxx_xyy[k] * tke_0;

            to_0_y_xxx_xz[k] = 2.0 * to_xxx_xyz[k] * tke_0;

            to_0_y_xxx_yy[k] = -2.0 * to_xxx_y[k] + 2.0 * to_xxx_yyy[k] * tke_0;

            to_0_y_xxx_yz[k] = -to_xxx_z[k] + 2.0 * to_xxx_yyz[k] * tke_0;

            to_0_y_xxx_zz[k] = 2.0 * to_xxx_yzz[k] * tke_0;
        }

        // Set up 66-72 components of targeted buffer : FD

        auto to_0_y_xxy_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 6);

        auto to_0_y_xxy_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 7);

        auto to_0_y_xxy_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 8);

        auto to_0_y_xxy_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 9);

        auto to_0_y_xxy_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 10);

        auto to_0_y_xxy_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_0_y_xxy_xx, to_0_y_xxy_xy, to_0_y_xxy_xz, to_0_y_xxy_yy, to_0_y_xxy_yz, to_0_y_xxy_zz, to_xxy_x, to_xxy_xxy, to_xxy_xyy, to_xxy_xyz, to_xxy_y, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxy_xx[k] = 2.0 * to_xxy_xxy[k] * tke_0;

            to_0_y_xxy_xy[k] = -to_xxy_x[k] + 2.0 * to_xxy_xyy[k] * tke_0;

            to_0_y_xxy_xz[k] = 2.0 * to_xxy_xyz[k] * tke_0;

            to_0_y_xxy_yy[k] = -2.0 * to_xxy_y[k] + 2.0 * to_xxy_yyy[k] * tke_0;

            to_0_y_xxy_yz[k] = -to_xxy_z[k] + 2.0 * to_xxy_yyz[k] * tke_0;

            to_0_y_xxy_zz[k] = 2.0 * to_xxy_yzz[k] * tke_0;
        }

        // Set up 72-78 components of targeted buffer : FD

        auto to_0_y_xxz_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 12);

        auto to_0_y_xxz_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 13);

        auto to_0_y_xxz_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 14);

        auto to_0_y_xxz_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 15);

        auto to_0_y_xxz_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 16);

        auto to_0_y_xxz_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_0_y_xxz_xx, to_0_y_xxz_xy, to_0_y_xxz_xz, to_0_y_xxz_yy, to_0_y_xxz_yz, to_0_y_xxz_zz, to_xxz_x, to_xxz_xxy, to_xxz_xyy, to_xxz_xyz, to_xxz_y, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxz_xx[k] = 2.0 * to_xxz_xxy[k] * tke_0;

            to_0_y_xxz_xy[k] = -to_xxz_x[k] + 2.0 * to_xxz_xyy[k] * tke_0;

            to_0_y_xxz_xz[k] = 2.0 * to_xxz_xyz[k] * tke_0;

            to_0_y_xxz_yy[k] = -2.0 * to_xxz_y[k] + 2.0 * to_xxz_yyy[k] * tke_0;

            to_0_y_xxz_yz[k] = -to_xxz_z[k] + 2.0 * to_xxz_yyz[k] * tke_0;

            to_0_y_xxz_zz[k] = 2.0 * to_xxz_yzz[k] * tke_0;
        }

        // Set up 78-84 components of targeted buffer : FD

        auto to_0_y_xyy_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 18);

        auto to_0_y_xyy_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 19);

        auto to_0_y_xyy_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 20);

        auto to_0_y_xyy_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 21);

        auto to_0_y_xyy_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 22);

        auto to_0_y_xyy_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_0_y_xyy_xx, to_0_y_xyy_xy, to_0_y_xyy_xz, to_0_y_xyy_yy, to_0_y_xyy_yz, to_0_y_xyy_zz, to_xyy_x, to_xyy_xxy, to_xyy_xyy, to_xyy_xyz, to_xyy_y, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyy_xx[k] = 2.0 * to_xyy_xxy[k] * tke_0;

            to_0_y_xyy_xy[k] = -to_xyy_x[k] + 2.0 * to_xyy_xyy[k] * tke_0;

            to_0_y_xyy_xz[k] = 2.0 * to_xyy_xyz[k] * tke_0;

            to_0_y_xyy_yy[k] = -2.0 * to_xyy_y[k] + 2.0 * to_xyy_yyy[k] * tke_0;

            to_0_y_xyy_yz[k] = -to_xyy_z[k] + 2.0 * to_xyy_yyz[k] * tke_0;

            to_0_y_xyy_zz[k] = 2.0 * to_xyy_yzz[k] * tke_0;
        }

        // Set up 84-90 components of targeted buffer : FD

        auto to_0_y_xyz_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 24);

        auto to_0_y_xyz_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 25);

        auto to_0_y_xyz_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 26);

        auto to_0_y_xyz_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 27);

        auto to_0_y_xyz_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 28);

        auto to_0_y_xyz_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_0_y_xyz_xx, to_0_y_xyz_xy, to_0_y_xyz_xz, to_0_y_xyz_yy, to_0_y_xyz_yz, to_0_y_xyz_zz, to_xyz_x, to_xyz_xxy, to_xyz_xyy, to_xyz_xyz, to_xyz_y, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyz_xx[k] = 2.0 * to_xyz_xxy[k] * tke_0;

            to_0_y_xyz_xy[k] = -to_xyz_x[k] + 2.0 * to_xyz_xyy[k] * tke_0;

            to_0_y_xyz_xz[k] = 2.0 * to_xyz_xyz[k] * tke_0;

            to_0_y_xyz_yy[k] = -2.0 * to_xyz_y[k] + 2.0 * to_xyz_yyy[k] * tke_0;

            to_0_y_xyz_yz[k] = -to_xyz_z[k] + 2.0 * to_xyz_yyz[k] * tke_0;

            to_0_y_xyz_zz[k] = 2.0 * to_xyz_yzz[k] * tke_0;
        }

        // Set up 90-96 components of targeted buffer : FD

        auto to_0_y_xzz_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 30);

        auto to_0_y_xzz_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 31);

        auto to_0_y_xzz_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 32);

        auto to_0_y_xzz_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 33);

        auto to_0_y_xzz_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 34);

        auto to_0_y_xzz_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_0_y_xzz_xx, to_0_y_xzz_xy, to_0_y_xzz_xz, to_0_y_xzz_yy, to_0_y_xzz_yz, to_0_y_xzz_zz, to_xzz_x, to_xzz_xxy, to_xzz_xyy, to_xzz_xyz, to_xzz_y, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzz_xx[k] = 2.0 * to_xzz_xxy[k] * tke_0;

            to_0_y_xzz_xy[k] = -to_xzz_x[k] + 2.0 * to_xzz_xyy[k] * tke_0;

            to_0_y_xzz_xz[k] = 2.0 * to_xzz_xyz[k] * tke_0;

            to_0_y_xzz_yy[k] = -2.0 * to_xzz_y[k] + 2.0 * to_xzz_yyy[k] * tke_0;

            to_0_y_xzz_yz[k] = -to_xzz_z[k] + 2.0 * to_xzz_yyz[k] * tke_0;

            to_0_y_xzz_zz[k] = 2.0 * to_xzz_yzz[k] * tke_0;
        }

        // Set up 96-102 components of targeted buffer : FD

        auto to_0_y_yyy_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 36);

        auto to_0_y_yyy_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 37);

        auto to_0_y_yyy_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 38);

        auto to_0_y_yyy_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 39);

        auto to_0_y_yyy_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 40);

        auto to_0_y_yyy_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_0_y_yyy_xx, to_0_y_yyy_xy, to_0_y_yyy_xz, to_0_y_yyy_yy, to_0_y_yyy_yz, to_0_y_yyy_zz, to_yyy_x, to_yyy_xxy, to_yyy_xyy, to_yyy_xyz, to_yyy_y, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyy_xx[k] = 2.0 * to_yyy_xxy[k] * tke_0;

            to_0_y_yyy_xy[k] = -to_yyy_x[k] + 2.0 * to_yyy_xyy[k] * tke_0;

            to_0_y_yyy_xz[k] = 2.0 * to_yyy_xyz[k] * tke_0;

            to_0_y_yyy_yy[k] = -2.0 * to_yyy_y[k] + 2.0 * to_yyy_yyy[k] * tke_0;

            to_0_y_yyy_yz[k] = -to_yyy_z[k] + 2.0 * to_yyy_yyz[k] * tke_0;

            to_0_y_yyy_zz[k] = 2.0 * to_yyy_yzz[k] * tke_0;
        }

        // Set up 102-108 components of targeted buffer : FD

        auto to_0_y_yyz_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 42);

        auto to_0_y_yyz_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 43);

        auto to_0_y_yyz_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 44);

        auto to_0_y_yyz_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 45);

        auto to_0_y_yyz_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 46);

        auto to_0_y_yyz_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_0_y_yyz_xx, to_0_y_yyz_xy, to_0_y_yyz_xz, to_0_y_yyz_yy, to_0_y_yyz_yz, to_0_y_yyz_zz, to_yyz_x, to_yyz_xxy, to_yyz_xyy, to_yyz_xyz, to_yyz_y, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyz_xx[k] = 2.0 * to_yyz_xxy[k] * tke_0;

            to_0_y_yyz_xy[k] = -to_yyz_x[k] + 2.0 * to_yyz_xyy[k] * tke_0;

            to_0_y_yyz_xz[k] = 2.0 * to_yyz_xyz[k] * tke_0;

            to_0_y_yyz_yy[k] = -2.0 * to_yyz_y[k] + 2.0 * to_yyz_yyy[k] * tke_0;

            to_0_y_yyz_yz[k] = -to_yyz_z[k] + 2.0 * to_yyz_yyz[k] * tke_0;

            to_0_y_yyz_zz[k] = 2.0 * to_yyz_yzz[k] * tke_0;
        }

        // Set up 108-114 components of targeted buffer : FD

        auto to_0_y_yzz_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 48);

        auto to_0_y_yzz_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 49);

        auto to_0_y_yzz_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 50);

        auto to_0_y_yzz_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 51);

        auto to_0_y_yzz_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 52);

        auto to_0_y_yzz_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_0_y_yzz_xx, to_0_y_yzz_xy, to_0_y_yzz_xz, to_0_y_yzz_yy, to_0_y_yzz_yz, to_0_y_yzz_zz, to_yzz_x, to_yzz_xxy, to_yzz_xyy, to_yzz_xyz, to_yzz_y, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzz_xx[k] = 2.0 * to_yzz_xxy[k] * tke_0;

            to_0_y_yzz_xy[k] = -to_yzz_x[k] + 2.0 * to_yzz_xyy[k] * tke_0;

            to_0_y_yzz_xz[k] = 2.0 * to_yzz_xyz[k] * tke_0;

            to_0_y_yzz_yy[k] = -2.0 * to_yzz_y[k] + 2.0 * to_yzz_yyy[k] * tke_0;

            to_0_y_yzz_yz[k] = -to_yzz_z[k] + 2.0 * to_yzz_yyz[k] * tke_0;

            to_0_y_yzz_zz[k] = 2.0 * to_yzz_yzz[k] * tke_0;
        }

        // Set up 114-120 components of targeted buffer : FD

        auto to_0_y_zzz_xx = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 54);

        auto to_0_y_zzz_xy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 55);

        auto to_0_y_zzz_xz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 56);

        auto to_0_y_zzz_yy = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 57);

        auto to_0_y_zzz_yz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 58);

        auto to_0_y_zzz_zz = pbuffer.data(idx_op_geom_001_fd + 1 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_0_y_zzz_xx, to_0_y_zzz_xy, to_0_y_zzz_xz, to_0_y_zzz_yy, to_0_y_zzz_yz, to_0_y_zzz_zz, to_zzz_x, to_zzz_xxy, to_zzz_xyy, to_zzz_xyz, to_zzz_y, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzz_xx[k] = 2.0 * to_zzz_xxy[k] * tke_0;

            to_0_y_zzz_xy[k] = -to_zzz_x[k] + 2.0 * to_zzz_xyy[k] * tke_0;

            to_0_y_zzz_xz[k] = 2.0 * to_zzz_xyz[k] * tke_0;

            to_0_y_zzz_yy[k] = -2.0 * to_zzz_y[k] + 2.0 * to_zzz_yyy[k] * tke_0;

            to_0_y_zzz_yz[k] = -to_zzz_z[k] + 2.0 * to_zzz_yyz[k] * tke_0;

            to_0_y_zzz_zz[k] = 2.0 * to_zzz_yzz[k] * tke_0;
        }

        // Set up 120-126 components of targeted buffer : FD

        auto to_0_z_xxx_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 0);

        auto to_0_z_xxx_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 1);

        auto to_0_z_xxx_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 2);

        auto to_0_z_xxx_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 3);

        auto to_0_z_xxx_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 4);

        auto to_0_z_xxx_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_0_z_xxx_xx, to_0_z_xxx_xy, to_0_z_xxx_xz, to_0_z_xxx_yy, to_0_z_xxx_yz, to_0_z_xxx_zz, to_xxx_x, to_xxx_xxz, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxx_xx[k] = 2.0 * to_xxx_xxz[k] * tke_0;

            to_0_z_xxx_xy[k] = 2.0 * to_xxx_xyz[k] * tke_0;

            to_0_z_xxx_xz[k] = -to_xxx_x[k] + 2.0 * to_xxx_xzz[k] * tke_0;

            to_0_z_xxx_yy[k] = 2.0 * to_xxx_yyz[k] * tke_0;

            to_0_z_xxx_yz[k] = -to_xxx_y[k] + 2.0 * to_xxx_yzz[k] * tke_0;

            to_0_z_xxx_zz[k] = -2.0 * to_xxx_z[k] + 2.0 * to_xxx_zzz[k] * tke_0;
        }

        // Set up 126-132 components of targeted buffer : FD

        auto to_0_z_xxy_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 6);

        auto to_0_z_xxy_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 7);

        auto to_0_z_xxy_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 8);

        auto to_0_z_xxy_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 9);

        auto to_0_z_xxy_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 10);

        auto to_0_z_xxy_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_0_z_xxy_xx, to_0_z_xxy_xy, to_0_z_xxy_xz, to_0_z_xxy_yy, to_0_z_xxy_yz, to_0_z_xxy_zz, to_xxy_x, to_xxy_xxz, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxy_xx[k] = 2.0 * to_xxy_xxz[k] * tke_0;

            to_0_z_xxy_xy[k] = 2.0 * to_xxy_xyz[k] * tke_0;

            to_0_z_xxy_xz[k] = -to_xxy_x[k] + 2.0 * to_xxy_xzz[k] * tke_0;

            to_0_z_xxy_yy[k] = 2.0 * to_xxy_yyz[k] * tke_0;

            to_0_z_xxy_yz[k] = -to_xxy_y[k] + 2.0 * to_xxy_yzz[k] * tke_0;

            to_0_z_xxy_zz[k] = -2.0 * to_xxy_z[k] + 2.0 * to_xxy_zzz[k] * tke_0;
        }

        // Set up 132-138 components of targeted buffer : FD

        auto to_0_z_xxz_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 12);

        auto to_0_z_xxz_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 13);

        auto to_0_z_xxz_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 14);

        auto to_0_z_xxz_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 15);

        auto to_0_z_xxz_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 16);

        auto to_0_z_xxz_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_0_z_xxz_xx, to_0_z_xxz_xy, to_0_z_xxz_xz, to_0_z_xxz_yy, to_0_z_xxz_yz, to_0_z_xxz_zz, to_xxz_x, to_xxz_xxz, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxz_xx[k] = 2.0 * to_xxz_xxz[k] * tke_0;

            to_0_z_xxz_xy[k] = 2.0 * to_xxz_xyz[k] * tke_0;

            to_0_z_xxz_xz[k] = -to_xxz_x[k] + 2.0 * to_xxz_xzz[k] * tke_0;

            to_0_z_xxz_yy[k] = 2.0 * to_xxz_yyz[k] * tke_0;

            to_0_z_xxz_yz[k] = -to_xxz_y[k] + 2.0 * to_xxz_yzz[k] * tke_0;

            to_0_z_xxz_zz[k] = -2.0 * to_xxz_z[k] + 2.0 * to_xxz_zzz[k] * tke_0;
        }

        // Set up 138-144 components of targeted buffer : FD

        auto to_0_z_xyy_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 18);

        auto to_0_z_xyy_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 19);

        auto to_0_z_xyy_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 20);

        auto to_0_z_xyy_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 21);

        auto to_0_z_xyy_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 22);

        auto to_0_z_xyy_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_0_z_xyy_xx, to_0_z_xyy_xy, to_0_z_xyy_xz, to_0_z_xyy_yy, to_0_z_xyy_yz, to_0_z_xyy_zz, to_xyy_x, to_xyy_xxz, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyy_xx[k] = 2.0 * to_xyy_xxz[k] * tke_0;

            to_0_z_xyy_xy[k] = 2.0 * to_xyy_xyz[k] * tke_0;

            to_0_z_xyy_xz[k] = -to_xyy_x[k] + 2.0 * to_xyy_xzz[k] * tke_0;

            to_0_z_xyy_yy[k] = 2.0 * to_xyy_yyz[k] * tke_0;

            to_0_z_xyy_yz[k] = -to_xyy_y[k] + 2.0 * to_xyy_yzz[k] * tke_0;

            to_0_z_xyy_zz[k] = -2.0 * to_xyy_z[k] + 2.0 * to_xyy_zzz[k] * tke_0;
        }

        // Set up 144-150 components of targeted buffer : FD

        auto to_0_z_xyz_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 24);

        auto to_0_z_xyz_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 25);

        auto to_0_z_xyz_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 26);

        auto to_0_z_xyz_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 27);

        auto to_0_z_xyz_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 28);

        auto to_0_z_xyz_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_0_z_xyz_xx, to_0_z_xyz_xy, to_0_z_xyz_xz, to_0_z_xyz_yy, to_0_z_xyz_yz, to_0_z_xyz_zz, to_xyz_x, to_xyz_xxz, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyz_xx[k] = 2.0 * to_xyz_xxz[k] * tke_0;

            to_0_z_xyz_xy[k] = 2.0 * to_xyz_xyz[k] * tke_0;

            to_0_z_xyz_xz[k] = -to_xyz_x[k] + 2.0 * to_xyz_xzz[k] * tke_0;

            to_0_z_xyz_yy[k] = 2.0 * to_xyz_yyz[k] * tke_0;

            to_0_z_xyz_yz[k] = -to_xyz_y[k] + 2.0 * to_xyz_yzz[k] * tke_0;

            to_0_z_xyz_zz[k] = -2.0 * to_xyz_z[k] + 2.0 * to_xyz_zzz[k] * tke_0;
        }

        // Set up 150-156 components of targeted buffer : FD

        auto to_0_z_xzz_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 30);

        auto to_0_z_xzz_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 31);

        auto to_0_z_xzz_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 32);

        auto to_0_z_xzz_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 33);

        auto to_0_z_xzz_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 34);

        auto to_0_z_xzz_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_0_z_xzz_xx, to_0_z_xzz_xy, to_0_z_xzz_xz, to_0_z_xzz_yy, to_0_z_xzz_yz, to_0_z_xzz_zz, to_xzz_x, to_xzz_xxz, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzz_xx[k] = 2.0 * to_xzz_xxz[k] * tke_0;

            to_0_z_xzz_xy[k] = 2.0 * to_xzz_xyz[k] * tke_0;

            to_0_z_xzz_xz[k] = -to_xzz_x[k] + 2.0 * to_xzz_xzz[k] * tke_0;

            to_0_z_xzz_yy[k] = 2.0 * to_xzz_yyz[k] * tke_0;

            to_0_z_xzz_yz[k] = -to_xzz_y[k] + 2.0 * to_xzz_yzz[k] * tke_0;

            to_0_z_xzz_zz[k] = -2.0 * to_xzz_z[k] + 2.0 * to_xzz_zzz[k] * tke_0;
        }

        // Set up 156-162 components of targeted buffer : FD

        auto to_0_z_yyy_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 36);

        auto to_0_z_yyy_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 37);

        auto to_0_z_yyy_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 38);

        auto to_0_z_yyy_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 39);

        auto to_0_z_yyy_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 40);

        auto to_0_z_yyy_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_0_z_yyy_xx, to_0_z_yyy_xy, to_0_z_yyy_xz, to_0_z_yyy_yy, to_0_z_yyy_yz, to_0_z_yyy_zz, to_yyy_x, to_yyy_xxz, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_yyz, to_yyy_yzz, to_yyy_z, to_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyy_xx[k] = 2.0 * to_yyy_xxz[k] * tke_0;

            to_0_z_yyy_xy[k] = 2.0 * to_yyy_xyz[k] * tke_0;

            to_0_z_yyy_xz[k] = -to_yyy_x[k] + 2.0 * to_yyy_xzz[k] * tke_0;

            to_0_z_yyy_yy[k] = 2.0 * to_yyy_yyz[k] * tke_0;

            to_0_z_yyy_yz[k] = -to_yyy_y[k] + 2.0 * to_yyy_yzz[k] * tke_0;

            to_0_z_yyy_zz[k] = -2.0 * to_yyy_z[k] + 2.0 * to_yyy_zzz[k] * tke_0;
        }

        // Set up 162-168 components of targeted buffer : FD

        auto to_0_z_yyz_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 42);

        auto to_0_z_yyz_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 43);

        auto to_0_z_yyz_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 44);

        auto to_0_z_yyz_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 45);

        auto to_0_z_yyz_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 46);

        auto to_0_z_yyz_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_0_z_yyz_xx, to_0_z_yyz_xy, to_0_z_yyz_xz, to_0_z_yyz_yy, to_0_z_yyz_yz, to_0_z_yyz_zz, to_yyz_x, to_yyz_xxz, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyz_xx[k] = 2.0 * to_yyz_xxz[k] * tke_0;

            to_0_z_yyz_xy[k] = 2.0 * to_yyz_xyz[k] * tke_0;

            to_0_z_yyz_xz[k] = -to_yyz_x[k] + 2.0 * to_yyz_xzz[k] * tke_0;

            to_0_z_yyz_yy[k] = 2.0 * to_yyz_yyz[k] * tke_0;

            to_0_z_yyz_yz[k] = -to_yyz_y[k] + 2.0 * to_yyz_yzz[k] * tke_0;

            to_0_z_yyz_zz[k] = -2.0 * to_yyz_z[k] + 2.0 * to_yyz_zzz[k] * tke_0;
        }

        // Set up 168-174 components of targeted buffer : FD

        auto to_0_z_yzz_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 48);

        auto to_0_z_yzz_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 49);

        auto to_0_z_yzz_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 50);

        auto to_0_z_yzz_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 51);

        auto to_0_z_yzz_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 52);

        auto to_0_z_yzz_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_0_z_yzz_xx, to_0_z_yzz_xy, to_0_z_yzz_xz, to_0_z_yzz_yy, to_0_z_yzz_yz, to_0_z_yzz_zz, to_yzz_x, to_yzz_xxz, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzz_xx[k] = 2.0 * to_yzz_xxz[k] * tke_0;

            to_0_z_yzz_xy[k] = 2.0 * to_yzz_xyz[k] * tke_0;

            to_0_z_yzz_xz[k] = -to_yzz_x[k] + 2.0 * to_yzz_xzz[k] * tke_0;

            to_0_z_yzz_yy[k] = 2.0 * to_yzz_yyz[k] * tke_0;

            to_0_z_yzz_yz[k] = -to_yzz_y[k] + 2.0 * to_yzz_yzz[k] * tke_0;

            to_0_z_yzz_zz[k] = -2.0 * to_yzz_z[k] + 2.0 * to_yzz_zzz[k] * tke_0;
        }

        // Set up 174-180 components of targeted buffer : FD

        auto to_0_z_zzz_xx = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 54);

        auto to_0_z_zzz_xy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 55);

        auto to_0_z_zzz_xz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 56);

        auto to_0_z_zzz_yy = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 57);

        auto to_0_z_zzz_yz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 58);

        auto to_0_z_zzz_zz = pbuffer.data(idx_op_geom_001_fd + 2 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_0_z_zzz_xx, to_0_z_zzz_xy, to_0_z_zzz_xz, to_0_z_zzz_yy, to_0_z_zzz_yz, to_0_z_zzz_zz, to_zzz_x, to_zzz_xxz, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_yyz, to_zzz_yzz, to_zzz_z, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzz_xx[k] = 2.0 * to_zzz_xxz[k] * tke_0;

            to_0_z_zzz_xy[k] = 2.0 * to_zzz_xyz[k] * tke_0;

            to_0_z_zzz_xz[k] = -to_zzz_x[k] + 2.0 * to_zzz_xzz[k] * tke_0;

            to_0_z_zzz_yy[k] = 2.0 * to_zzz_yyz[k] * tke_0;

            to_0_z_zzz_yz[k] = -to_zzz_y[k] + 2.0 * to_zzz_yzz[k] * tke_0;

            to_0_z_zzz_zz[k] = -2.0 * to_zzz_z[k] + 2.0 * to_zzz_zzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

