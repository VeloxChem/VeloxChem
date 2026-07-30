#include "GeometricalDerivatives020ForFS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_fs(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_fs,
                         const int idx_op_ps,
                         const int idx_op_dp,
                         const int idx_op_fs,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const int idx_op_hs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : PS

    auto tr_x_0 = pbuffer.data(idx_op_ps);

    auto tr_y_0 = pbuffer.data(idx_op_ps + 1);

    auto tr_z_0 = pbuffer.data(idx_op_ps + 2);

    // Set up components of auxiliary buffer : DP

    auto tr_xx_x = pbuffer.data(idx_op_dp);

    auto tr_xx_y = pbuffer.data(idx_op_dp + 1);

    auto tr_xx_z = pbuffer.data(idx_op_dp + 2);

    auto tr_xy_x = pbuffer.data(idx_op_dp + 3);

    auto tr_xy_y = pbuffer.data(idx_op_dp + 4);

    auto tr_xy_z = pbuffer.data(idx_op_dp + 5);

    auto tr_xz_x = pbuffer.data(idx_op_dp + 6);

    auto tr_xz_y = pbuffer.data(idx_op_dp + 7);

    auto tr_xz_z = pbuffer.data(idx_op_dp + 8);

    auto tr_yy_x = pbuffer.data(idx_op_dp + 9);

    auto tr_yy_y = pbuffer.data(idx_op_dp + 10);

    auto tr_yy_z = pbuffer.data(idx_op_dp + 11);

    auto tr_yz_x = pbuffer.data(idx_op_dp + 12);

    auto tr_yz_y = pbuffer.data(idx_op_dp + 13);

    auto tr_yz_z = pbuffer.data(idx_op_dp + 14);

    auto tr_zz_x = pbuffer.data(idx_op_dp + 15);

    auto tr_zz_y = pbuffer.data(idx_op_dp + 16);

    auto tr_zz_z = pbuffer.data(idx_op_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto tr_xxx_0 = pbuffer.data(idx_op_fs);

    auto tr_xxy_0 = pbuffer.data(idx_op_fs + 1);

    auto tr_xxz_0 = pbuffer.data(idx_op_fs + 2);

    auto tr_xyy_0 = pbuffer.data(idx_op_fs + 3);

    auto tr_xyz_0 = pbuffer.data(idx_op_fs + 4);

    auto tr_xzz_0 = pbuffer.data(idx_op_fs + 5);

    auto tr_yyy_0 = pbuffer.data(idx_op_fs + 6);

    auto tr_yyz_0 = pbuffer.data(idx_op_fs + 7);

    auto tr_yzz_0 = pbuffer.data(idx_op_fs + 8);

    auto tr_zzz_0 = pbuffer.data(idx_op_fs + 9);

    // Set up components of auxiliary buffer : FD

    auto tr_xxx_xx = pbuffer.data(idx_op_fd);

    auto tr_xxx_xy = pbuffer.data(idx_op_fd + 1);

    auto tr_xxx_xz = pbuffer.data(idx_op_fd + 2);

    auto tr_xxx_yy = pbuffer.data(idx_op_fd + 3);

    auto tr_xxx_yz = pbuffer.data(idx_op_fd + 4);

    auto tr_xxx_zz = pbuffer.data(idx_op_fd + 5);

    auto tr_xxy_xx = pbuffer.data(idx_op_fd + 6);

    auto tr_xxy_xy = pbuffer.data(idx_op_fd + 7);

    auto tr_xxy_xz = pbuffer.data(idx_op_fd + 8);

    auto tr_xxy_yy = pbuffer.data(idx_op_fd + 9);

    auto tr_xxy_yz = pbuffer.data(idx_op_fd + 10);

    auto tr_xxy_zz = pbuffer.data(idx_op_fd + 11);

    auto tr_xxz_xx = pbuffer.data(idx_op_fd + 12);

    auto tr_xxz_xy = pbuffer.data(idx_op_fd + 13);

    auto tr_xxz_xz = pbuffer.data(idx_op_fd + 14);

    auto tr_xxz_yy = pbuffer.data(idx_op_fd + 15);

    auto tr_xxz_yz = pbuffer.data(idx_op_fd + 16);

    auto tr_xxz_zz = pbuffer.data(idx_op_fd + 17);

    auto tr_xyy_xx = pbuffer.data(idx_op_fd + 18);

    auto tr_xyy_xy = pbuffer.data(idx_op_fd + 19);

    auto tr_xyy_xz = pbuffer.data(idx_op_fd + 20);

    auto tr_xyy_yy = pbuffer.data(idx_op_fd + 21);

    auto tr_xyy_yz = pbuffer.data(idx_op_fd + 22);

    auto tr_xyy_zz = pbuffer.data(idx_op_fd + 23);

    auto tr_xyz_xx = pbuffer.data(idx_op_fd + 24);

    auto tr_xyz_xy = pbuffer.data(idx_op_fd + 25);

    auto tr_xyz_xz = pbuffer.data(idx_op_fd + 26);

    auto tr_xyz_yy = pbuffer.data(idx_op_fd + 27);

    auto tr_xyz_yz = pbuffer.data(idx_op_fd + 28);

    auto tr_xyz_zz = pbuffer.data(idx_op_fd + 29);

    auto tr_xzz_xx = pbuffer.data(idx_op_fd + 30);

    auto tr_xzz_xy = pbuffer.data(idx_op_fd + 31);

    auto tr_xzz_xz = pbuffer.data(idx_op_fd + 32);

    auto tr_xzz_yy = pbuffer.data(idx_op_fd + 33);

    auto tr_xzz_yz = pbuffer.data(idx_op_fd + 34);

    auto tr_xzz_zz = pbuffer.data(idx_op_fd + 35);

    auto tr_yyy_xx = pbuffer.data(idx_op_fd + 36);

    auto tr_yyy_xy = pbuffer.data(idx_op_fd + 37);

    auto tr_yyy_xz = pbuffer.data(idx_op_fd + 38);

    auto tr_yyy_yy = pbuffer.data(idx_op_fd + 39);

    auto tr_yyy_yz = pbuffer.data(idx_op_fd + 40);

    auto tr_yyy_zz = pbuffer.data(idx_op_fd + 41);

    auto tr_yyz_xx = pbuffer.data(idx_op_fd + 42);

    auto tr_yyz_xy = pbuffer.data(idx_op_fd + 43);

    auto tr_yyz_xz = pbuffer.data(idx_op_fd + 44);

    auto tr_yyz_yy = pbuffer.data(idx_op_fd + 45);

    auto tr_yyz_yz = pbuffer.data(idx_op_fd + 46);

    auto tr_yyz_zz = pbuffer.data(idx_op_fd + 47);

    auto tr_yzz_xx = pbuffer.data(idx_op_fd + 48);

    auto tr_yzz_xy = pbuffer.data(idx_op_fd + 49);

    auto tr_yzz_xz = pbuffer.data(idx_op_fd + 50);

    auto tr_yzz_yy = pbuffer.data(idx_op_fd + 51);

    auto tr_yzz_yz = pbuffer.data(idx_op_fd + 52);

    auto tr_yzz_zz = pbuffer.data(idx_op_fd + 53);

    auto tr_zzz_xx = pbuffer.data(idx_op_fd + 54);

    auto tr_zzz_xy = pbuffer.data(idx_op_fd + 55);

    auto tr_zzz_xz = pbuffer.data(idx_op_fd + 56);

    auto tr_zzz_yy = pbuffer.data(idx_op_fd + 57);

    auto tr_zzz_yz = pbuffer.data(idx_op_fd + 58);

    auto tr_zzz_zz = pbuffer.data(idx_op_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto tr_xxxx_x = pbuffer.data(idx_op_gp);

    auto tr_xxxx_y = pbuffer.data(idx_op_gp + 1);

    auto tr_xxxx_z = pbuffer.data(idx_op_gp + 2);

    auto tr_xxxy_x = pbuffer.data(idx_op_gp + 3);

    auto tr_xxxy_y = pbuffer.data(idx_op_gp + 4);

    auto tr_xxxy_z = pbuffer.data(idx_op_gp + 5);

    auto tr_xxxz_x = pbuffer.data(idx_op_gp + 6);

    auto tr_xxxz_y = pbuffer.data(idx_op_gp + 7);

    auto tr_xxxz_z = pbuffer.data(idx_op_gp + 8);

    auto tr_xxyy_x = pbuffer.data(idx_op_gp + 9);

    auto tr_xxyy_y = pbuffer.data(idx_op_gp + 10);

    auto tr_xxyy_z = pbuffer.data(idx_op_gp + 11);

    auto tr_xxyz_x = pbuffer.data(idx_op_gp + 12);

    auto tr_xxyz_y = pbuffer.data(idx_op_gp + 13);

    auto tr_xxyz_z = pbuffer.data(idx_op_gp + 14);

    auto tr_xxzz_x = pbuffer.data(idx_op_gp + 15);

    auto tr_xxzz_y = pbuffer.data(idx_op_gp + 16);

    auto tr_xxzz_z = pbuffer.data(idx_op_gp + 17);

    auto tr_xyyy_x = pbuffer.data(idx_op_gp + 18);

    auto tr_xyyy_y = pbuffer.data(idx_op_gp + 19);

    auto tr_xyyy_z = pbuffer.data(idx_op_gp + 20);

    auto tr_xyyz_x = pbuffer.data(idx_op_gp + 21);

    auto tr_xyyz_y = pbuffer.data(idx_op_gp + 22);

    auto tr_xyyz_z = pbuffer.data(idx_op_gp + 23);

    auto tr_xyzz_x = pbuffer.data(idx_op_gp + 24);

    auto tr_xyzz_y = pbuffer.data(idx_op_gp + 25);

    auto tr_xyzz_z = pbuffer.data(idx_op_gp + 26);

    auto tr_xzzz_x = pbuffer.data(idx_op_gp + 27);

    auto tr_xzzz_y = pbuffer.data(idx_op_gp + 28);

    auto tr_xzzz_z = pbuffer.data(idx_op_gp + 29);

    auto tr_yyyy_x = pbuffer.data(idx_op_gp + 30);

    auto tr_yyyy_y = pbuffer.data(idx_op_gp + 31);

    auto tr_yyyy_z = pbuffer.data(idx_op_gp + 32);

    auto tr_yyyz_x = pbuffer.data(idx_op_gp + 33);

    auto tr_yyyz_y = pbuffer.data(idx_op_gp + 34);

    auto tr_yyyz_z = pbuffer.data(idx_op_gp + 35);

    auto tr_yyzz_x = pbuffer.data(idx_op_gp + 36);

    auto tr_yyzz_y = pbuffer.data(idx_op_gp + 37);

    auto tr_yyzz_z = pbuffer.data(idx_op_gp + 38);

    auto tr_yzzz_x = pbuffer.data(idx_op_gp + 39);

    auto tr_yzzz_y = pbuffer.data(idx_op_gp + 40);

    auto tr_yzzz_z = pbuffer.data(idx_op_gp + 41);

    auto tr_zzzz_x = pbuffer.data(idx_op_gp + 42);

    auto tr_zzzz_y = pbuffer.data(idx_op_gp + 43);

    auto tr_zzzz_z = pbuffer.data(idx_op_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto tr_xxxxx_0 = pbuffer.data(idx_op_hs);

    auto tr_xxxxy_0 = pbuffer.data(idx_op_hs + 1);

    auto tr_xxxxz_0 = pbuffer.data(idx_op_hs + 2);

    auto tr_xxxyy_0 = pbuffer.data(idx_op_hs + 3);

    auto tr_xxxyz_0 = pbuffer.data(idx_op_hs + 4);

    auto tr_xxxzz_0 = pbuffer.data(idx_op_hs + 5);

    auto tr_xxyyy_0 = pbuffer.data(idx_op_hs + 6);

    auto tr_xxyyz_0 = pbuffer.data(idx_op_hs + 7);

    auto tr_xxyzz_0 = pbuffer.data(idx_op_hs + 8);

    auto tr_xxzzz_0 = pbuffer.data(idx_op_hs + 9);

    auto tr_xyyyy_0 = pbuffer.data(idx_op_hs + 10);

    auto tr_xyyyz_0 = pbuffer.data(idx_op_hs + 11);

    auto tr_xyyzz_0 = pbuffer.data(idx_op_hs + 12);

    auto tr_xyzzz_0 = pbuffer.data(idx_op_hs + 13);

    auto tr_xzzzz_0 = pbuffer.data(idx_op_hs + 14);

    auto tr_yyyyy_0 = pbuffer.data(idx_op_hs + 15);

    auto tr_yyyyz_0 = pbuffer.data(idx_op_hs + 16);

    auto tr_yyyzz_0 = pbuffer.data(idx_op_hs + 17);

    auto tr_yyzzz_0 = pbuffer.data(idx_op_hs + 18);

    auto tr_yzzzz_0 = pbuffer.data(idx_op_hs + 19);

    auto tr_zzzzz_0 = pbuffer.data(idx_op_hs + 20);

    // Set up components of targeted buffer : FS

    auto tr_0_0_xx_xxx_0 = pbuffer.data(idx_op_geom_020_fs);

    auto tr_0_0_xx_xxy_0 = pbuffer.data(idx_op_geom_020_fs + 1);

    auto tr_0_0_xx_xxz_0 = pbuffer.data(idx_op_geom_020_fs + 2);

    auto tr_0_0_xx_xyy_0 = pbuffer.data(idx_op_geom_020_fs + 3);

    auto tr_0_0_xx_xyz_0 = pbuffer.data(idx_op_geom_020_fs + 4);

    auto tr_0_0_xx_xzz_0 = pbuffer.data(idx_op_geom_020_fs + 5);

    auto tr_0_0_xx_yyy_0 = pbuffer.data(idx_op_geom_020_fs + 6);

    auto tr_0_0_xx_yyz_0 = pbuffer.data(idx_op_geom_020_fs + 7);

    auto tr_0_0_xx_yzz_0 = pbuffer.data(idx_op_geom_020_fs + 8);

    auto tr_0_0_xx_zzz_0 = pbuffer.data(idx_op_geom_020_fs + 9);

    auto tr_0_0_xy_xxx_0 = pbuffer.data(idx_op_geom_020_fs + 10);

    auto tr_0_0_xy_xxy_0 = pbuffer.data(idx_op_geom_020_fs + 11);

    auto tr_0_0_xy_xxz_0 = pbuffer.data(idx_op_geom_020_fs + 12);

    auto tr_0_0_xy_xyy_0 = pbuffer.data(idx_op_geom_020_fs + 13);

    auto tr_0_0_xy_xyz_0 = pbuffer.data(idx_op_geom_020_fs + 14);

    auto tr_0_0_xy_xzz_0 = pbuffer.data(idx_op_geom_020_fs + 15);

    auto tr_0_0_xy_yyy_0 = pbuffer.data(idx_op_geom_020_fs + 16);

    auto tr_0_0_xy_yyz_0 = pbuffer.data(idx_op_geom_020_fs + 17);

    auto tr_0_0_xy_yzz_0 = pbuffer.data(idx_op_geom_020_fs + 18);

    auto tr_0_0_xy_zzz_0 = pbuffer.data(idx_op_geom_020_fs + 19);

    auto tr_0_0_xz_xxx_0 = pbuffer.data(idx_op_geom_020_fs + 20);

    auto tr_0_0_xz_xxy_0 = pbuffer.data(idx_op_geom_020_fs + 21);

    auto tr_0_0_xz_xxz_0 = pbuffer.data(idx_op_geom_020_fs + 22);

    auto tr_0_0_xz_xyy_0 = pbuffer.data(idx_op_geom_020_fs + 23);

    auto tr_0_0_xz_xyz_0 = pbuffer.data(idx_op_geom_020_fs + 24);

    auto tr_0_0_xz_xzz_0 = pbuffer.data(idx_op_geom_020_fs + 25);

    auto tr_0_0_xz_yyy_0 = pbuffer.data(idx_op_geom_020_fs + 26);

    auto tr_0_0_xz_yyz_0 = pbuffer.data(idx_op_geom_020_fs + 27);

    auto tr_0_0_xz_yzz_0 = pbuffer.data(idx_op_geom_020_fs + 28);

    auto tr_0_0_xz_zzz_0 = pbuffer.data(idx_op_geom_020_fs + 29);

    auto tr_0_0_yy_xxx_0 = pbuffer.data(idx_op_geom_020_fs + 30);

    auto tr_0_0_yy_xxy_0 = pbuffer.data(idx_op_geom_020_fs + 31);

    auto tr_0_0_yy_xxz_0 = pbuffer.data(idx_op_geom_020_fs + 32);

    auto tr_0_0_yy_xyy_0 = pbuffer.data(idx_op_geom_020_fs + 33);

    auto tr_0_0_yy_xyz_0 = pbuffer.data(idx_op_geom_020_fs + 34);

    auto tr_0_0_yy_xzz_0 = pbuffer.data(idx_op_geom_020_fs + 35);

    auto tr_0_0_yy_yyy_0 = pbuffer.data(idx_op_geom_020_fs + 36);

    auto tr_0_0_yy_yyz_0 = pbuffer.data(idx_op_geom_020_fs + 37);

    auto tr_0_0_yy_yzz_0 = pbuffer.data(idx_op_geom_020_fs + 38);

    auto tr_0_0_yy_zzz_0 = pbuffer.data(idx_op_geom_020_fs + 39);

    auto tr_0_0_yz_xxx_0 = pbuffer.data(idx_op_geom_020_fs + 40);

    auto tr_0_0_yz_xxy_0 = pbuffer.data(idx_op_geom_020_fs + 41);

    auto tr_0_0_yz_xxz_0 = pbuffer.data(idx_op_geom_020_fs + 42);

    auto tr_0_0_yz_xyy_0 = pbuffer.data(idx_op_geom_020_fs + 43);

    auto tr_0_0_yz_xyz_0 = pbuffer.data(idx_op_geom_020_fs + 44);

    auto tr_0_0_yz_xzz_0 = pbuffer.data(idx_op_geom_020_fs + 45);

    auto tr_0_0_yz_yyy_0 = pbuffer.data(idx_op_geom_020_fs + 46);

    auto tr_0_0_yz_yyz_0 = pbuffer.data(idx_op_geom_020_fs + 47);

    auto tr_0_0_yz_yzz_0 = pbuffer.data(idx_op_geom_020_fs + 48);

    auto tr_0_0_yz_zzz_0 = pbuffer.data(idx_op_geom_020_fs + 49);

    auto tr_0_0_zz_xxx_0 = pbuffer.data(idx_op_geom_020_fs + 50);

    auto tr_0_0_zz_xxy_0 = pbuffer.data(idx_op_geom_020_fs + 51);

    auto tr_0_0_zz_xxz_0 = pbuffer.data(idx_op_geom_020_fs + 52);

    auto tr_0_0_zz_xyy_0 = pbuffer.data(idx_op_geom_020_fs + 53);

    auto tr_0_0_zz_xyz_0 = pbuffer.data(idx_op_geom_020_fs + 54);

    auto tr_0_0_zz_xzz_0 = pbuffer.data(idx_op_geom_020_fs + 55);

    auto tr_0_0_zz_yyy_0 = pbuffer.data(idx_op_geom_020_fs + 56);

    auto tr_0_0_zz_yyz_0 = pbuffer.data(idx_op_geom_020_fs + 57);

    auto tr_0_0_zz_yzz_0 = pbuffer.data(idx_op_geom_020_fs + 58);

    auto tr_0_0_zz_zzz_0 = pbuffer.data(idx_op_geom_020_fs + 59);

    #pragma omp simd aligned(tr_0_0_xx_xxx_0, tr_0_0_xx_xxy_0, tr_0_0_xx_xxz_0, tr_0_0_xx_xyy_0, tr_0_0_xx_xyz_0, tr_0_0_xx_xzz_0, tr_0_0_xx_yyy_0, tr_0_0_xx_yyz_0, tr_0_0_xx_yzz_0, tr_0_0_xx_zzz_0, tr_0_0_xy_xxx_0, tr_0_0_xy_xxy_0, tr_0_0_xy_xxz_0, tr_0_0_xy_xyy_0, tr_0_0_xy_xyz_0, tr_0_0_xy_xzz_0, tr_0_0_xy_yyy_0, tr_0_0_xy_yyz_0, tr_0_0_xy_yzz_0, tr_0_0_xy_zzz_0, tr_0_0_xz_xxx_0, tr_0_0_xz_xxy_0, tr_0_0_xz_xxz_0, tr_0_0_xz_xyy_0, tr_0_0_xz_xyz_0, tr_0_0_xz_xzz_0, tr_0_0_xz_yyy_0, tr_0_0_xz_yyz_0, tr_0_0_xz_yzz_0, tr_0_0_xz_zzz_0, tr_0_0_yy_xxx_0, tr_0_0_yy_xxy_0, tr_0_0_yy_xxz_0, tr_0_0_yy_xyy_0, tr_0_0_yy_xyz_0, tr_0_0_yy_xzz_0, tr_0_0_yy_yyy_0, tr_0_0_yy_yyz_0, tr_0_0_yy_yzz_0, tr_0_0_yy_zzz_0, tr_0_0_yz_xxx_0, tr_0_0_yz_xxy_0, tr_0_0_yz_xxz_0, tr_0_0_yz_xyy_0, tr_0_0_yz_xyz_0, tr_0_0_yz_xzz_0, tr_0_0_yz_yyy_0, tr_0_0_yz_yyz_0, tr_0_0_yz_yzz_0, tr_0_0_yz_zzz_0, tr_0_0_zz_xxx_0, tr_0_0_zz_xxy_0, tr_0_0_zz_xxz_0, tr_0_0_zz_xyy_0, tr_0_0_zz_xyz_0, tr_0_0_zz_xzz_0, tr_0_0_zz_yyy_0, tr_0_0_zz_yyz_0, tr_0_0_zz_yzz_0, tr_0_0_zz_zzz_0, tr_x_0, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxx_0, tr_xxxxy_0, tr_xxxxz_0, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyy_0, tr_xxxyz_0, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxxzz_0, tr_xxy_0, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyy_0, tr_xxyyz_0, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxyzz_0, tr_xxz_0, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xxzzz_0, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyy_0, tr_xyyyz_0, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyyzz_0, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xyzzz_0, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_xzzzz_0, tr_y_0, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, tr_yyyyy_0, tr_yyyyz_0, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyyzz_0, tr_yyz_0, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yyzzz_0, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_yzzzz_0, tr_z_0, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, tr_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxx_0[i] = 6.0 * tr_x_0[i] - 12.0 * tr_xx_x[i] * tke_0 - 14.0 * tr_xxx_0[i] * tbe_0 - 2.0 * tr_xxx_0[i] * tke_0 + 4.0 * tr_xxx_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxx_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_0[i] = 2.0 * tr_y_0[i] - 8.0 * tr_xy_x[i] * tke_0 - 10.0 * tr_xxy_0[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tke_0 + 4.0 * tr_xxy_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_0[i] = 2.0 * tr_z_0[i] - 8.0 * tr_xz_x[i] * tke_0 - 10.0 * tr_xxz_0[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tke_0 + 4.0 * tr_xxz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_0[i] = -4.0 * tr_yy_x[i] * tke_0 - 6.0 * tr_xyy_0[i] * tbe_0 - 2.0 * tr_xyy_0[i] * tke_0 + 4.0 * tr_xyy_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_0[i] = -4.0 * tr_yz_x[i] * tke_0 - 6.0 * tr_xyz_0[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tke_0 + 4.0 * tr_xyz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_0[i] = -4.0 * tr_zz_x[i] * tke_0 - 6.0 * tr_xzz_0[i] * tbe_0 - 2.0 * tr_xzz_0[i] * tke_0 + 4.0 * tr_xzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_0[i] = -2.0 * tr_yyy_0[i] * tbe_0 - 2.0 * tr_yyy_0[i] * tke_0 + 4.0 * tr_yyy_xx[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_0[i] = -2.0 * tr_yyz_0[i] * tbe_0 - 2.0 * tr_yyz_0[i] * tke_0 + 4.0 * tr_yyz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_0[i] = -2.0 * tr_yzz_0[i] * tbe_0 - 2.0 * tr_yzz_0[i] * tke_0 + 4.0 * tr_yzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_0[i] = -2.0 * tr_zzz_0[i] * tbe_0 - 2.0 * tr_zzz_0[i] * tke_0 + 4.0 * tr_zzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_0[i] = -6.0 * tr_xx_y[i] * tke_0 - 6.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxx_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxxy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_0[i] = 2.0 * tr_x_0[i] - 4.0 * tr_xy_y[i] * tke_0 - 4.0 * tr_xyy_0[i] * tbe_0 - 2.0 * tr_xx_x[i] * tke_0 + 4.0 * tr_xxy_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxyy_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_0[i] * tbe_0 + 4.0 * tr_xxxy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_0[i] = -4.0 * tr_xz_y[i] * tke_0 - 4.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xxz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_0[i] = 2.0 * tr_y_0[i] - 2.0 * tr_yy_y[i] * tke_0 - 2.0 * tr_yyy_0[i] * tbe_0 - 4.0 * tr_xy_x[i] * tke_0 + 4.0 * tr_xyy_xy[i] * tke_0 * tke_0 + 4.0 * tr_xyyy_x[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_0[i] = tr_z_0[i] - 2.0 * tr_yz_y[i] * tke_0 - 2.0 * tr_yyz_0[i] * tbe_0 - 2.0 * tr_xz_x[i] * tke_0 + 4.0 * tr_xyz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_0[i] = -2.0 * tr_zz_y[i] * tke_0 - 2.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_xzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_0[i] = -6.0 * tr_yy_x[i] * tke_0 + 4.0 * tr_yyy_xy[i] * tke_0 * tke_0 + 4.0 * tr_yyyy_x[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_0[i] = -4.0 * tr_yz_x[i] * tke_0 + 4.0 * tr_yyz_xy[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_0[i] = -2.0 * tr_zz_x[i] * tke_0 + 4.0 * tr_yzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_0[i] = 4.0 * tr_zzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_0[i] = -6.0 * tr_xx_z[i] * tke_0 - 6.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxx_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_0[i] = -4.0 * tr_xy_z[i] * tke_0 - 4.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xxy_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_0[i] = 2.0 * tr_x_0[i] - 4.0 * tr_xz_z[i] * tke_0 - 4.0 * tr_xzz_0[i] * tbe_0 - 2.0 * tr_xx_x[i] * tke_0 + 4.0 * tr_xxz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_0[i] * tbe_0 + 4.0 * tr_xxxz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_0[i] = -2.0 * tr_yy_z[i] * tke_0 - 2.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_xyy_xz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_0[i] = tr_y_0[i] - 2.0 * tr_yz_z[i] * tke_0 - 2.0 * tr_yzz_0[i] * tbe_0 - 2.0 * tr_xy_x[i] * tke_0 + 4.0 * tr_xyz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_0[i] = 2.0 * tr_z_0[i] - 2.0 * tr_zz_z[i] * tke_0 - 2.0 * tr_zzz_0[i] * tbe_0 - 4.0 * tr_xz_x[i] * tke_0 + 4.0 * tr_xzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_0[i] = 4.0 * tr_yyy_xz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_0[i] = -2.0 * tr_yy_x[i] * tke_0 + 4.0 * tr_yyz_xz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_0[i] = -4.0 * tr_yz_x[i] * tke_0 + 4.0 * tr_yzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_0[i] = -6.0 * tr_zz_x[i] * tke_0 + 4.0 * tr_zzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_x[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_0[i] = -2.0 * tr_xxx_0[i] * tbe_0 - 2.0 * tr_xxx_0[i] * tke_0 + 4.0 * tr_xxx_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_0[i] = -4.0 * tr_xx_y[i] * tke_0 - 6.0 * tr_xxy_0[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tke_0 + 4.0 * tr_xxy_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_0[i] = -2.0 * tr_xxz_0[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tke_0 + 4.0 * tr_xxz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_0[i] = 2.0 * tr_x_0[i] - 8.0 * tr_xy_y[i] * tke_0 - 10.0 * tr_xyy_0[i] * tbe_0 - 2.0 * tr_xyy_0[i] * tke_0 + 4.0 * tr_xyy_yy[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_0[i] = -4.0 * tr_xz_y[i] * tke_0 - 6.0 * tr_xyz_0[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tke_0 + 4.0 * tr_xyz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_0[i] = -2.0 * tr_xzz_0[i] * tbe_0 - 2.0 * tr_xzz_0[i] * tke_0 + 4.0 * tr_xzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_0[i] = 6.0 * tr_y_0[i] - 12.0 * tr_yy_y[i] * tke_0 - 14.0 * tr_yyy_0[i] * tbe_0 - 2.0 * tr_yyy_0[i] * tke_0 + 4.0 * tr_yyy_yy[i] * tke_0 * tke_0 + 8.0 * tr_yyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_0[i] = 2.0 * tr_z_0[i] - 8.0 * tr_yz_y[i] * tke_0 - 10.0 * tr_yyz_0[i] * tbe_0 - 2.0 * tr_yyz_0[i] * tke_0 + 4.0 * tr_yyz_yy[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_0[i] = -4.0 * tr_zz_y[i] * tke_0 - 6.0 * tr_yzz_0[i] * tbe_0 - 2.0 * tr_yzz_0[i] * tke_0 + 4.0 * tr_yzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_0[i] = -2.0 * tr_zzz_0[i] * tbe_0 - 2.0 * tr_zzz_0[i] * tke_0 + 4.0 * tr_zzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_0[i] = 4.0 * tr_xxx_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_0[i] = -2.0 * tr_xx_z[i] * tke_0 - 2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxy_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_0[i] = -2.0 * tr_xx_y[i] * tke_0 + 4.0 * tr_xxz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_0[i] = -4.0 * tr_xy_z[i] * tke_0 - 4.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyy_yz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_0[i] = tr_x_0[i] - 2.0 * tr_xz_z[i] * tke_0 - 2.0 * tr_xzz_0[i] * tbe_0 - 2.0 * tr_xy_y[i] * tke_0 + 4.0 * tr_xyz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_0[i] = -4.0 * tr_xz_y[i] * tke_0 + 4.0 * tr_xzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_y[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_0[i] = -6.0 * tr_yy_z[i] * tke_0 - 6.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyy_yz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_0[i] = 2.0 * tr_y_0[i] - 4.0 * tr_yz_z[i] * tke_0 - 4.0 * tr_yzz_0[i] * tbe_0 - 2.0 * tr_yy_y[i] * tke_0 + 4.0 * tr_yyz_yz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_0[i] * tbe_0 + 4.0 * tr_yyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_0[i] = 2.0 * tr_z_0[i] - 2.0 * tr_zz_z[i] * tke_0 - 2.0 * tr_zzz_0[i] * tbe_0 - 4.0 * tr_yz_y[i] * tke_0 + 4.0 * tr_yzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_y[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_0[i] = -6.0 * tr_zz_y[i] * tke_0 + 4.0 * tr_zzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_y[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_yzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_0[i] = -2.0 * tr_xxx_0[i] * tbe_0 - 2.0 * tr_xxx_0[i] * tke_0 + 4.0 * tr_xxx_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_0[i] = -2.0 * tr_xxy_0[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tke_0 + 4.0 * tr_xxy_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_0[i] = -4.0 * tr_xx_z[i] * tke_0 - 6.0 * tr_xxz_0[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tke_0 + 4.0 * tr_xxz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_0[i] = -2.0 * tr_xyy_0[i] * tbe_0 - 2.0 * tr_xyy_0[i] * tke_0 + 4.0 * tr_xyy_zz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_0[i] = -4.0 * tr_xy_z[i] * tke_0 - 6.0 * tr_xyz_0[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tke_0 + 4.0 * tr_xyz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_0[i] = 2.0 * tr_x_0[i] - 8.0 * tr_xz_z[i] * tke_0 - 10.0 * tr_xzz_0[i] * tbe_0 - 2.0 * tr_xzz_0[i] * tke_0 + 4.0 * tr_xzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_0[i] = -2.0 * tr_yyy_0[i] * tbe_0 - 2.0 * tr_yyy_0[i] * tke_0 + 4.0 * tr_yyy_zz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_0[i] = -4.0 * tr_yy_z[i] * tke_0 - 6.0 * tr_yyz_0[i] * tbe_0 - 2.0 * tr_yyz_0[i] * tke_0 + 4.0 * tr_yyz_zz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_0[i] = 2.0 * tr_y_0[i] - 8.0 * tr_yz_z[i] * tke_0 - 10.0 * tr_yzz_0[i] * tbe_0 - 2.0 * tr_yzz_0[i] * tke_0 + 4.0 * tr_yzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_0[i] = 6.0 * tr_z_0[i] - 12.0 * tr_zz_z[i] * tke_0 - 14.0 * tr_zzz_0[i] * tbe_0 - 2.0 * tr_zzz_0[i] * tke_0 + 4.0 * tr_zzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_zzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_0[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

