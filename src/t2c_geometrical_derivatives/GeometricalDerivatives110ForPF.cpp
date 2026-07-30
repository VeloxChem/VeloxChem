#include "GeometricalDerivatives110ForPF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_pf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_pf,
                         const int idx_op_sd,
                         const int idx_op_sg,
                         const int idx_op_pf,
                         const int idx_op_dd,
                         const int idx_op_dg,
                         const int idx_op_ff,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SD

    auto tr_0_xx = pbuffer.data(idx_op_sd);

    auto tr_0_xy = pbuffer.data(idx_op_sd + 1);

    auto tr_0_xz = pbuffer.data(idx_op_sd + 2);

    auto tr_0_yy = pbuffer.data(idx_op_sd + 3);

    auto tr_0_yz = pbuffer.data(idx_op_sd + 4);

    auto tr_0_zz = pbuffer.data(idx_op_sd + 5);

    // Set up components of auxiliary buffer : SG

    auto tr_0_xxxx = pbuffer.data(idx_op_sg);

    auto tr_0_xxxy = pbuffer.data(idx_op_sg + 1);

    auto tr_0_xxxz = pbuffer.data(idx_op_sg + 2);

    auto tr_0_xxyy = pbuffer.data(idx_op_sg + 3);

    auto tr_0_xxyz = pbuffer.data(idx_op_sg + 4);

    auto tr_0_xxzz = pbuffer.data(idx_op_sg + 5);

    auto tr_0_xyyy = pbuffer.data(idx_op_sg + 6);

    auto tr_0_xyyz = pbuffer.data(idx_op_sg + 7);

    auto tr_0_xyzz = pbuffer.data(idx_op_sg + 8);

    auto tr_0_xzzz = pbuffer.data(idx_op_sg + 9);

    auto tr_0_yyyy = pbuffer.data(idx_op_sg + 10);

    auto tr_0_yyyz = pbuffer.data(idx_op_sg + 11);

    auto tr_0_yyzz = pbuffer.data(idx_op_sg + 12);

    auto tr_0_yzzz = pbuffer.data(idx_op_sg + 13);

    auto tr_0_zzzz = pbuffer.data(idx_op_sg + 14);

    // Set up components of auxiliary buffer : PF

    auto tr_x_xxx = pbuffer.data(idx_op_pf);

    auto tr_x_xxy = pbuffer.data(idx_op_pf + 1);

    auto tr_x_xxz = pbuffer.data(idx_op_pf + 2);

    auto tr_x_xyy = pbuffer.data(idx_op_pf + 3);

    auto tr_x_xyz = pbuffer.data(idx_op_pf + 4);

    auto tr_x_xzz = pbuffer.data(idx_op_pf + 5);

    auto tr_x_yyy = pbuffer.data(idx_op_pf + 6);

    auto tr_x_yyz = pbuffer.data(idx_op_pf + 7);

    auto tr_x_yzz = pbuffer.data(idx_op_pf + 8);

    auto tr_x_zzz = pbuffer.data(idx_op_pf + 9);

    auto tr_y_xxx = pbuffer.data(idx_op_pf + 10);

    auto tr_y_xxy = pbuffer.data(idx_op_pf + 11);

    auto tr_y_xxz = pbuffer.data(idx_op_pf + 12);

    auto tr_y_xyy = pbuffer.data(idx_op_pf + 13);

    auto tr_y_xyz = pbuffer.data(idx_op_pf + 14);

    auto tr_y_xzz = pbuffer.data(idx_op_pf + 15);

    auto tr_y_yyy = pbuffer.data(idx_op_pf + 16);

    auto tr_y_yyz = pbuffer.data(idx_op_pf + 17);

    auto tr_y_yzz = pbuffer.data(idx_op_pf + 18);

    auto tr_y_zzz = pbuffer.data(idx_op_pf + 19);

    auto tr_z_xxx = pbuffer.data(idx_op_pf + 20);

    auto tr_z_xxy = pbuffer.data(idx_op_pf + 21);

    auto tr_z_xxz = pbuffer.data(idx_op_pf + 22);

    auto tr_z_xyy = pbuffer.data(idx_op_pf + 23);

    auto tr_z_xyz = pbuffer.data(idx_op_pf + 24);

    auto tr_z_xzz = pbuffer.data(idx_op_pf + 25);

    auto tr_z_yyy = pbuffer.data(idx_op_pf + 26);

    auto tr_z_yyz = pbuffer.data(idx_op_pf + 27);

    auto tr_z_yzz = pbuffer.data(idx_op_pf + 28);

    auto tr_z_zzz = pbuffer.data(idx_op_pf + 29);

    // Set up components of auxiliary buffer : DD

    auto tr_xx_xx = pbuffer.data(idx_op_dd);

    auto tr_xx_xy = pbuffer.data(idx_op_dd + 1);

    auto tr_xx_xz = pbuffer.data(idx_op_dd + 2);

    auto tr_xx_yy = pbuffer.data(idx_op_dd + 3);

    auto tr_xx_yz = pbuffer.data(idx_op_dd + 4);

    auto tr_xx_zz = pbuffer.data(idx_op_dd + 5);

    auto tr_xy_xx = pbuffer.data(idx_op_dd + 6);

    auto tr_xy_xy = pbuffer.data(idx_op_dd + 7);

    auto tr_xy_xz = pbuffer.data(idx_op_dd + 8);

    auto tr_xy_yy = pbuffer.data(idx_op_dd + 9);

    auto tr_xy_yz = pbuffer.data(idx_op_dd + 10);

    auto tr_xy_zz = pbuffer.data(idx_op_dd + 11);

    auto tr_xz_xx = pbuffer.data(idx_op_dd + 12);

    auto tr_xz_xy = pbuffer.data(idx_op_dd + 13);

    auto tr_xz_xz = pbuffer.data(idx_op_dd + 14);

    auto tr_xz_yy = pbuffer.data(idx_op_dd + 15);

    auto tr_xz_yz = pbuffer.data(idx_op_dd + 16);

    auto tr_xz_zz = pbuffer.data(idx_op_dd + 17);

    auto tr_yy_xx = pbuffer.data(idx_op_dd + 18);

    auto tr_yy_xy = pbuffer.data(idx_op_dd + 19);

    auto tr_yy_xz = pbuffer.data(idx_op_dd + 20);

    auto tr_yy_yy = pbuffer.data(idx_op_dd + 21);

    auto tr_yy_yz = pbuffer.data(idx_op_dd + 22);

    auto tr_yy_zz = pbuffer.data(idx_op_dd + 23);

    auto tr_yz_xx = pbuffer.data(idx_op_dd + 24);

    auto tr_yz_xy = pbuffer.data(idx_op_dd + 25);

    auto tr_yz_xz = pbuffer.data(idx_op_dd + 26);

    auto tr_yz_yy = pbuffer.data(idx_op_dd + 27);

    auto tr_yz_yz = pbuffer.data(idx_op_dd + 28);

    auto tr_yz_zz = pbuffer.data(idx_op_dd + 29);

    auto tr_zz_xx = pbuffer.data(idx_op_dd + 30);

    auto tr_zz_xy = pbuffer.data(idx_op_dd + 31);

    auto tr_zz_xz = pbuffer.data(idx_op_dd + 32);

    auto tr_zz_yy = pbuffer.data(idx_op_dd + 33);

    auto tr_zz_yz = pbuffer.data(idx_op_dd + 34);

    auto tr_zz_zz = pbuffer.data(idx_op_dd + 35);

    // Set up components of auxiliary buffer : DG

    auto tr_xx_xxxx = pbuffer.data(idx_op_dg);

    auto tr_xx_xxxy = pbuffer.data(idx_op_dg + 1);

    auto tr_xx_xxxz = pbuffer.data(idx_op_dg + 2);

    auto tr_xx_xxyy = pbuffer.data(idx_op_dg + 3);

    auto tr_xx_xxyz = pbuffer.data(idx_op_dg + 4);

    auto tr_xx_xxzz = pbuffer.data(idx_op_dg + 5);

    auto tr_xx_xyyy = pbuffer.data(idx_op_dg + 6);

    auto tr_xx_xyyz = pbuffer.data(idx_op_dg + 7);

    auto tr_xx_xyzz = pbuffer.data(idx_op_dg + 8);

    auto tr_xx_xzzz = pbuffer.data(idx_op_dg + 9);

    auto tr_xx_yyyy = pbuffer.data(idx_op_dg + 10);

    auto tr_xx_yyyz = pbuffer.data(idx_op_dg + 11);

    auto tr_xx_yyzz = pbuffer.data(idx_op_dg + 12);

    auto tr_xx_yzzz = pbuffer.data(idx_op_dg + 13);

    auto tr_xx_zzzz = pbuffer.data(idx_op_dg + 14);

    auto tr_xy_xxxx = pbuffer.data(idx_op_dg + 15);

    auto tr_xy_xxxy = pbuffer.data(idx_op_dg + 16);

    auto tr_xy_xxxz = pbuffer.data(idx_op_dg + 17);

    auto tr_xy_xxyy = pbuffer.data(idx_op_dg + 18);

    auto tr_xy_xxyz = pbuffer.data(idx_op_dg + 19);

    auto tr_xy_xxzz = pbuffer.data(idx_op_dg + 20);

    auto tr_xy_xyyy = pbuffer.data(idx_op_dg + 21);

    auto tr_xy_xyyz = pbuffer.data(idx_op_dg + 22);

    auto tr_xy_xyzz = pbuffer.data(idx_op_dg + 23);

    auto tr_xy_xzzz = pbuffer.data(idx_op_dg + 24);

    auto tr_xy_yyyy = pbuffer.data(idx_op_dg + 25);

    auto tr_xy_yyyz = pbuffer.data(idx_op_dg + 26);

    auto tr_xy_yyzz = pbuffer.data(idx_op_dg + 27);

    auto tr_xy_yzzz = pbuffer.data(idx_op_dg + 28);

    auto tr_xy_zzzz = pbuffer.data(idx_op_dg + 29);

    auto tr_xz_xxxx = pbuffer.data(idx_op_dg + 30);

    auto tr_xz_xxxy = pbuffer.data(idx_op_dg + 31);

    auto tr_xz_xxxz = pbuffer.data(idx_op_dg + 32);

    auto tr_xz_xxyy = pbuffer.data(idx_op_dg + 33);

    auto tr_xz_xxyz = pbuffer.data(idx_op_dg + 34);

    auto tr_xz_xxzz = pbuffer.data(idx_op_dg + 35);

    auto tr_xz_xyyy = pbuffer.data(idx_op_dg + 36);

    auto tr_xz_xyyz = pbuffer.data(idx_op_dg + 37);

    auto tr_xz_xyzz = pbuffer.data(idx_op_dg + 38);

    auto tr_xz_xzzz = pbuffer.data(idx_op_dg + 39);

    auto tr_xz_yyyy = pbuffer.data(idx_op_dg + 40);

    auto tr_xz_yyyz = pbuffer.data(idx_op_dg + 41);

    auto tr_xz_yyzz = pbuffer.data(idx_op_dg + 42);

    auto tr_xz_yzzz = pbuffer.data(idx_op_dg + 43);

    auto tr_xz_zzzz = pbuffer.data(idx_op_dg + 44);

    auto tr_yy_xxxx = pbuffer.data(idx_op_dg + 45);

    auto tr_yy_xxxy = pbuffer.data(idx_op_dg + 46);

    auto tr_yy_xxxz = pbuffer.data(idx_op_dg + 47);

    auto tr_yy_xxyy = pbuffer.data(idx_op_dg + 48);

    auto tr_yy_xxyz = pbuffer.data(idx_op_dg + 49);

    auto tr_yy_xxzz = pbuffer.data(idx_op_dg + 50);

    auto tr_yy_xyyy = pbuffer.data(idx_op_dg + 51);

    auto tr_yy_xyyz = pbuffer.data(idx_op_dg + 52);

    auto tr_yy_xyzz = pbuffer.data(idx_op_dg + 53);

    auto tr_yy_xzzz = pbuffer.data(idx_op_dg + 54);

    auto tr_yy_yyyy = pbuffer.data(idx_op_dg + 55);

    auto tr_yy_yyyz = pbuffer.data(idx_op_dg + 56);

    auto tr_yy_yyzz = pbuffer.data(idx_op_dg + 57);

    auto tr_yy_yzzz = pbuffer.data(idx_op_dg + 58);

    auto tr_yy_zzzz = pbuffer.data(idx_op_dg + 59);

    auto tr_yz_xxxx = pbuffer.data(idx_op_dg + 60);

    auto tr_yz_xxxy = pbuffer.data(idx_op_dg + 61);

    auto tr_yz_xxxz = pbuffer.data(idx_op_dg + 62);

    auto tr_yz_xxyy = pbuffer.data(idx_op_dg + 63);

    auto tr_yz_xxyz = pbuffer.data(idx_op_dg + 64);

    auto tr_yz_xxzz = pbuffer.data(idx_op_dg + 65);

    auto tr_yz_xyyy = pbuffer.data(idx_op_dg + 66);

    auto tr_yz_xyyz = pbuffer.data(idx_op_dg + 67);

    auto tr_yz_xyzz = pbuffer.data(idx_op_dg + 68);

    auto tr_yz_xzzz = pbuffer.data(idx_op_dg + 69);

    auto tr_yz_yyyy = pbuffer.data(idx_op_dg + 70);

    auto tr_yz_yyyz = pbuffer.data(idx_op_dg + 71);

    auto tr_yz_yyzz = pbuffer.data(idx_op_dg + 72);

    auto tr_yz_yzzz = pbuffer.data(idx_op_dg + 73);

    auto tr_yz_zzzz = pbuffer.data(idx_op_dg + 74);

    auto tr_zz_xxxx = pbuffer.data(idx_op_dg + 75);

    auto tr_zz_xxxy = pbuffer.data(idx_op_dg + 76);

    auto tr_zz_xxxz = pbuffer.data(idx_op_dg + 77);

    auto tr_zz_xxyy = pbuffer.data(idx_op_dg + 78);

    auto tr_zz_xxyz = pbuffer.data(idx_op_dg + 79);

    auto tr_zz_xxzz = pbuffer.data(idx_op_dg + 80);

    auto tr_zz_xyyy = pbuffer.data(idx_op_dg + 81);

    auto tr_zz_xyyz = pbuffer.data(idx_op_dg + 82);

    auto tr_zz_xyzz = pbuffer.data(idx_op_dg + 83);

    auto tr_zz_xzzz = pbuffer.data(idx_op_dg + 84);

    auto tr_zz_yyyy = pbuffer.data(idx_op_dg + 85);

    auto tr_zz_yyyz = pbuffer.data(idx_op_dg + 86);

    auto tr_zz_yyzz = pbuffer.data(idx_op_dg + 87);

    auto tr_zz_yzzz = pbuffer.data(idx_op_dg + 88);

    auto tr_zz_zzzz = pbuffer.data(idx_op_dg + 89);

    // Set up components of auxiliary buffer : FF

    auto tr_xxx_xxx = pbuffer.data(idx_op_ff);

    auto tr_xxx_xxy = pbuffer.data(idx_op_ff + 1);

    auto tr_xxx_xxz = pbuffer.data(idx_op_ff + 2);

    auto tr_xxx_xyy = pbuffer.data(idx_op_ff + 3);

    auto tr_xxx_xyz = pbuffer.data(idx_op_ff + 4);

    auto tr_xxx_xzz = pbuffer.data(idx_op_ff + 5);

    auto tr_xxx_yyy = pbuffer.data(idx_op_ff + 6);

    auto tr_xxx_yyz = pbuffer.data(idx_op_ff + 7);

    auto tr_xxx_yzz = pbuffer.data(idx_op_ff + 8);

    auto tr_xxx_zzz = pbuffer.data(idx_op_ff + 9);

    auto tr_xxy_xxx = pbuffer.data(idx_op_ff + 10);

    auto tr_xxy_xxy = pbuffer.data(idx_op_ff + 11);

    auto tr_xxy_xxz = pbuffer.data(idx_op_ff + 12);

    auto tr_xxy_xyy = pbuffer.data(idx_op_ff + 13);

    auto tr_xxy_xyz = pbuffer.data(idx_op_ff + 14);

    auto tr_xxy_xzz = pbuffer.data(idx_op_ff + 15);

    auto tr_xxy_yyy = pbuffer.data(idx_op_ff + 16);

    auto tr_xxy_yyz = pbuffer.data(idx_op_ff + 17);

    auto tr_xxy_yzz = pbuffer.data(idx_op_ff + 18);

    auto tr_xxy_zzz = pbuffer.data(idx_op_ff + 19);

    auto tr_xxz_xxx = pbuffer.data(idx_op_ff + 20);

    auto tr_xxz_xxy = pbuffer.data(idx_op_ff + 21);

    auto tr_xxz_xxz = pbuffer.data(idx_op_ff + 22);

    auto tr_xxz_xyy = pbuffer.data(idx_op_ff + 23);

    auto tr_xxz_xyz = pbuffer.data(idx_op_ff + 24);

    auto tr_xxz_xzz = pbuffer.data(idx_op_ff + 25);

    auto tr_xxz_yyy = pbuffer.data(idx_op_ff + 26);

    auto tr_xxz_yyz = pbuffer.data(idx_op_ff + 27);

    auto tr_xxz_yzz = pbuffer.data(idx_op_ff + 28);

    auto tr_xxz_zzz = pbuffer.data(idx_op_ff + 29);

    auto tr_xyy_xxx = pbuffer.data(idx_op_ff + 30);

    auto tr_xyy_xxy = pbuffer.data(idx_op_ff + 31);

    auto tr_xyy_xxz = pbuffer.data(idx_op_ff + 32);

    auto tr_xyy_xyy = pbuffer.data(idx_op_ff + 33);

    auto tr_xyy_xyz = pbuffer.data(idx_op_ff + 34);

    auto tr_xyy_xzz = pbuffer.data(idx_op_ff + 35);

    auto tr_xyy_yyy = pbuffer.data(idx_op_ff + 36);

    auto tr_xyy_yyz = pbuffer.data(idx_op_ff + 37);

    auto tr_xyy_yzz = pbuffer.data(idx_op_ff + 38);

    auto tr_xyy_zzz = pbuffer.data(idx_op_ff + 39);

    auto tr_xyz_xxx = pbuffer.data(idx_op_ff + 40);

    auto tr_xyz_xxy = pbuffer.data(idx_op_ff + 41);

    auto tr_xyz_xxz = pbuffer.data(idx_op_ff + 42);

    auto tr_xyz_xyy = pbuffer.data(idx_op_ff + 43);

    auto tr_xyz_xyz = pbuffer.data(idx_op_ff + 44);

    auto tr_xyz_xzz = pbuffer.data(idx_op_ff + 45);

    auto tr_xyz_yyy = pbuffer.data(idx_op_ff + 46);

    auto tr_xyz_yyz = pbuffer.data(idx_op_ff + 47);

    auto tr_xyz_yzz = pbuffer.data(idx_op_ff + 48);

    auto tr_xyz_zzz = pbuffer.data(idx_op_ff + 49);

    auto tr_xzz_xxx = pbuffer.data(idx_op_ff + 50);

    auto tr_xzz_xxy = pbuffer.data(idx_op_ff + 51);

    auto tr_xzz_xxz = pbuffer.data(idx_op_ff + 52);

    auto tr_xzz_xyy = pbuffer.data(idx_op_ff + 53);

    auto tr_xzz_xyz = pbuffer.data(idx_op_ff + 54);

    auto tr_xzz_xzz = pbuffer.data(idx_op_ff + 55);

    auto tr_xzz_yyy = pbuffer.data(idx_op_ff + 56);

    auto tr_xzz_yyz = pbuffer.data(idx_op_ff + 57);

    auto tr_xzz_yzz = pbuffer.data(idx_op_ff + 58);

    auto tr_xzz_zzz = pbuffer.data(idx_op_ff + 59);

    auto tr_yyy_xxx = pbuffer.data(idx_op_ff + 60);

    auto tr_yyy_xxy = pbuffer.data(idx_op_ff + 61);

    auto tr_yyy_xxz = pbuffer.data(idx_op_ff + 62);

    auto tr_yyy_xyy = pbuffer.data(idx_op_ff + 63);

    auto tr_yyy_xyz = pbuffer.data(idx_op_ff + 64);

    auto tr_yyy_xzz = pbuffer.data(idx_op_ff + 65);

    auto tr_yyy_yyy = pbuffer.data(idx_op_ff + 66);

    auto tr_yyy_yyz = pbuffer.data(idx_op_ff + 67);

    auto tr_yyy_yzz = pbuffer.data(idx_op_ff + 68);

    auto tr_yyy_zzz = pbuffer.data(idx_op_ff + 69);

    auto tr_yyz_xxx = pbuffer.data(idx_op_ff + 70);

    auto tr_yyz_xxy = pbuffer.data(idx_op_ff + 71);

    auto tr_yyz_xxz = pbuffer.data(idx_op_ff + 72);

    auto tr_yyz_xyy = pbuffer.data(idx_op_ff + 73);

    auto tr_yyz_xyz = pbuffer.data(idx_op_ff + 74);

    auto tr_yyz_xzz = pbuffer.data(idx_op_ff + 75);

    auto tr_yyz_yyy = pbuffer.data(idx_op_ff + 76);

    auto tr_yyz_yyz = pbuffer.data(idx_op_ff + 77);

    auto tr_yyz_yzz = pbuffer.data(idx_op_ff + 78);

    auto tr_yyz_zzz = pbuffer.data(idx_op_ff + 79);

    auto tr_yzz_xxx = pbuffer.data(idx_op_ff + 80);

    auto tr_yzz_xxy = pbuffer.data(idx_op_ff + 81);

    auto tr_yzz_xxz = pbuffer.data(idx_op_ff + 82);

    auto tr_yzz_xyy = pbuffer.data(idx_op_ff + 83);

    auto tr_yzz_xyz = pbuffer.data(idx_op_ff + 84);

    auto tr_yzz_xzz = pbuffer.data(idx_op_ff + 85);

    auto tr_yzz_yyy = pbuffer.data(idx_op_ff + 86);

    auto tr_yzz_yyz = pbuffer.data(idx_op_ff + 87);

    auto tr_yzz_yzz = pbuffer.data(idx_op_ff + 88);

    auto tr_yzz_zzz = pbuffer.data(idx_op_ff + 89);

    auto tr_zzz_xxx = pbuffer.data(idx_op_ff + 90);

    auto tr_zzz_xxy = pbuffer.data(idx_op_ff + 91);

    auto tr_zzz_xxz = pbuffer.data(idx_op_ff + 92);

    auto tr_zzz_xyy = pbuffer.data(idx_op_ff + 93);

    auto tr_zzz_xyz = pbuffer.data(idx_op_ff + 94);

    auto tr_zzz_xzz = pbuffer.data(idx_op_ff + 95);

    auto tr_zzz_yyy = pbuffer.data(idx_op_ff + 96);

    auto tr_zzz_yyz = pbuffer.data(idx_op_ff + 97);

    auto tr_zzz_yzz = pbuffer.data(idx_op_ff + 98);

    auto tr_zzz_zzz = pbuffer.data(idx_op_ff + 99);

    // Set up 0-10 components of targeted buffer : PF

    auto tr_x_0_x_x_xxx = pbuffer.data(idx_op_geom_110_pf);

    auto tr_x_0_x_x_xxy = pbuffer.data(idx_op_geom_110_pf + 1);

    auto tr_x_0_x_x_xxz = pbuffer.data(idx_op_geom_110_pf + 2);

    auto tr_x_0_x_x_xyy = pbuffer.data(idx_op_geom_110_pf + 3);

    auto tr_x_0_x_x_xyz = pbuffer.data(idx_op_geom_110_pf + 4);

    auto tr_x_0_x_x_xzz = pbuffer.data(idx_op_geom_110_pf + 5);

    auto tr_x_0_x_x_yyy = pbuffer.data(idx_op_geom_110_pf + 6);

    auto tr_x_0_x_x_yyz = pbuffer.data(idx_op_geom_110_pf + 7);

    auto tr_x_0_x_x_yzz = pbuffer.data(idx_op_geom_110_pf + 8);

    auto tr_x_0_x_x_zzz = pbuffer.data(idx_op_geom_110_pf + 9);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_0_x_x_xxx, tr_x_0_x_x_xxy, tr_x_0_x_x_xxz, tr_x_0_x_x_xyy, tr_x_0_x_x_xyz, tr_x_0_x_x_xzz, tr_x_0_x_x_yyy, tr_x_0_x_x_yyz, tr_x_0_x_x_yzz, tr_x_0_x_x_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xx, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_x_xxx[i] = 3.0 * tr_0_xx[i] - 2.0 * tr_0_xxxx[i] * tke_0 - 6.0 * tr_x_xxx[i] * tbe_0 - 6.0 * tr_xx_xx[i] * tbe_0 + 4.0 * tr_xx_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xxy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xxxy[i] * tke_0 - 6.0 * tr_x_xxy[i] * tbe_0 - 4.0 * tr_xx_xy[i] * tbe_0 + 4.0 * tr_xx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xxz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xxxz[i] * tke_0 - 6.0 * tr_x_xxz[i] * tbe_0 - 4.0 * tr_xx_xz[i] * tbe_0 + 4.0 * tr_xx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xyy[i] = tr_0_yy[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 6.0 * tr_x_xyy[i] * tbe_0 - 2.0 * tr_xx_yy[i] * tbe_0 + 4.0 * tr_xx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xyz[i] = tr_0_yz[i] - 2.0 * tr_0_xxyz[i] * tke_0 - 6.0 * tr_x_xyz[i] * tbe_0 - 2.0 * tr_xx_yz[i] * tbe_0 + 4.0 * tr_xx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xzz[i] = tr_0_zz[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 6.0 * tr_x_xzz[i] * tbe_0 - 2.0 * tr_xx_zz[i] * tbe_0 + 4.0 * tr_xx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_yyy[i] = -2.0 * tr_0_xyyy[i] * tke_0 - 6.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_xx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_yyz[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 6.0 * tr_x_yyz[i] * tbe_0 + 4.0 * tr_xx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_yzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 - 6.0 * tr_x_yzz[i] * tbe_0 + 4.0 * tr_xx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_zzz[i] = -2.0 * tr_0_xzzz[i] * tke_0 - 6.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_xx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 10-20 components of targeted buffer : PF

    auto tr_x_0_x_y_xxx = pbuffer.data(idx_op_geom_110_pf + 10);

    auto tr_x_0_x_y_xxy = pbuffer.data(idx_op_geom_110_pf + 11);

    auto tr_x_0_x_y_xxz = pbuffer.data(idx_op_geom_110_pf + 12);

    auto tr_x_0_x_y_xyy = pbuffer.data(idx_op_geom_110_pf + 13);

    auto tr_x_0_x_y_xyz = pbuffer.data(idx_op_geom_110_pf + 14);

    auto tr_x_0_x_y_xzz = pbuffer.data(idx_op_geom_110_pf + 15);

    auto tr_x_0_x_y_yyy = pbuffer.data(idx_op_geom_110_pf + 16);

    auto tr_x_0_x_y_yyz = pbuffer.data(idx_op_geom_110_pf + 17);

    auto tr_x_0_x_y_yzz = pbuffer.data(idx_op_geom_110_pf + 18);

    auto tr_x_0_x_y_zzz = pbuffer.data(idx_op_geom_110_pf + 19);

    #pragma omp simd aligned(tr_x_0_x_y_xxx, tr_x_0_x_y_xxy, tr_x_0_x_y_xxz, tr_x_0_x_y_xyy, tr_x_0_x_y_xyz, tr_x_0_x_y_xzz, tr_x_0_x_y_yyy, tr_x_0_x_y_yyz, tr_x_0_x_y_yzz, tr_x_0_x_y_zzz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_y_xxx[i] = -2.0 * tr_y_xxx[i] * tbe_0 - 6.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xxy[i] = -2.0 * tr_y_xxy[i] * tbe_0 - 4.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xxz[i] = -2.0 * tr_y_xxz[i] * tbe_0 - 4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xyy[i] = -2.0 * tr_y_xyy[i] * tbe_0 - 2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xyz[i] = -2.0 * tr_y_xyz[i] * tbe_0 - 2.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xzz[i] = -2.0 * tr_y_xzz[i] * tbe_0 - 2.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_yyy[i] = -2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_yyz[i] = -2.0 * tr_y_yyz[i] * tbe_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_yzz[i] = -2.0 * tr_y_yzz[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_zzz[i] = -2.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 20-30 components of targeted buffer : PF

    auto tr_x_0_x_z_xxx = pbuffer.data(idx_op_geom_110_pf + 20);

    auto tr_x_0_x_z_xxy = pbuffer.data(idx_op_geom_110_pf + 21);

    auto tr_x_0_x_z_xxz = pbuffer.data(idx_op_geom_110_pf + 22);

    auto tr_x_0_x_z_xyy = pbuffer.data(idx_op_geom_110_pf + 23);

    auto tr_x_0_x_z_xyz = pbuffer.data(idx_op_geom_110_pf + 24);

    auto tr_x_0_x_z_xzz = pbuffer.data(idx_op_geom_110_pf + 25);

    auto tr_x_0_x_z_yyy = pbuffer.data(idx_op_geom_110_pf + 26);

    auto tr_x_0_x_z_yyz = pbuffer.data(idx_op_geom_110_pf + 27);

    auto tr_x_0_x_z_yzz = pbuffer.data(idx_op_geom_110_pf + 28);

    auto tr_x_0_x_z_zzz = pbuffer.data(idx_op_geom_110_pf + 29);

    #pragma omp simd aligned(tr_x_0_x_z_xxx, tr_x_0_x_z_xxy, tr_x_0_x_z_xxz, tr_x_0_x_z_xyy, tr_x_0_x_z_xyz, tr_x_0_x_z_xzz, tr_x_0_x_z_yyy, tr_x_0_x_z_yyz, tr_x_0_x_z_yzz, tr_x_0_x_z_zzz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_z_xxx[i] = -2.0 * tr_z_xxx[i] * tbe_0 - 6.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xxy[i] = -2.0 * tr_z_xxy[i] * tbe_0 - 4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xxz[i] = -2.0 * tr_z_xxz[i] * tbe_0 - 4.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xyy[i] = -2.0 * tr_z_xyy[i] * tbe_0 - 2.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xyz[i] = -2.0 * tr_z_xyz[i] * tbe_0 - 2.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xzz[i] = -2.0 * tr_z_xzz[i] * tbe_0 - 2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_yyy[i] = -2.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_yyz[i] = -2.0 * tr_z_yyz[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_yzz[i] = -2.0 * tr_z_yzz[i] * tbe_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_zzz[i] = -2.0 * tr_z_zzz[i] * tbe_0 + 4.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-40 components of targeted buffer : PF

    auto tr_x_0_y_x_xxx = pbuffer.data(idx_op_geom_110_pf + 30);

    auto tr_x_0_y_x_xxy = pbuffer.data(idx_op_geom_110_pf + 31);

    auto tr_x_0_y_x_xxz = pbuffer.data(idx_op_geom_110_pf + 32);

    auto tr_x_0_y_x_xyy = pbuffer.data(idx_op_geom_110_pf + 33);

    auto tr_x_0_y_x_xyz = pbuffer.data(idx_op_geom_110_pf + 34);

    auto tr_x_0_y_x_xzz = pbuffer.data(idx_op_geom_110_pf + 35);

    auto tr_x_0_y_x_yyy = pbuffer.data(idx_op_geom_110_pf + 36);

    auto tr_x_0_y_x_yyz = pbuffer.data(idx_op_geom_110_pf + 37);

    auto tr_x_0_y_x_yzz = pbuffer.data(idx_op_geom_110_pf + 38);

    auto tr_x_0_y_x_zzz = pbuffer.data(idx_op_geom_110_pf + 39);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxy, tr_0_xxyy, tr_0_xxyz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_x_0_y_x_xxx, tr_x_0_y_x_xxy, tr_x_0_y_x_xxz, tr_x_0_y_x_xyy, tr_x_0_y_x_xyz, tr_x_0_y_x_xzz, tr_x_0_y_x_yyy, tr_x_0_y_x_yyz, tr_x_0_y_x_yzz, tr_x_0_y_x_zzz, tr_xx_xx, tr_xx_xxxy, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_yy, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_x_xxx[i] = -2.0 * tr_0_xxxy[i] * tke_0 - 2.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_xx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xxy[i] = tr_0_xx[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_y_xxy[i] * tbe_0 - 2.0 * tr_xx_xx[i] * tbe_0 + 4.0 * tr_xx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xxz[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_y_xxz[i] * tbe_0 + 4.0 * tr_xx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xyy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xyyy[i] * tke_0 - 2.0 * tr_y_xyy[i] * tbe_0 - 4.0 * tr_xx_xy[i] * tbe_0 + 4.0 * tr_xx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xyz[i] = tr_0_xz[i] - 2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_y_xyz[i] * tbe_0 - 2.0 * tr_xx_xz[i] * tbe_0 + 4.0 * tr_xx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_y_xzz[i] * tbe_0 + 4.0 * tr_xx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_yyy[i] = 3.0 * tr_0_yy[i] - 2.0 * tr_0_yyyy[i] * tke_0 - 2.0 * tr_y_yyy[i] * tbe_0 - 6.0 * tr_xx_yy[i] * tbe_0 + 4.0 * tr_xx_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_yyz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yyyz[i] * tke_0 - 2.0 * tr_y_yyz[i] * tbe_0 - 4.0 * tr_xx_yz[i] * tbe_0 + 4.0 * tr_xx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_yzz[i] = tr_0_zz[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_y_yzz[i] * tbe_0 - 2.0 * tr_xx_zz[i] * tbe_0 + 4.0 * tr_xx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_zzz[i] = -2.0 * tr_0_yzzz[i] * tke_0 - 2.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_xx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 40-50 components of targeted buffer : PF

    auto tr_x_0_y_y_xxx = pbuffer.data(idx_op_geom_110_pf + 40);

    auto tr_x_0_y_y_xxy = pbuffer.data(idx_op_geom_110_pf + 41);

    auto tr_x_0_y_y_xxz = pbuffer.data(idx_op_geom_110_pf + 42);

    auto tr_x_0_y_y_xyy = pbuffer.data(idx_op_geom_110_pf + 43);

    auto tr_x_0_y_y_xyz = pbuffer.data(idx_op_geom_110_pf + 44);

    auto tr_x_0_y_y_xzz = pbuffer.data(idx_op_geom_110_pf + 45);

    auto tr_x_0_y_y_yyy = pbuffer.data(idx_op_geom_110_pf + 46);

    auto tr_x_0_y_y_yyz = pbuffer.data(idx_op_geom_110_pf + 47);

    auto tr_x_0_y_y_yzz = pbuffer.data(idx_op_geom_110_pf + 48);

    auto tr_x_0_y_y_zzz = pbuffer.data(idx_op_geom_110_pf + 49);

    #pragma omp simd aligned(tr_x_0_y_y_xxx, tr_x_0_y_y_xxy, tr_x_0_y_y_xxz, tr_x_0_y_y_xyy, tr_x_0_y_y_xyz, tr_x_0_y_y_xzz, tr_x_0_y_y_yyy, tr_x_0_y_y_yyz, tr_x_0_y_y_yzz, tr_x_0_y_y_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_y_xxx[i] = -2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xxy[i] = -2.0 * tr_x_xxy[i] * tbe_0 - 2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xxz[i] = -2.0 * tr_x_xxz[i] * tbe_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xyy[i] = -2.0 * tr_x_xyy[i] * tbe_0 - 4.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xyz[i] = -2.0 * tr_x_xyz[i] * tbe_0 - 2.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xzz[i] = -2.0 * tr_x_xzz[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_yyy[i] = -2.0 * tr_x_yyy[i] * tbe_0 - 6.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_yyz[i] = -2.0 * tr_x_yyz[i] * tbe_0 - 4.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_yzz[i] = -2.0 * tr_x_yzz[i] * tbe_0 - 2.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_zzz[i] = -2.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 50-60 components of targeted buffer : PF

    auto tr_x_0_y_z_xxx = pbuffer.data(idx_op_geom_110_pf + 50);

    auto tr_x_0_y_z_xxy = pbuffer.data(idx_op_geom_110_pf + 51);

    auto tr_x_0_y_z_xxz = pbuffer.data(idx_op_geom_110_pf + 52);

    auto tr_x_0_y_z_xyy = pbuffer.data(idx_op_geom_110_pf + 53);

    auto tr_x_0_y_z_xyz = pbuffer.data(idx_op_geom_110_pf + 54);

    auto tr_x_0_y_z_xzz = pbuffer.data(idx_op_geom_110_pf + 55);

    auto tr_x_0_y_z_yyy = pbuffer.data(idx_op_geom_110_pf + 56);

    auto tr_x_0_y_z_yyz = pbuffer.data(idx_op_geom_110_pf + 57);

    auto tr_x_0_y_z_yzz = pbuffer.data(idx_op_geom_110_pf + 58);

    auto tr_x_0_y_z_zzz = pbuffer.data(idx_op_geom_110_pf + 59);

    #pragma omp simd aligned(tr_x_0_y_z_xxx, tr_x_0_y_z_xxy, tr_x_0_y_z_xxz, tr_x_0_y_z_xyy, tr_x_0_y_z_xyz, tr_x_0_y_z_xzz, tr_x_0_y_z_yyy, tr_x_0_y_z_yyz, tr_x_0_y_z_yzz, tr_x_0_y_z_zzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_z_xxx[i] = 4.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xxy[i] = -2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xxz[i] = 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xyy[i] = -4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xyz[i] = -2.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xzz[i] = 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_yyy[i] = -6.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_yyz[i] = -4.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_yzz[i] = -2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_zzz[i] = 4.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-70 components of targeted buffer : PF

    auto tr_x_0_z_x_xxx = pbuffer.data(idx_op_geom_110_pf + 60);

    auto tr_x_0_z_x_xxy = pbuffer.data(idx_op_geom_110_pf + 61);

    auto tr_x_0_z_x_xxz = pbuffer.data(idx_op_geom_110_pf + 62);

    auto tr_x_0_z_x_xyy = pbuffer.data(idx_op_geom_110_pf + 63);

    auto tr_x_0_z_x_xyz = pbuffer.data(idx_op_geom_110_pf + 64);

    auto tr_x_0_z_x_xzz = pbuffer.data(idx_op_geom_110_pf + 65);

    auto tr_x_0_z_x_yyy = pbuffer.data(idx_op_geom_110_pf + 66);

    auto tr_x_0_z_x_yyz = pbuffer.data(idx_op_geom_110_pf + 67);

    auto tr_x_0_z_x_yzz = pbuffer.data(idx_op_geom_110_pf + 68);

    auto tr_x_0_z_x_zzz = pbuffer.data(idx_op_geom_110_pf + 69);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxz, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_x_0_z_x_xxx, tr_x_0_z_x_xxy, tr_x_0_z_x_xxz, tr_x_0_z_x_xyy, tr_x_0_z_x_xyz, tr_x_0_z_x_xzz, tr_x_0_z_x_yyy, tr_x_0_z_x_yyz, tr_x_0_z_x_yzz, tr_x_0_z_x_zzz, tr_xx_xx, tr_xx_xxxz, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xx_zzzz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_x_xxx[i] = -2.0 * tr_0_xxxz[i] * tke_0 - 2.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_xx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xxy[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_z_xxy[i] * tbe_0 + 4.0 * tr_xx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xxz[i] = tr_0_xx[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_z_xxz[i] * tbe_0 - 2.0 * tr_xx_xx[i] * tbe_0 + 4.0 * tr_xx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xyy[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_z_xyy[i] * tbe_0 + 4.0 * tr_xx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xyz[i] = tr_0_xy[i] - 2.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_z_xyz[i] * tbe_0 - 2.0 * tr_xx_xy[i] * tbe_0 + 4.0 * tr_xx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xzz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xzzz[i] * tke_0 - 2.0 * tr_z_xzz[i] * tbe_0 - 4.0 * tr_xx_xz[i] * tbe_0 + 4.0 * tr_xx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_yyy[i] = -2.0 * tr_0_yyyz[i] * tke_0 - 2.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_xx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_yyz[i] = tr_0_yy[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_z_yyz[i] * tbe_0 - 2.0 * tr_xx_yy[i] * tbe_0 + 4.0 * tr_xx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_yzz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yzzz[i] * tke_0 - 2.0 * tr_z_yzz[i] * tbe_0 - 4.0 * tr_xx_yz[i] * tbe_0 + 4.0 * tr_xx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_zzz[i] = 3.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tke_0 - 2.0 * tr_z_zzz[i] * tbe_0 - 6.0 * tr_xx_zz[i] * tbe_0 + 4.0 * tr_xx_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 70-80 components of targeted buffer : PF

    auto tr_x_0_z_y_xxx = pbuffer.data(idx_op_geom_110_pf + 70);

    auto tr_x_0_z_y_xxy = pbuffer.data(idx_op_geom_110_pf + 71);

    auto tr_x_0_z_y_xxz = pbuffer.data(idx_op_geom_110_pf + 72);

    auto tr_x_0_z_y_xyy = pbuffer.data(idx_op_geom_110_pf + 73);

    auto tr_x_0_z_y_xyz = pbuffer.data(idx_op_geom_110_pf + 74);

    auto tr_x_0_z_y_xzz = pbuffer.data(idx_op_geom_110_pf + 75);

    auto tr_x_0_z_y_yyy = pbuffer.data(idx_op_geom_110_pf + 76);

    auto tr_x_0_z_y_yyz = pbuffer.data(idx_op_geom_110_pf + 77);

    auto tr_x_0_z_y_yzz = pbuffer.data(idx_op_geom_110_pf + 78);

    auto tr_x_0_z_y_zzz = pbuffer.data(idx_op_geom_110_pf + 79);

    #pragma omp simd aligned(tr_x_0_z_y_xxx, tr_x_0_z_y_xxy, tr_x_0_z_y_xxz, tr_x_0_z_y_xyy, tr_x_0_z_y_xyz, tr_x_0_z_y_xzz, tr_x_0_z_y_yyy, tr_x_0_z_y_yyz, tr_x_0_z_y_yzz, tr_x_0_z_y_zzz, tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_y_xxx[i] = 4.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xxy[i] = 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xxz[i] = -2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xyy[i] = 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xyz[i] = -2.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xzz[i] = -4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_yyy[i] = 4.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_yyz[i] = -2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_yzz[i] = -4.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_zzz[i] = -6.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 80-90 components of targeted buffer : PF

    auto tr_x_0_z_z_xxx = pbuffer.data(idx_op_geom_110_pf + 80);

    auto tr_x_0_z_z_xxy = pbuffer.data(idx_op_geom_110_pf + 81);

    auto tr_x_0_z_z_xxz = pbuffer.data(idx_op_geom_110_pf + 82);

    auto tr_x_0_z_z_xyy = pbuffer.data(idx_op_geom_110_pf + 83);

    auto tr_x_0_z_z_xyz = pbuffer.data(idx_op_geom_110_pf + 84);

    auto tr_x_0_z_z_xzz = pbuffer.data(idx_op_geom_110_pf + 85);

    auto tr_x_0_z_z_yyy = pbuffer.data(idx_op_geom_110_pf + 86);

    auto tr_x_0_z_z_yyz = pbuffer.data(idx_op_geom_110_pf + 87);

    auto tr_x_0_z_z_yzz = pbuffer.data(idx_op_geom_110_pf + 88);

    auto tr_x_0_z_z_zzz = pbuffer.data(idx_op_geom_110_pf + 89);

    #pragma omp simd aligned(tr_x_0_z_z_xxx, tr_x_0_z_z_xxy, tr_x_0_z_z_xxz, tr_x_0_z_z_xyy, tr_x_0_z_z_xyz, tr_x_0_z_z_xzz, tr_x_0_z_z_yyy, tr_x_0_z_z_yyz, tr_x_0_z_z_yzz, tr_x_0_z_z_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_z_xxx[i] = -2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xxy[i] = -2.0 * tr_x_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xxz[i] = -2.0 * tr_x_xxz[i] * tbe_0 - 2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xyy[i] = -2.0 * tr_x_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xyz[i] = -2.0 * tr_x_xyz[i] * tbe_0 - 2.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xzz[i] = -2.0 * tr_x_xzz[i] * tbe_0 - 4.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_yyy[i] = -2.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_yyz[i] = -2.0 * tr_x_yyz[i] * tbe_0 - 2.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_yzz[i] = -2.0 * tr_x_yzz[i] * tbe_0 - 4.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_zzz[i] = -2.0 * tr_x_zzz[i] * tbe_0 - 6.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-100 components of targeted buffer : PF

    auto tr_y_0_x_x_xxx = pbuffer.data(idx_op_geom_110_pf + 90);

    auto tr_y_0_x_x_xxy = pbuffer.data(idx_op_geom_110_pf + 91);

    auto tr_y_0_x_x_xxz = pbuffer.data(idx_op_geom_110_pf + 92);

    auto tr_y_0_x_x_xyy = pbuffer.data(idx_op_geom_110_pf + 93);

    auto tr_y_0_x_x_xyz = pbuffer.data(idx_op_geom_110_pf + 94);

    auto tr_y_0_x_x_xzz = pbuffer.data(idx_op_geom_110_pf + 95);

    auto tr_y_0_x_x_yyy = pbuffer.data(idx_op_geom_110_pf + 96);

    auto tr_y_0_x_x_yyz = pbuffer.data(idx_op_geom_110_pf + 97);

    auto tr_y_0_x_x_yzz = pbuffer.data(idx_op_geom_110_pf + 98);

    auto tr_y_0_x_x_zzz = pbuffer.data(idx_op_geom_110_pf + 99);

    #pragma omp simd aligned(tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_y_0_x_x_xxx, tr_y_0_x_x_xxy, tr_y_0_x_x_xxz, tr_y_0_x_x_xyy, tr_y_0_x_x_xyz, tr_y_0_x_x_xzz, tr_y_0_x_x_yyy, tr_y_0_x_x_yyz, tr_y_0_x_x_yzz, tr_y_0_x_x_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_x_xxx[i] = -2.0 * tr_y_xxx[i] * tbe_0 - 6.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xxy[i] = -2.0 * tr_y_xxy[i] * tbe_0 - 4.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xxz[i] = -2.0 * tr_y_xxz[i] * tbe_0 - 4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xyy[i] = -2.0 * tr_y_xyy[i] * tbe_0 - 2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xyz[i] = -2.0 * tr_y_xyz[i] * tbe_0 - 2.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xzz[i] = -2.0 * tr_y_xzz[i] * tbe_0 - 2.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_yyy[i] = -2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_yyz[i] = -2.0 * tr_y_yyz[i] * tbe_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_yzz[i] = -2.0 * tr_y_yzz[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_zzz[i] = -2.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 100-110 components of targeted buffer : PF

    auto tr_y_0_x_y_xxx = pbuffer.data(idx_op_geom_110_pf + 100);

    auto tr_y_0_x_y_xxy = pbuffer.data(idx_op_geom_110_pf + 101);

    auto tr_y_0_x_y_xxz = pbuffer.data(idx_op_geom_110_pf + 102);

    auto tr_y_0_x_y_xyy = pbuffer.data(idx_op_geom_110_pf + 103);

    auto tr_y_0_x_y_xyz = pbuffer.data(idx_op_geom_110_pf + 104);

    auto tr_y_0_x_y_xzz = pbuffer.data(idx_op_geom_110_pf + 105);

    auto tr_y_0_x_y_yyy = pbuffer.data(idx_op_geom_110_pf + 106);

    auto tr_y_0_x_y_yyz = pbuffer.data(idx_op_geom_110_pf + 107);

    auto tr_y_0_x_y_yzz = pbuffer.data(idx_op_geom_110_pf + 108);

    auto tr_y_0_x_y_zzz = pbuffer.data(idx_op_geom_110_pf + 109);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_y_0_x_y_xxx, tr_y_0_x_y_xxy, tr_y_0_x_y_xxz, tr_y_0_x_y_xyy, tr_y_0_x_y_xyz, tr_y_0_x_y_xzz, tr_y_0_x_y_yyy, tr_y_0_x_y_yyz, tr_y_0_x_y_yzz, tr_y_0_x_y_zzz, tr_yy_xx, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_y_xxx[i] = 3.0 * tr_0_xx[i] - 2.0 * tr_0_xxxx[i] * tke_0 - 6.0 * tr_yy_xx[i] * tbe_0 + 4.0 * tr_yy_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xxy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xxxy[i] * tke_0 - 4.0 * tr_yy_xy[i] * tbe_0 + 4.0 * tr_yy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxy[i] * tbe_0 + 4.0 * tr_xyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xxz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xxxz[i] * tke_0 - 4.0 * tr_yy_xz[i] * tbe_0 + 4.0 * tr_yy_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xyy[i] = tr_0_yy[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_yy_yy[i] * tbe_0 + 4.0 * tr_yy_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyy[i] * tbe_0 + 4.0 * tr_xyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xyz[i] = tr_0_yz[i] - 2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_yy_yz[i] * tbe_0 + 4.0 * tr_yy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyz[i] * tbe_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xzz[i] = tr_0_zz[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_yy_zz[i] * tbe_0 + 4.0 * tr_yy_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzz[i] * tbe_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_yyy[i] = -2.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_yy_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_xyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_yyz[i] = -2.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_yy_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyz[i] * tbe_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_yzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_yy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzz[i] * tbe_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_zzz[i] = -2.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_yy_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_xyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 110-120 components of targeted buffer : PF

    auto tr_y_0_x_z_xxx = pbuffer.data(idx_op_geom_110_pf + 110);

    auto tr_y_0_x_z_xxy = pbuffer.data(idx_op_geom_110_pf + 111);

    auto tr_y_0_x_z_xxz = pbuffer.data(idx_op_geom_110_pf + 112);

    auto tr_y_0_x_z_xyy = pbuffer.data(idx_op_geom_110_pf + 113);

    auto tr_y_0_x_z_xyz = pbuffer.data(idx_op_geom_110_pf + 114);

    auto tr_y_0_x_z_xzz = pbuffer.data(idx_op_geom_110_pf + 115);

    auto tr_y_0_x_z_yyy = pbuffer.data(idx_op_geom_110_pf + 116);

    auto tr_y_0_x_z_yyz = pbuffer.data(idx_op_geom_110_pf + 117);

    auto tr_y_0_x_z_yzz = pbuffer.data(idx_op_geom_110_pf + 118);

    auto tr_y_0_x_z_zzz = pbuffer.data(idx_op_geom_110_pf + 119);

    #pragma omp simd aligned(tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_y_0_x_z_xxx, tr_y_0_x_z_xxy, tr_y_0_x_z_xxz, tr_y_0_x_z_xyy, tr_y_0_x_z_xyz, tr_y_0_x_z_xzz, tr_y_0_x_z_yyy, tr_y_0_x_z_yyz, tr_y_0_x_z_yzz, tr_y_0_x_z_zzz, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_z_xxx[i] = -6.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xxy[i] = -4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xxz[i] = -4.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xyy[i] = -2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xyz[i] = -2.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xzz[i] = -2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_yyy[i] = 4.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_yyz[i] = 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_yzz[i] = 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_zzz[i] = 4.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-130 components of targeted buffer : PF

    auto tr_y_0_y_x_xxx = pbuffer.data(idx_op_geom_110_pf + 120);

    auto tr_y_0_y_x_xxy = pbuffer.data(idx_op_geom_110_pf + 121);

    auto tr_y_0_y_x_xxz = pbuffer.data(idx_op_geom_110_pf + 122);

    auto tr_y_0_y_x_xyy = pbuffer.data(idx_op_geom_110_pf + 123);

    auto tr_y_0_y_x_xyz = pbuffer.data(idx_op_geom_110_pf + 124);

    auto tr_y_0_y_x_xzz = pbuffer.data(idx_op_geom_110_pf + 125);

    auto tr_y_0_y_x_yyy = pbuffer.data(idx_op_geom_110_pf + 126);

    auto tr_y_0_y_x_yyz = pbuffer.data(idx_op_geom_110_pf + 127);

    auto tr_y_0_y_x_yzz = pbuffer.data(idx_op_geom_110_pf + 128);

    auto tr_y_0_y_x_zzz = pbuffer.data(idx_op_geom_110_pf + 129);

    #pragma omp simd aligned(tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_y_0_y_x_xxx, tr_y_0_y_x_xxy, tr_y_0_y_x_xxz, tr_y_0_y_x_xyy, tr_y_0_y_x_xyz, tr_y_0_y_x_xzz, tr_y_0_y_x_yyy, tr_y_0_y_x_yyz, tr_y_0_y_x_yzz, tr_y_0_y_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_x_xxx[i] = -2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xxy[i] = -2.0 * tr_x_xxy[i] * tbe_0 - 2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xxz[i] = -2.0 * tr_x_xxz[i] * tbe_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xyy[i] = -2.0 * tr_x_xyy[i] * tbe_0 - 4.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xyz[i] = -2.0 * tr_x_xyz[i] * tbe_0 - 2.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xzz[i] = -2.0 * tr_x_xzz[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_yyy[i] = -2.0 * tr_x_yyy[i] * tbe_0 - 6.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_yyz[i] = -2.0 * tr_x_yyz[i] * tbe_0 - 4.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_yzz[i] = -2.0 * tr_x_yzz[i] * tbe_0 - 2.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_zzz[i] = -2.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 130-140 components of targeted buffer : PF

    auto tr_y_0_y_y_xxx = pbuffer.data(idx_op_geom_110_pf + 130);

    auto tr_y_0_y_y_xxy = pbuffer.data(idx_op_geom_110_pf + 131);

    auto tr_y_0_y_y_xxz = pbuffer.data(idx_op_geom_110_pf + 132);

    auto tr_y_0_y_y_xyy = pbuffer.data(idx_op_geom_110_pf + 133);

    auto tr_y_0_y_y_xyz = pbuffer.data(idx_op_geom_110_pf + 134);

    auto tr_y_0_y_y_xzz = pbuffer.data(idx_op_geom_110_pf + 135);

    auto tr_y_0_y_y_yyy = pbuffer.data(idx_op_geom_110_pf + 136);

    auto tr_y_0_y_y_yyz = pbuffer.data(idx_op_geom_110_pf + 137);

    auto tr_y_0_y_y_yzz = pbuffer.data(idx_op_geom_110_pf + 138);

    auto tr_y_0_y_y_zzz = pbuffer.data(idx_op_geom_110_pf + 139);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxy, tr_0_xxyy, tr_0_xxyz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_y_0_y_y_xxx, tr_y_0_y_y_xxy, tr_y_0_y_y_xxz, tr_y_0_y_y_xyy, tr_y_0_y_y_xyz, tr_y_0_y_y_xzz, tr_y_0_y_y_yyy, tr_y_0_y_y_yyz, tr_y_0_y_y_yzz, tr_y_0_y_y_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yy_xx, tr_yy_xxxy, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_yy, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_y_xxx[i] = -2.0 * tr_0_xxxy[i] * tke_0 - 6.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xxy[i] = tr_0_xx[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 6.0 * tr_y_xxy[i] * tbe_0 - 2.0 * tr_yy_xx[i] * tbe_0 + 4.0 * tr_yy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xxz[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 6.0 * tr_y_xxz[i] * tbe_0 + 4.0 * tr_yy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xyy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xyyy[i] * tke_0 - 6.0 * tr_y_xyy[i] * tbe_0 - 4.0 * tr_yy_xy[i] * tbe_0 + 4.0 * tr_yy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xyz[i] = tr_0_xz[i] - 2.0 * tr_0_xyyz[i] * tke_0 - 6.0 * tr_y_xyz[i] * tbe_0 - 2.0 * tr_yy_xz[i] * tbe_0 + 4.0 * tr_yy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 - 6.0 * tr_y_xzz[i] * tbe_0 + 4.0 * tr_yy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_yyy[i] = 3.0 * tr_0_yy[i] - 2.0 * tr_0_yyyy[i] * tke_0 - 6.0 * tr_y_yyy[i] * tbe_0 - 6.0 * tr_yy_yy[i] * tbe_0 + 4.0 * tr_yy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_yyz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yyyz[i] * tke_0 - 6.0 * tr_y_yyz[i] * tbe_0 - 4.0 * tr_yy_yz[i] * tbe_0 + 4.0 * tr_yy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_yzz[i] = tr_0_zz[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 6.0 * tr_y_yzz[i] * tbe_0 - 2.0 * tr_yy_zz[i] * tbe_0 + 4.0 * tr_yy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_zzz[i] = -2.0 * tr_0_yzzz[i] * tke_0 - 6.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_yy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 140-150 components of targeted buffer : PF

    auto tr_y_0_y_z_xxx = pbuffer.data(idx_op_geom_110_pf + 140);

    auto tr_y_0_y_z_xxy = pbuffer.data(idx_op_geom_110_pf + 141);

    auto tr_y_0_y_z_xxz = pbuffer.data(idx_op_geom_110_pf + 142);

    auto tr_y_0_y_z_xyy = pbuffer.data(idx_op_geom_110_pf + 143);

    auto tr_y_0_y_z_xyz = pbuffer.data(idx_op_geom_110_pf + 144);

    auto tr_y_0_y_z_xzz = pbuffer.data(idx_op_geom_110_pf + 145);

    auto tr_y_0_y_z_yyy = pbuffer.data(idx_op_geom_110_pf + 146);

    auto tr_y_0_y_z_yyz = pbuffer.data(idx_op_geom_110_pf + 147);

    auto tr_y_0_y_z_yzz = pbuffer.data(idx_op_geom_110_pf + 148);

    auto tr_y_0_y_z_zzz = pbuffer.data(idx_op_geom_110_pf + 149);

    #pragma omp simd aligned(tr_y_0_y_z_xxx, tr_y_0_y_z_xxy, tr_y_0_y_z_xxz, tr_y_0_y_z_xyy, tr_y_0_y_z_xyz, tr_y_0_y_z_xzz, tr_y_0_y_z_yyy, tr_y_0_y_z_yyz, tr_y_0_y_z_yzz, tr_y_0_y_z_zzz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_z_xxx[i] = -2.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xxy[i] = -2.0 * tr_z_xxy[i] * tbe_0 - 2.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xxz[i] = -2.0 * tr_z_xxz[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xyy[i] = -2.0 * tr_z_xyy[i] * tbe_0 - 4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xyz[i] = -2.0 * tr_z_xyz[i] * tbe_0 - 2.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xzz[i] = -2.0 * tr_z_xzz[i] * tbe_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_yyy[i] = -2.0 * tr_z_yyy[i] * tbe_0 - 6.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_yyz[i] = -2.0 * tr_z_yyz[i] * tbe_0 - 4.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_yzz[i] = -2.0 * tr_z_yzz[i] * tbe_0 - 2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_zzz[i] = -2.0 * tr_z_zzz[i] * tbe_0 + 4.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-160 components of targeted buffer : PF

    auto tr_y_0_z_x_xxx = pbuffer.data(idx_op_geom_110_pf + 150);

    auto tr_y_0_z_x_xxy = pbuffer.data(idx_op_geom_110_pf + 151);

    auto tr_y_0_z_x_xxz = pbuffer.data(idx_op_geom_110_pf + 152);

    auto tr_y_0_z_x_xyy = pbuffer.data(idx_op_geom_110_pf + 153);

    auto tr_y_0_z_x_xyz = pbuffer.data(idx_op_geom_110_pf + 154);

    auto tr_y_0_z_x_xzz = pbuffer.data(idx_op_geom_110_pf + 155);

    auto tr_y_0_z_x_yyy = pbuffer.data(idx_op_geom_110_pf + 156);

    auto tr_y_0_z_x_yyz = pbuffer.data(idx_op_geom_110_pf + 157);

    auto tr_y_0_z_x_yzz = pbuffer.data(idx_op_geom_110_pf + 158);

    auto tr_y_0_z_x_zzz = pbuffer.data(idx_op_geom_110_pf + 159);

    #pragma omp simd aligned(tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_y_0_z_x_xxx, tr_y_0_z_x_xxy, tr_y_0_z_x_xxz, tr_y_0_z_x_xyy, tr_y_0_z_x_xyz, tr_y_0_z_x_xzz, tr_y_0_z_x_yyy, tr_y_0_z_x_yyz, tr_y_0_z_x_yzz, tr_y_0_z_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_x_xxx[i] = 4.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xxy[i] = 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xxz[i] = -2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xyy[i] = 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xyz[i] = -2.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xzz[i] = -4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_yyy[i] = 4.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_yyz[i] = -2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_yzz[i] = -4.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_zzz[i] = -6.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 160-170 components of targeted buffer : PF

    auto tr_y_0_z_y_xxx = pbuffer.data(idx_op_geom_110_pf + 160);

    auto tr_y_0_z_y_xxy = pbuffer.data(idx_op_geom_110_pf + 161);

    auto tr_y_0_z_y_xxz = pbuffer.data(idx_op_geom_110_pf + 162);

    auto tr_y_0_z_y_xyy = pbuffer.data(idx_op_geom_110_pf + 163);

    auto tr_y_0_z_y_xyz = pbuffer.data(idx_op_geom_110_pf + 164);

    auto tr_y_0_z_y_xzz = pbuffer.data(idx_op_geom_110_pf + 165);

    auto tr_y_0_z_y_yyy = pbuffer.data(idx_op_geom_110_pf + 166);

    auto tr_y_0_z_y_yyz = pbuffer.data(idx_op_geom_110_pf + 167);

    auto tr_y_0_z_y_yzz = pbuffer.data(idx_op_geom_110_pf + 168);

    auto tr_y_0_z_y_zzz = pbuffer.data(idx_op_geom_110_pf + 169);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxz, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_y_0_z_y_xxx, tr_y_0_z_y_xxy, tr_y_0_z_y_xxz, tr_y_0_z_y_xyy, tr_y_0_z_y_xyz, tr_y_0_z_y_xzz, tr_y_0_z_y_yyy, tr_y_0_z_y_yyz, tr_y_0_z_y_yzz, tr_y_0_z_y_zzz, tr_yy_xx, tr_yy_xxxz, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_y_xxx[i] = -2.0 * tr_0_xxxz[i] * tke_0 - 2.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xxy[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_z_xxy[i] * tbe_0 + 4.0 * tr_yy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xxz[i] = tr_0_xx[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_z_xxz[i] * tbe_0 - 2.0 * tr_yy_xx[i] * tbe_0 + 4.0 * tr_yy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xyy[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_z_xyy[i] * tbe_0 + 4.0 * tr_yy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xyz[i] = tr_0_xy[i] - 2.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_z_xyz[i] * tbe_0 - 2.0 * tr_yy_xy[i] * tbe_0 + 4.0 * tr_yy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xzz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xzzz[i] * tke_0 - 2.0 * tr_z_xzz[i] * tbe_0 - 4.0 * tr_yy_xz[i] * tbe_0 + 4.0 * tr_yy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_yyy[i] = -2.0 * tr_0_yyyz[i] * tke_0 - 2.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_yy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_yyz[i] = tr_0_yy[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_z_yyz[i] * tbe_0 - 2.0 * tr_yy_yy[i] * tbe_0 + 4.0 * tr_yy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_yzz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yzzz[i] * tke_0 - 2.0 * tr_z_yzz[i] * tbe_0 - 4.0 * tr_yy_yz[i] * tbe_0 + 4.0 * tr_yy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_zzz[i] = 3.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tke_0 - 2.0 * tr_z_zzz[i] * tbe_0 - 6.0 * tr_yy_zz[i] * tbe_0 + 4.0 * tr_yy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 170-180 components of targeted buffer : PF

    auto tr_y_0_z_z_xxx = pbuffer.data(idx_op_geom_110_pf + 170);

    auto tr_y_0_z_z_xxy = pbuffer.data(idx_op_geom_110_pf + 171);

    auto tr_y_0_z_z_xxz = pbuffer.data(idx_op_geom_110_pf + 172);

    auto tr_y_0_z_z_xyy = pbuffer.data(idx_op_geom_110_pf + 173);

    auto tr_y_0_z_z_xyz = pbuffer.data(idx_op_geom_110_pf + 174);

    auto tr_y_0_z_z_xzz = pbuffer.data(idx_op_geom_110_pf + 175);

    auto tr_y_0_z_z_yyy = pbuffer.data(idx_op_geom_110_pf + 176);

    auto tr_y_0_z_z_yyz = pbuffer.data(idx_op_geom_110_pf + 177);

    auto tr_y_0_z_z_yzz = pbuffer.data(idx_op_geom_110_pf + 178);

    auto tr_y_0_z_z_zzz = pbuffer.data(idx_op_geom_110_pf + 179);

    #pragma omp simd aligned(tr_y_0_z_z_xxx, tr_y_0_z_z_xxy, tr_y_0_z_z_xxz, tr_y_0_z_z_xyy, tr_y_0_z_z_xyz, tr_y_0_z_z_xzz, tr_y_0_z_z_yyy, tr_y_0_z_z_yyz, tr_y_0_z_z_yzz, tr_y_0_z_z_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_z_xxx[i] = -2.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xxy[i] = -2.0 * tr_y_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xxz[i] = -2.0 * tr_y_xxz[i] * tbe_0 - 2.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xyy[i] = -2.0 * tr_y_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xyz[i] = -2.0 * tr_y_xyz[i] * tbe_0 - 2.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xzz[i] = -2.0 * tr_y_xzz[i] * tbe_0 - 4.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_yyy[i] = -2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_yyz[i] = -2.0 * tr_y_yyz[i] * tbe_0 - 2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_yzz[i] = -2.0 * tr_y_yzz[i] * tbe_0 - 4.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_zzz[i] = -2.0 * tr_y_zzz[i] * tbe_0 - 6.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-190 components of targeted buffer : PF

    auto tr_z_0_x_x_xxx = pbuffer.data(idx_op_geom_110_pf + 180);

    auto tr_z_0_x_x_xxy = pbuffer.data(idx_op_geom_110_pf + 181);

    auto tr_z_0_x_x_xxz = pbuffer.data(idx_op_geom_110_pf + 182);

    auto tr_z_0_x_x_xyy = pbuffer.data(idx_op_geom_110_pf + 183);

    auto tr_z_0_x_x_xyz = pbuffer.data(idx_op_geom_110_pf + 184);

    auto tr_z_0_x_x_xzz = pbuffer.data(idx_op_geom_110_pf + 185);

    auto tr_z_0_x_x_yyy = pbuffer.data(idx_op_geom_110_pf + 186);

    auto tr_z_0_x_x_yyz = pbuffer.data(idx_op_geom_110_pf + 187);

    auto tr_z_0_x_x_yzz = pbuffer.data(idx_op_geom_110_pf + 188);

    auto tr_z_0_x_x_zzz = pbuffer.data(idx_op_geom_110_pf + 189);

    #pragma omp simd aligned(tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_z_0_x_x_xxx, tr_z_0_x_x_xxy, tr_z_0_x_x_xxz, tr_z_0_x_x_xyy, tr_z_0_x_x_xyz, tr_z_0_x_x_xzz, tr_z_0_x_x_yyy, tr_z_0_x_x_yyz, tr_z_0_x_x_yzz, tr_z_0_x_x_zzz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_x_xxx[i] = -2.0 * tr_z_xxx[i] * tbe_0 - 6.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xxy[i] = -2.0 * tr_z_xxy[i] * tbe_0 - 4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xxz[i] = -2.0 * tr_z_xxz[i] * tbe_0 - 4.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xyy[i] = -2.0 * tr_z_xyy[i] * tbe_0 - 2.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xyz[i] = -2.0 * tr_z_xyz[i] * tbe_0 - 2.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xzz[i] = -2.0 * tr_z_xzz[i] * tbe_0 - 2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_yyy[i] = -2.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_yyz[i] = -2.0 * tr_z_yyz[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_yzz[i] = -2.0 * tr_z_yzz[i] * tbe_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_zzz[i] = -2.0 * tr_z_zzz[i] * tbe_0 + 4.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 190-200 components of targeted buffer : PF

    auto tr_z_0_x_y_xxx = pbuffer.data(idx_op_geom_110_pf + 190);

    auto tr_z_0_x_y_xxy = pbuffer.data(idx_op_geom_110_pf + 191);

    auto tr_z_0_x_y_xxz = pbuffer.data(idx_op_geom_110_pf + 192);

    auto tr_z_0_x_y_xyy = pbuffer.data(idx_op_geom_110_pf + 193);

    auto tr_z_0_x_y_xyz = pbuffer.data(idx_op_geom_110_pf + 194);

    auto tr_z_0_x_y_xzz = pbuffer.data(idx_op_geom_110_pf + 195);

    auto tr_z_0_x_y_yyy = pbuffer.data(idx_op_geom_110_pf + 196);

    auto tr_z_0_x_y_yyz = pbuffer.data(idx_op_geom_110_pf + 197);

    auto tr_z_0_x_y_yzz = pbuffer.data(idx_op_geom_110_pf + 198);

    auto tr_z_0_x_y_zzz = pbuffer.data(idx_op_geom_110_pf + 199);

    #pragma omp simd aligned(tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_z_0_x_y_xxx, tr_z_0_x_y_xxy, tr_z_0_x_y_xxz, tr_z_0_x_y_xyy, tr_z_0_x_y_xyz, tr_z_0_x_y_xzz, tr_z_0_x_y_yyy, tr_z_0_x_y_yyz, tr_z_0_x_y_yzz, tr_z_0_x_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_y_xxx[i] = -6.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xxy[i] = -4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xxz[i] = -4.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xyy[i] = -2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xyz[i] = -2.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xzz[i] = -2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_yyy[i] = 4.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_yyz[i] = 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_yzz[i] = 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_zzz[i] = 4.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 200-210 components of targeted buffer : PF

    auto tr_z_0_x_z_xxx = pbuffer.data(idx_op_geom_110_pf + 200);

    auto tr_z_0_x_z_xxy = pbuffer.data(idx_op_geom_110_pf + 201);

    auto tr_z_0_x_z_xxz = pbuffer.data(idx_op_geom_110_pf + 202);

    auto tr_z_0_x_z_xyy = pbuffer.data(idx_op_geom_110_pf + 203);

    auto tr_z_0_x_z_xyz = pbuffer.data(idx_op_geom_110_pf + 204);

    auto tr_z_0_x_z_xzz = pbuffer.data(idx_op_geom_110_pf + 205);

    auto tr_z_0_x_z_yyy = pbuffer.data(idx_op_geom_110_pf + 206);

    auto tr_z_0_x_z_yyz = pbuffer.data(idx_op_geom_110_pf + 207);

    auto tr_z_0_x_z_yzz = pbuffer.data(idx_op_geom_110_pf + 208);

    auto tr_z_0_x_z_zzz = pbuffer.data(idx_op_geom_110_pf + 209);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_z_0_x_z_xxx, tr_z_0_x_z_xxy, tr_z_0_x_z_xxz, tr_z_0_x_z_xyy, tr_z_0_x_z_xyz, tr_z_0_x_z_xzz, tr_z_0_x_z_yyy, tr_z_0_x_z_yyz, tr_z_0_x_z_yzz, tr_z_0_x_z_zzz, tr_zz_xx, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_z_xxx[i] = 3.0 * tr_0_xx[i] - 2.0 * tr_0_xxxx[i] * tke_0 - 6.0 * tr_zz_xx[i] * tbe_0 + 4.0 * tr_zz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xxy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xxxy[i] * tke_0 - 4.0 * tr_zz_xy[i] * tbe_0 + 4.0 * tr_zz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xxz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xxxz[i] * tke_0 - 4.0 * tr_zz_xz[i] * tbe_0 + 4.0 * tr_zz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxz[i] * tbe_0 + 4.0 * tr_xzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xyy[i] = tr_0_yy[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_zz_yy[i] * tbe_0 + 4.0 * tr_zz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyy[i] * tbe_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xyz[i] = tr_0_yz[i] - 2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_zz_yz[i] * tbe_0 + 4.0 * tr_zz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyz[i] * tbe_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xzz[i] = tr_0_zz[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_zz_zz[i] * tbe_0 + 4.0 * tr_zz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzz[i] * tbe_0 + 4.0 * tr_xzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_yyy[i] = -2.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_zz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_xzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_yyz[i] = -2.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_zz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyz[i] * tbe_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_yzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_zz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzz[i] * tbe_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_zzz[i] = -2.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_zz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_xzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-220 components of targeted buffer : PF

    auto tr_z_0_y_x_xxx = pbuffer.data(idx_op_geom_110_pf + 210);

    auto tr_z_0_y_x_xxy = pbuffer.data(idx_op_geom_110_pf + 211);

    auto tr_z_0_y_x_xxz = pbuffer.data(idx_op_geom_110_pf + 212);

    auto tr_z_0_y_x_xyy = pbuffer.data(idx_op_geom_110_pf + 213);

    auto tr_z_0_y_x_xyz = pbuffer.data(idx_op_geom_110_pf + 214);

    auto tr_z_0_y_x_xzz = pbuffer.data(idx_op_geom_110_pf + 215);

    auto tr_z_0_y_x_yyy = pbuffer.data(idx_op_geom_110_pf + 216);

    auto tr_z_0_y_x_yyz = pbuffer.data(idx_op_geom_110_pf + 217);

    auto tr_z_0_y_x_yzz = pbuffer.data(idx_op_geom_110_pf + 218);

    auto tr_z_0_y_x_zzz = pbuffer.data(idx_op_geom_110_pf + 219);

    #pragma omp simd aligned(tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_z_0_y_x_xxx, tr_z_0_y_x_xxy, tr_z_0_y_x_xxz, tr_z_0_y_x_xyy, tr_z_0_y_x_xyz, tr_z_0_y_x_xzz, tr_z_0_y_x_yyy, tr_z_0_y_x_yyz, tr_z_0_y_x_yzz, tr_z_0_y_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_x_xxx[i] = 4.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xxy[i] = -2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xxz[i] = 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xyy[i] = -4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xyz[i] = -2.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xzz[i] = 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_yyy[i] = -6.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_yyz[i] = -4.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_yzz[i] = -2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_zzz[i] = 4.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 220-230 components of targeted buffer : PF

    auto tr_z_0_y_y_xxx = pbuffer.data(idx_op_geom_110_pf + 220);

    auto tr_z_0_y_y_xxy = pbuffer.data(idx_op_geom_110_pf + 221);

    auto tr_z_0_y_y_xxz = pbuffer.data(idx_op_geom_110_pf + 222);

    auto tr_z_0_y_y_xyy = pbuffer.data(idx_op_geom_110_pf + 223);

    auto tr_z_0_y_y_xyz = pbuffer.data(idx_op_geom_110_pf + 224);

    auto tr_z_0_y_y_xzz = pbuffer.data(idx_op_geom_110_pf + 225);

    auto tr_z_0_y_y_yyy = pbuffer.data(idx_op_geom_110_pf + 226);

    auto tr_z_0_y_y_yyz = pbuffer.data(idx_op_geom_110_pf + 227);

    auto tr_z_0_y_y_yzz = pbuffer.data(idx_op_geom_110_pf + 228);

    auto tr_z_0_y_y_zzz = pbuffer.data(idx_op_geom_110_pf + 229);

    #pragma omp simd aligned(tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_z_0_y_y_xxx, tr_z_0_y_y_xxy, tr_z_0_y_y_xxz, tr_z_0_y_y_xyy, tr_z_0_y_y_xyz, tr_z_0_y_y_xzz, tr_z_0_y_y_yyy, tr_z_0_y_y_yyz, tr_z_0_y_y_yzz, tr_z_0_y_y_zzz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_y_xxx[i] = -2.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xxy[i] = -2.0 * tr_z_xxy[i] * tbe_0 - 2.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xxz[i] = -2.0 * tr_z_xxz[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xyy[i] = -2.0 * tr_z_xyy[i] * tbe_0 - 4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xyz[i] = -2.0 * tr_z_xyz[i] * tbe_0 - 2.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xzz[i] = -2.0 * tr_z_xzz[i] * tbe_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_yyy[i] = -2.0 * tr_z_yyy[i] * tbe_0 - 6.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_yyz[i] = -2.0 * tr_z_yyz[i] * tbe_0 - 4.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_yzz[i] = -2.0 * tr_z_yzz[i] * tbe_0 - 2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_zzz[i] = -2.0 * tr_z_zzz[i] * tbe_0 + 4.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 230-240 components of targeted buffer : PF

    auto tr_z_0_y_z_xxx = pbuffer.data(idx_op_geom_110_pf + 230);

    auto tr_z_0_y_z_xxy = pbuffer.data(idx_op_geom_110_pf + 231);

    auto tr_z_0_y_z_xxz = pbuffer.data(idx_op_geom_110_pf + 232);

    auto tr_z_0_y_z_xyy = pbuffer.data(idx_op_geom_110_pf + 233);

    auto tr_z_0_y_z_xyz = pbuffer.data(idx_op_geom_110_pf + 234);

    auto tr_z_0_y_z_xzz = pbuffer.data(idx_op_geom_110_pf + 235);

    auto tr_z_0_y_z_yyy = pbuffer.data(idx_op_geom_110_pf + 236);

    auto tr_z_0_y_z_yyz = pbuffer.data(idx_op_geom_110_pf + 237);

    auto tr_z_0_y_z_yzz = pbuffer.data(idx_op_geom_110_pf + 238);

    auto tr_z_0_y_z_zzz = pbuffer.data(idx_op_geom_110_pf + 239);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxy, tr_0_xxyy, tr_0_xxyz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_z_0_y_z_xxx, tr_z_0_y_z_xxy, tr_z_0_y_z_xxz, tr_z_0_y_z_xyy, tr_z_0_y_z_xyz, tr_z_0_y_z_xzz, tr_z_0_y_z_yyy, tr_z_0_y_z_yyz, tr_z_0_y_z_yzz, tr_z_0_y_z_zzz, tr_zz_xx, tr_zz_xxxy, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_yy, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_z_xxx[i] = -2.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_zz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xxy[i] = tr_0_xx[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_zz_xx[i] * tbe_0 + 4.0 * tr_zz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xxz[i] = -2.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_zz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxz[i] * tbe_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xyy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xyyy[i] * tke_0 - 4.0 * tr_zz_xy[i] * tbe_0 + 4.0 * tr_zz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyy[i] * tbe_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xyz[i] = tr_0_xz[i] - 2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_zz_xz[i] * tbe_0 + 4.0 * tr_zz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyz[i] * tbe_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_zz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xzz[i] * tbe_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_yyy[i] = 3.0 * tr_0_yy[i] - 2.0 * tr_0_yyyy[i] * tke_0 - 6.0 * tr_zz_yy[i] * tbe_0 + 4.0 * tr_zz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_yzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_yyz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yyyz[i] * tke_0 - 4.0 * tr_zz_yz[i] * tbe_0 + 4.0 * tr_zz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyz[i] * tbe_0 + 4.0 * tr_yzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_yzz[i] = tr_0_zz[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_zz_zz[i] * tbe_0 + 4.0 * tr_zz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yzz[i] * tbe_0 + 4.0 * tr_yzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_zzz[i] = -2.0 * tr_0_yzzz[i] * tke_0 + 4.0 * tr_zz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_yzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-250 components of targeted buffer : PF

    auto tr_z_0_z_x_xxx = pbuffer.data(idx_op_geom_110_pf + 240);

    auto tr_z_0_z_x_xxy = pbuffer.data(idx_op_geom_110_pf + 241);

    auto tr_z_0_z_x_xxz = pbuffer.data(idx_op_geom_110_pf + 242);

    auto tr_z_0_z_x_xyy = pbuffer.data(idx_op_geom_110_pf + 243);

    auto tr_z_0_z_x_xyz = pbuffer.data(idx_op_geom_110_pf + 244);

    auto tr_z_0_z_x_xzz = pbuffer.data(idx_op_geom_110_pf + 245);

    auto tr_z_0_z_x_yyy = pbuffer.data(idx_op_geom_110_pf + 246);

    auto tr_z_0_z_x_yyz = pbuffer.data(idx_op_geom_110_pf + 247);

    auto tr_z_0_z_x_yzz = pbuffer.data(idx_op_geom_110_pf + 248);

    auto tr_z_0_z_x_zzz = pbuffer.data(idx_op_geom_110_pf + 249);

    #pragma omp simd aligned(tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_z_0_z_x_xxx, tr_z_0_z_x_xxy, tr_z_0_z_x_xxz, tr_z_0_z_x_xyy, tr_z_0_z_x_xyz, tr_z_0_z_x_xzz, tr_z_0_z_x_yyy, tr_z_0_z_x_yyz, tr_z_0_z_x_yzz, tr_z_0_z_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_x_xxx[i] = -2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xxy[i] = -2.0 * tr_x_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xxz[i] = -2.0 * tr_x_xxz[i] * tbe_0 - 2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xyy[i] = -2.0 * tr_x_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xyz[i] = -2.0 * tr_x_xyz[i] * tbe_0 - 2.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xzz[i] = -2.0 * tr_x_xzz[i] * tbe_0 - 4.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_yyy[i] = -2.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_yyz[i] = -2.0 * tr_x_yyz[i] * tbe_0 - 2.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_yzz[i] = -2.0 * tr_x_yzz[i] * tbe_0 - 4.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_zzz[i] = -2.0 * tr_x_zzz[i] * tbe_0 - 6.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 250-260 components of targeted buffer : PF

    auto tr_z_0_z_y_xxx = pbuffer.data(idx_op_geom_110_pf + 250);

    auto tr_z_0_z_y_xxy = pbuffer.data(idx_op_geom_110_pf + 251);

    auto tr_z_0_z_y_xxz = pbuffer.data(idx_op_geom_110_pf + 252);

    auto tr_z_0_z_y_xyy = pbuffer.data(idx_op_geom_110_pf + 253);

    auto tr_z_0_z_y_xyz = pbuffer.data(idx_op_geom_110_pf + 254);

    auto tr_z_0_z_y_xzz = pbuffer.data(idx_op_geom_110_pf + 255);

    auto tr_z_0_z_y_yyy = pbuffer.data(idx_op_geom_110_pf + 256);

    auto tr_z_0_z_y_yyz = pbuffer.data(idx_op_geom_110_pf + 257);

    auto tr_z_0_z_y_yzz = pbuffer.data(idx_op_geom_110_pf + 258);

    auto tr_z_0_z_y_zzz = pbuffer.data(idx_op_geom_110_pf + 259);

    #pragma omp simd aligned(tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_z_0_z_y_xxx, tr_z_0_z_y_xxy, tr_z_0_z_y_xxz, tr_z_0_z_y_xyy, tr_z_0_z_y_xyz, tr_z_0_z_y_xzz, tr_z_0_z_y_yyy, tr_z_0_z_y_yyz, tr_z_0_z_y_yzz, tr_z_0_z_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_y_xxx[i] = -2.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xxy[i] = -2.0 * tr_y_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xxz[i] = -2.0 * tr_y_xxz[i] * tbe_0 - 2.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xyy[i] = -2.0 * tr_y_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xyz[i] = -2.0 * tr_y_xyz[i] * tbe_0 - 2.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xzz[i] = -2.0 * tr_y_xzz[i] * tbe_0 - 4.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_yyy[i] = -2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_yyz[i] = -2.0 * tr_y_yyz[i] * tbe_0 - 2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_yzz[i] = -2.0 * tr_y_yzz[i] * tbe_0 - 4.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_zzz[i] = -2.0 * tr_y_zzz[i] * tbe_0 - 6.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 260-270 components of targeted buffer : PF

    auto tr_z_0_z_z_xxx = pbuffer.data(idx_op_geom_110_pf + 260);

    auto tr_z_0_z_z_xxy = pbuffer.data(idx_op_geom_110_pf + 261);

    auto tr_z_0_z_z_xxz = pbuffer.data(idx_op_geom_110_pf + 262);

    auto tr_z_0_z_z_xyy = pbuffer.data(idx_op_geom_110_pf + 263);

    auto tr_z_0_z_z_xyz = pbuffer.data(idx_op_geom_110_pf + 264);

    auto tr_z_0_z_z_xzz = pbuffer.data(idx_op_geom_110_pf + 265);

    auto tr_z_0_z_z_yyy = pbuffer.data(idx_op_geom_110_pf + 266);

    auto tr_z_0_z_z_yyz = pbuffer.data(idx_op_geom_110_pf + 267);

    auto tr_z_0_z_z_yzz = pbuffer.data(idx_op_geom_110_pf + 268);

    auto tr_z_0_z_z_zzz = pbuffer.data(idx_op_geom_110_pf + 269);

    #pragma omp simd aligned(tr_0_xx, tr_0_xxxz, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_z_0_z_z_xxx, tr_z_0_z_z_xxy, tr_z_0_z_z_xxz, tr_z_0_z_z_xyy, tr_z_0_z_z_xyz, tr_z_0_z_z_xzz, tr_z_0_z_z_yyy, tr_z_0_z_z_yyz, tr_z_0_z_z_yzz, tr_z_0_z_z_zzz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, tr_zz_xx, tr_zz_xxxz, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_z_xxx[i] = -2.0 * tr_0_xxxz[i] * tke_0 - 6.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_zz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xxy[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 6.0 * tr_z_xxy[i] * tbe_0 + 4.0 * tr_zz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xxz[i] = tr_0_xx[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 6.0 * tr_z_xxz[i] * tbe_0 - 2.0 * tr_zz_xx[i] * tbe_0 + 4.0 * tr_zz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xyy[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 6.0 * tr_z_xyy[i] * tbe_0 + 4.0 * tr_zz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xyz[i] = tr_0_xy[i] - 2.0 * tr_0_xyzz[i] * tke_0 - 6.0 * tr_z_xyz[i] * tbe_0 - 2.0 * tr_zz_xy[i] * tbe_0 + 4.0 * tr_zz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xzz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xzzz[i] * tke_0 - 6.0 * tr_z_xzz[i] * tbe_0 - 4.0 * tr_zz_xz[i] * tbe_0 + 4.0 * tr_zz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_yyy[i] = -2.0 * tr_0_yyyz[i] * tke_0 - 6.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_zz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_yyz[i] = tr_0_yy[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 6.0 * tr_z_yyz[i] * tbe_0 - 2.0 * tr_zz_yy[i] * tbe_0 + 4.0 * tr_zz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_yzz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yzzz[i] * tke_0 - 6.0 * tr_z_yzz[i] * tbe_0 - 4.0 * tr_zz_yz[i] * tbe_0 + 4.0 * tr_zz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_zzz[i] = 3.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tke_0 - 6.0 * tr_z_zzz[i] * tbe_0 - 6.0 * tr_zz_zz[i] * tbe_0 + 4.0 * tr_zz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_zzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

