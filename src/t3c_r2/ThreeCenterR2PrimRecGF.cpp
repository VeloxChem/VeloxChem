#include "ThreeCenterR2PrimRecGF.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_gf(CSimdArray<double>& pbuffer, 
                const size_t idx_g_gf,
                const size_t idx_df,
                const size_t idx_fd,
                const size_t idx_ff,
                const size_t idx_gp,
                const size_t idx_gd,
                const size_t idx_gf,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_df);

    auto ts_xx_xxy = pbuffer.data(idx_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_df + 9);

    auto ts_xy_xxx = pbuffer.data(idx_df + 10);

    auto ts_xy_xxy = pbuffer.data(idx_df + 11);

    auto ts_xy_xxz = pbuffer.data(idx_df + 12);

    auto ts_xy_xyy = pbuffer.data(idx_df + 13);

    auto ts_xy_xyz = pbuffer.data(idx_df + 14);

    auto ts_xy_xzz = pbuffer.data(idx_df + 15);

    auto ts_xy_yyy = pbuffer.data(idx_df + 16);

    auto ts_xy_yyz = pbuffer.data(idx_df + 17);

    auto ts_xy_yzz = pbuffer.data(idx_df + 18);

    auto ts_xy_zzz = pbuffer.data(idx_df + 19);

    auto ts_xz_xxx = pbuffer.data(idx_df + 20);

    auto ts_xz_xxy = pbuffer.data(idx_df + 21);

    auto ts_xz_xxz = pbuffer.data(idx_df + 22);

    auto ts_xz_xyy = pbuffer.data(idx_df + 23);

    auto ts_xz_xyz = pbuffer.data(idx_df + 24);

    auto ts_xz_xzz = pbuffer.data(idx_df + 25);

    auto ts_xz_yyy = pbuffer.data(idx_df + 26);

    auto ts_xz_yyz = pbuffer.data(idx_df + 27);

    auto ts_xz_yzz = pbuffer.data(idx_df + 28);

    auto ts_xz_zzz = pbuffer.data(idx_df + 29);

    auto ts_yy_xxx = pbuffer.data(idx_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_df + 39);

    auto ts_yz_xxx = pbuffer.data(idx_df + 40);

    auto ts_yz_xxy = pbuffer.data(idx_df + 41);

    auto ts_yz_xxz = pbuffer.data(idx_df + 42);

    auto ts_yz_xyy = pbuffer.data(idx_df + 43);

    auto ts_yz_xyz = pbuffer.data(idx_df + 44);

    auto ts_yz_xzz = pbuffer.data(idx_df + 45);

    auto ts_yz_yyy = pbuffer.data(idx_df + 46);

    auto ts_yz_yyz = pbuffer.data(idx_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_df + 48);

    auto ts_yz_zzz = pbuffer.data(idx_df + 49);

    auto ts_zz_xxx = pbuffer.data(idx_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_df + 59);

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_fd);

    auto ts_xxx_xy = pbuffer.data(idx_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_fd + 6);

    auto ts_xxy_xy = pbuffer.data(idx_fd + 7);

    auto ts_xxy_xz = pbuffer.data(idx_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_fd + 10);

    auto ts_xxy_zz = pbuffer.data(idx_fd + 11);

    auto ts_xxz_xx = pbuffer.data(idx_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_fd + 14);

    auto ts_xxz_yy = pbuffer.data(idx_fd + 15);

    auto ts_xxz_yz = pbuffer.data(idx_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_fd + 17);

    auto ts_xyy_xx = pbuffer.data(idx_fd + 18);

    auto ts_xyy_xy = pbuffer.data(idx_fd + 19);

    auto ts_xyy_xz = pbuffer.data(idx_fd + 20);

    auto ts_xyy_yy = pbuffer.data(idx_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_fd + 23);

    auto ts_xyz_xx = pbuffer.data(idx_fd + 24);

    auto ts_xyz_xy = pbuffer.data(idx_fd + 25);

    auto ts_xyz_xz = pbuffer.data(idx_fd + 26);

    auto ts_xyz_yy = pbuffer.data(idx_fd + 27);

    auto ts_xyz_yz = pbuffer.data(idx_fd + 28);

    auto ts_xyz_zz = pbuffer.data(idx_fd + 29);

    auto ts_xzz_xx = pbuffer.data(idx_fd + 30);

    auto ts_xzz_xy = pbuffer.data(idx_fd + 31);

    auto ts_xzz_xz = pbuffer.data(idx_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_fd + 41);

    auto ts_yyz_xx = pbuffer.data(idx_fd + 42);

    auto ts_yyz_xy = pbuffer.data(idx_fd + 43);

    auto ts_yyz_xz = pbuffer.data(idx_fd + 44);

    auto ts_yyz_yy = pbuffer.data(idx_fd + 45);

    auto ts_yyz_yz = pbuffer.data(idx_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_fd + 47);

    auto ts_yzz_xx = pbuffer.data(idx_fd + 48);

    auto ts_yzz_xy = pbuffer.data(idx_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_fd + 59);

    // Set up components of auxiliary buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ff + 9);

    auto ts_xxy_xxx = pbuffer.data(idx_ff + 10);

    auto ts_xxy_xxy = pbuffer.data(idx_ff + 11);

    auto ts_xxy_xxz = pbuffer.data(idx_ff + 12);

    auto ts_xxy_xyy = pbuffer.data(idx_ff + 13);

    auto ts_xxy_xyz = pbuffer.data(idx_ff + 14);

    auto ts_xxy_xzz = pbuffer.data(idx_ff + 15);

    auto ts_xxy_yyy = pbuffer.data(idx_ff + 16);

    auto ts_xxy_yyz = pbuffer.data(idx_ff + 17);

    auto ts_xxy_yzz = pbuffer.data(idx_ff + 18);

    auto ts_xxy_zzz = pbuffer.data(idx_ff + 19);

    auto ts_xxz_xxx = pbuffer.data(idx_ff + 20);

    auto ts_xxz_xxy = pbuffer.data(idx_ff + 21);

    auto ts_xxz_xxz = pbuffer.data(idx_ff + 22);

    auto ts_xxz_xyy = pbuffer.data(idx_ff + 23);

    auto ts_xxz_xyz = pbuffer.data(idx_ff + 24);

    auto ts_xxz_xzz = pbuffer.data(idx_ff + 25);

    auto ts_xxz_yyy = pbuffer.data(idx_ff + 26);

    auto ts_xxz_yyz = pbuffer.data(idx_ff + 27);

    auto ts_xxz_yzz = pbuffer.data(idx_ff + 28);

    auto ts_xxz_zzz = pbuffer.data(idx_ff + 29);

    auto ts_xyy_xxx = pbuffer.data(idx_ff + 30);

    auto ts_xyy_xxy = pbuffer.data(idx_ff + 31);

    auto ts_xyy_xxz = pbuffer.data(idx_ff + 32);

    auto ts_xyy_xyy = pbuffer.data(idx_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ff + 34);

    auto ts_xyy_xzz = pbuffer.data(idx_ff + 35);

    auto ts_xyy_yyy = pbuffer.data(idx_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ff + 38);

    auto ts_xyy_zzz = pbuffer.data(idx_ff + 39);

    auto ts_xyz_xxx = pbuffer.data(idx_ff + 40);

    auto ts_xyz_xxy = pbuffer.data(idx_ff + 41);

    auto ts_xyz_xxz = pbuffer.data(idx_ff + 42);

    auto ts_xyz_xyy = pbuffer.data(idx_ff + 43);

    auto ts_xyz_xyz = pbuffer.data(idx_ff + 44);

    auto ts_xyz_xzz = pbuffer.data(idx_ff + 45);

    auto ts_xyz_yyy = pbuffer.data(idx_ff + 46);

    auto ts_xyz_yyz = pbuffer.data(idx_ff + 47);

    auto ts_xyz_yzz = pbuffer.data(idx_ff + 48);

    auto ts_xyz_zzz = pbuffer.data(idx_ff + 49);

    auto ts_xzz_xxx = pbuffer.data(idx_ff + 50);

    auto ts_xzz_xxy = pbuffer.data(idx_ff + 51);

    auto ts_xzz_xxz = pbuffer.data(idx_ff + 52);

    auto ts_xzz_xyy = pbuffer.data(idx_ff + 53);

    auto ts_xzz_xyz = pbuffer.data(idx_ff + 54);

    auto ts_xzz_xzz = pbuffer.data(idx_ff + 55);

    auto ts_xzz_yyy = pbuffer.data(idx_ff + 56);

    auto ts_xzz_yyz = pbuffer.data(idx_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ff + 59);

    auto ts_yyy_xxx = pbuffer.data(idx_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ff + 69);

    auto ts_yyz_xxx = pbuffer.data(idx_ff + 70);

    auto ts_yyz_xxy = pbuffer.data(idx_ff + 71);

    auto ts_yyz_xxz = pbuffer.data(idx_ff + 72);

    auto ts_yyz_xyy = pbuffer.data(idx_ff + 73);

    auto ts_yyz_xyz = pbuffer.data(idx_ff + 74);

    auto ts_yyz_xzz = pbuffer.data(idx_ff + 75);

    auto ts_yyz_yyy = pbuffer.data(idx_ff + 76);

    auto ts_yyz_yyz = pbuffer.data(idx_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ff + 79);

    auto ts_yzz_xxx = pbuffer.data(idx_ff + 80);

    auto ts_yzz_xxy = pbuffer.data(idx_ff + 81);

    auto ts_yzz_xxz = pbuffer.data(idx_ff + 82);

    auto ts_yzz_xyy = pbuffer.data(idx_ff + 83);

    auto ts_yzz_xyz = pbuffer.data(idx_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ff + 86);

    auto ts_yzz_yyz = pbuffer.data(idx_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ff + 89);

    auto ts_zzz_xxx = pbuffer.data(idx_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ff + 99);

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_gp);

    auto ts_xxxx_y = pbuffer.data(idx_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_gp + 4);

    auto ts_xxxy_z = pbuffer.data(idx_gp + 5);

    auto ts_xxxz_x = pbuffer.data(idx_gp + 6);

    auto ts_xxxz_y = pbuffer.data(idx_gp + 7);

    auto ts_xxxz_z = pbuffer.data(idx_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_gp + 11);

    auto ts_xxyz_x = pbuffer.data(idx_gp + 12);

    auto ts_xxyz_y = pbuffer.data(idx_gp + 13);

    auto ts_xxyz_z = pbuffer.data(idx_gp + 14);

    auto ts_xxzz_x = pbuffer.data(idx_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_gp + 17);

    auto ts_xyyy_x = pbuffer.data(idx_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_gp + 20);

    auto ts_xyyz_x = pbuffer.data(idx_gp + 21);

    auto ts_xyyz_y = pbuffer.data(idx_gp + 22);

    auto ts_xyyz_z = pbuffer.data(idx_gp + 23);

    auto ts_xyzz_x = pbuffer.data(idx_gp + 24);

    auto ts_xyzz_y = pbuffer.data(idx_gp + 25);

    auto ts_xyzz_z = pbuffer.data(idx_gp + 26);

    auto ts_xzzz_x = pbuffer.data(idx_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_gp + 32);

    auto ts_yyyz_x = pbuffer.data(idx_gp + 33);

    auto ts_yyyz_y = pbuffer.data(idx_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto ts_xxxy_zz = pbuffer.data(idx_gd + 11);

    auto ts_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto ts_xxxz_yy = pbuffer.data(idx_gd + 15);

    auto ts_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto ts_xxyz_xx = pbuffer.data(idx_gd + 24);

    auto ts_xxyz_xy = pbuffer.data(idx_gd + 25);

    auto ts_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto ts_xxyz_yy = pbuffer.data(idx_gd + 27);

    auto ts_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto ts_xxyz_zz = pbuffer.data(idx_gd + 29);

    auto ts_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto ts_xyyy_xx = pbuffer.data(idx_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto ts_xyyy_xz = pbuffer.data(idx_gd + 38);

    auto ts_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto ts_xyyz_xx = pbuffer.data(idx_gd + 42);

    auto ts_xyyz_xy = pbuffer.data(idx_gd + 43);

    auto ts_xyyz_xz = pbuffer.data(idx_gd + 44);

    auto ts_xyyz_yy = pbuffer.data(idx_gd + 45);

    auto ts_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto ts_xyzz_xx = pbuffer.data(idx_gd + 48);

    auto ts_xyzz_xy = pbuffer.data(idx_gd + 49);

    auto ts_xyzz_xz = pbuffer.data(idx_gd + 50);

    auto ts_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto ts_xyzz_zz = pbuffer.data(idx_gd + 53);

    auto ts_xzzz_xx = pbuffer.data(idx_gd + 54);

    auto ts_xzzz_xy = pbuffer.data(idx_gd + 55);

    auto ts_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto ts_yyyz_xx = pbuffer.data(idx_gd + 66);

    auto ts_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_gf + 9);

    auto ts_xxxy_xxx = pbuffer.data(idx_gf + 10);

    auto ts_xxxy_xxy = pbuffer.data(idx_gf + 11);

    auto ts_xxxy_xxz = pbuffer.data(idx_gf + 12);

    auto ts_xxxy_xyy = pbuffer.data(idx_gf + 13);

    auto ts_xxxy_xyz = pbuffer.data(idx_gf + 14);

    auto ts_xxxy_xzz = pbuffer.data(idx_gf + 15);

    auto ts_xxxy_yyy = pbuffer.data(idx_gf + 16);

    auto ts_xxxy_yyz = pbuffer.data(idx_gf + 17);

    auto ts_xxxy_yzz = pbuffer.data(idx_gf + 18);

    auto ts_xxxy_zzz = pbuffer.data(idx_gf + 19);

    auto ts_xxxz_xxx = pbuffer.data(idx_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_gf + 21);

    auto ts_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto ts_xxxz_xyy = pbuffer.data(idx_gf + 23);

    auto ts_xxxz_xyz = pbuffer.data(idx_gf + 24);

    auto ts_xxxz_xzz = pbuffer.data(idx_gf + 25);

    auto ts_xxxz_yyy = pbuffer.data(idx_gf + 26);

    auto ts_xxxz_yyz = pbuffer.data(idx_gf + 27);

    auto ts_xxxz_yzz = pbuffer.data(idx_gf + 28);

    auto ts_xxxz_zzz = pbuffer.data(idx_gf + 29);

    auto ts_xxyy_xxx = pbuffer.data(idx_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_gf + 39);

    auto ts_xxyz_xxx = pbuffer.data(idx_gf + 40);

    auto ts_xxyz_xxy = pbuffer.data(idx_gf + 41);

    auto ts_xxyz_xxz = pbuffer.data(idx_gf + 42);

    auto ts_xxyz_xyy = pbuffer.data(idx_gf + 43);

    auto ts_xxyz_xyz = pbuffer.data(idx_gf + 44);

    auto ts_xxyz_xzz = pbuffer.data(idx_gf + 45);

    auto ts_xxyz_yyy = pbuffer.data(idx_gf + 46);

    auto ts_xxyz_yyz = pbuffer.data(idx_gf + 47);

    auto ts_xxyz_yzz = pbuffer.data(idx_gf + 48);

    auto ts_xxyz_zzz = pbuffer.data(idx_gf + 49);

    auto ts_xxzz_xxx = pbuffer.data(idx_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_gf + 59);

    auto ts_xyyy_xxx = pbuffer.data(idx_gf + 60);

    auto ts_xyyy_xxy = pbuffer.data(idx_gf + 61);

    auto ts_xyyy_xxz = pbuffer.data(idx_gf + 62);

    auto ts_xyyy_xyy = pbuffer.data(idx_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_gf + 64);

    auto ts_xyyy_xzz = pbuffer.data(idx_gf + 65);

    auto ts_xyyy_yyy = pbuffer.data(idx_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_gf + 69);

    auto ts_xyyz_xxx = pbuffer.data(idx_gf + 70);

    auto ts_xyyz_xxy = pbuffer.data(idx_gf + 71);

    auto ts_xyyz_xxz = pbuffer.data(idx_gf + 72);

    auto ts_xyyz_xyy = pbuffer.data(idx_gf + 73);

    auto ts_xyyz_xyz = pbuffer.data(idx_gf + 74);

    auto ts_xyyz_xzz = pbuffer.data(idx_gf + 75);

    auto ts_xyyz_yyy = pbuffer.data(idx_gf + 76);

    auto ts_xyyz_yyz = pbuffer.data(idx_gf + 77);

    auto ts_xyyz_yzz = pbuffer.data(idx_gf + 78);

    auto ts_xyyz_zzz = pbuffer.data(idx_gf + 79);

    auto ts_xyzz_xxx = pbuffer.data(idx_gf + 80);

    auto ts_xyzz_xxy = pbuffer.data(idx_gf + 81);

    auto ts_xyzz_xxz = pbuffer.data(idx_gf + 82);

    auto ts_xyzz_xyy = pbuffer.data(idx_gf + 83);

    auto ts_xyzz_xyz = pbuffer.data(idx_gf + 84);

    auto ts_xyzz_xzz = pbuffer.data(idx_gf + 85);

    auto ts_xyzz_yyy = pbuffer.data(idx_gf + 86);

    auto ts_xyzz_yyz = pbuffer.data(idx_gf + 87);

    auto ts_xyzz_yzz = pbuffer.data(idx_gf + 88);

    auto ts_xyzz_zzz = pbuffer.data(idx_gf + 89);

    auto ts_xzzz_xxx = pbuffer.data(idx_gf + 90);

    auto ts_xzzz_xxy = pbuffer.data(idx_gf + 91);

    auto ts_xzzz_xxz = pbuffer.data(idx_gf + 92);

    auto ts_xzzz_xyy = pbuffer.data(idx_gf + 93);

    auto ts_xzzz_xyz = pbuffer.data(idx_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_gf + 109);

    auto ts_yyyz_xxx = pbuffer.data(idx_gf + 110);

    auto ts_yyyz_xxy = pbuffer.data(idx_gf + 111);

    auto ts_yyyz_xxz = pbuffer.data(idx_gf + 112);

    auto ts_yyyz_xyy = pbuffer.data(idx_gf + 113);

    auto ts_yyyz_xyz = pbuffer.data(idx_gf + 114);

    auto ts_yyyz_xzz = pbuffer.data(idx_gf + 115);

    auto ts_yyyz_yyy = pbuffer.data(idx_gf + 116);

    auto ts_yyyz_yyz = pbuffer.data(idx_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_gf + 119);

    auto ts_yyzz_xxx = pbuffer.data(idx_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_gf + 129);

    auto ts_yzzz_xxx = pbuffer.data(idx_gf + 130);

    auto ts_yzzz_xxy = pbuffer.data(idx_gf + 131);

    auto ts_yzzz_xxz = pbuffer.data(idx_gf + 132);

    auto ts_yzzz_xyy = pbuffer.data(idx_gf + 133);

    auto ts_yzzz_xyz = pbuffer.data(idx_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_gf + 149);

    // Set up 0-10 components of targeted buffer : GF

    auto gr_xxxx_xxx = pbuffer.data(idx_g_gf);

    auto gr_xxxx_xxy = pbuffer.data(idx_g_gf + 1);

    auto gr_xxxx_xxz = pbuffer.data(idx_g_gf + 2);

    auto gr_xxxx_xyy = pbuffer.data(idx_g_gf + 3);

    auto gr_xxxx_xyz = pbuffer.data(idx_g_gf + 4);

    auto gr_xxxx_xzz = pbuffer.data(idx_g_gf + 5);

    auto gr_xxxx_yyy = pbuffer.data(idx_g_gf + 6);

    auto gr_xxxx_yyz = pbuffer.data(idx_g_gf + 7);

    auto gr_xxxx_yzz = pbuffer.data(idx_g_gf + 8);

    auto gr_xxxx_zzz = pbuffer.data(idx_g_gf + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_xxx, gr_xxxx_xxy, gr_xxxx_xxz, gr_xxxx_xyy, gr_xxxx_xyz, gr_xxxx_xzz, gr_xxxx_yyy, gr_xxxx_yyz, gr_xxxx_yzz, gr_xxxx_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xy, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xz, ts_xxxx_xzz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yz, ts_xxxx_yzz, ts_xxxx_z, ts_xxxx_zz, ts_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxx_xxx[i] = 12.0 * ts_xx_xxx[i] * gfe_0 + 24.0 * ts_xxx_xx[i] * gfe_0 + 8.0 * ts_xxx_xxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_x[i] * gfe_0 + 6.0 * ts_xxxx_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxx_xxx[i] * gfe_0 + ts_xxxx_xxx[i] * rgc2_0;

        gr_xxxx_xxy[i] = 12.0 * ts_xx_xxy[i] * gfe_0 + 16.0 * ts_xxx_xy[i] * gfe_0 + 8.0 * ts_xxx_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_y[i] * gfe_0 + 4.0 * ts_xxxx_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_xxy[i] * gfe_0 + ts_xxxx_xxy[i] * rgc2_0;

        gr_xxxx_xxz[i] = 12.0 * ts_xx_xxz[i] * gfe_0 + 16.0 * ts_xxx_xz[i] * gfe_0 + 8.0 * ts_xxx_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_z[i] * gfe_0 + 4.0 * ts_xxxx_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xxz[i] * gfe_0 + ts_xxxx_xxz[i] * rgc2_0;

        gr_xxxx_xyy[i] = 12.0 * ts_xx_xyy[i] * gfe_0 + 8.0 * ts_xxx_yy[i] * gfe_0 + 8.0 * ts_xxx_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_x[i] * gfe_0 + 4.0 * ts_xxxx_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_xyy[i] * gfe_0 + ts_xxxx_xyy[i] * rgc2_0;

        gr_xxxx_xyz[i] = 12.0 * ts_xx_xyz[i] * gfe_0 + 8.0 * ts_xxx_yz[i] * gfe_0 + 8.0 * ts_xxx_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xyz[i] * gfe_0 + ts_xxxx_xyz[i] * rgc2_0;

        gr_xxxx_xzz[i] = 12.0 * ts_xx_xzz[i] * gfe_0 + 8.0 * ts_xxx_zz[i] * gfe_0 + 8.0 * ts_xxx_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_x[i] * gfe_0 + 4.0 * ts_xxxx_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xzz[i] * gfe_0 + ts_xxxx_xzz[i] * rgc2_0;

        gr_xxxx_yyy[i] = 12.0 * ts_xx_yyy[i] * gfe_0 + 8.0 * ts_xxx_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_y[i] * gfe_0 + 6.0 * ts_xxxx_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_yyy[i] * gfe_0 + ts_xxxx_yyy[i] * rgc2_0;

        gr_xxxx_yyz[i] = 12.0 * ts_xx_yyz[i] * gfe_0 + 8.0 * ts_xxx_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_z[i] * gfe_0 + 4.0 * ts_xxxx_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_yyz[i] * gfe_0 + ts_xxxx_yyz[i] * rgc2_0;

        gr_xxxx_yzz[i] = 12.0 * ts_xx_yzz[i] * gfe_0 + 8.0 * ts_xxx_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_y[i] * gfe_0 + 4.0 * ts_xxxx_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_yzz[i] * gfe_0 + ts_xxxx_yzz[i] * rgc2_0;

        gr_xxxx_zzz[i] = 12.0 * ts_xx_zzz[i] * gfe_0 + 8.0 * ts_xxx_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxxx_z[i] * gfe_0 + 6.0 * ts_xxxx_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_zzz[i] * gfe_0 + ts_xxxx_zzz[i] * rgc2_0;
    }

    // Set up 10-20 components of targeted buffer : GF

    auto gr_xxxy_xxx = pbuffer.data(idx_g_gf + 10);

    auto gr_xxxy_xxy = pbuffer.data(idx_g_gf + 11);

    auto gr_xxxy_xxz = pbuffer.data(idx_g_gf + 12);

    auto gr_xxxy_xyy = pbuffer.data(idx_g_gf + 13);

    auto gr_xxxy_xyz = pbuffer.data(idx_g_gf + 14);

    auto gr_xxxy_xzz = pbuffer.data(idx_g_gf + 15);

    auto gr_xxxy_yyy = pbuffer.data(idx_g_gf + 16);

    auto gr_xxxy_yyz = pbuffer.data(idx_g_gf + 17);

    auto gr_xxxy_yzz = pbuffer.data(idx_g_gf + 18);

    auto gr_xxxy_zzz = pbuffer.data(idx_g_gf + 19);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_xxx, gr_xxxy_xxy, gr_xxxy_xxz, gr_xxxy_xyy, gr_xxxy_xyz, gr_xxxy_xzz, gr_xxxy_yyy, gr_xxxy_yyz, gr_xxxy_yzz, gr_xxxy_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xxx, ts_xxxy_xxy, ts_xxxy_xxz, ts_xxxy_xy, ts_xxxy_xyy, ts_xxxy_xyz, ts_xxxy_xz, ts_xxxy_xzz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yz, ts_xxxy_yzz, ts_xxxy_z, ts_xxxy_zz, ts_xxxy_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxy_xxx[i] = 6.0 * ts_xy_xxx[i] * gfe_0 + 18.0 * ts_xxy_xx[i] * gfe_0 + 6.0 * ts_xxy_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxy_x[i] * gfe_0 + 6.0 * ts_xxxy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxy_xxx[i] * gfe_0 + ts_xxxy_xxx[i] * rgc2_0;

        gr_xxxy_xxy[i] = 6.0 * ts_xy_xxy[i] * gfe_0 + 12.0 * ts_xxy_xy[i] * gfe_0 + 6.0 * ts_xxy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe_0 + 2.0 * ts_xxx_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_y[i] * gfe_0 + 4.0 * ts_xxxy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_xxy[i] * gfe_0 + ts_xxxy_xxy[i] * rgc2_0;

        gr_xxxy_xxz[i] = 6.0 * ts_xy_xxz[i] * gfe_0 + 12.0 * ts_xxy_xz[i] * gfe_0 + 6.0 * ts_xxy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_z[i] * gfe_0 + 4.0 * ts_xxxy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xxz[i] * gfe_0 + ts_xxxy_xxz[i] * rgc2_0;

        gr_xxxy_xyy[i] = 6.0 * ts_xy_xyy[i] * gfe_0 + 6.0 * ts_xxy_yy[i] * gfe_0 + 6.0 * ts_xxy_xyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_xy[i] * gfe_0 + 2.0 * ts_xxx_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_x[i] * gfe_0 + 4.0 * ts_xxxy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_xyy[i] * gfe_0 + ts_xxxy_xyy[i] * rgc2_0;

        gr_xxxy_xyz[i] = 6.0 * ts_xy_xyz[i] * gfe_0 + 6.0 * ts_xxy_yz[i] * gfe_0 + 6.0 * ts_xxy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xz[i] * gfe_0 + 2.0 * ts_xxx_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xyz[i] * gfe_0 + ts_xxxy_xyz[i] * rgc2_0;

        gr_xxxy_xzz[i] = 6.0 * ts_xy_xzz[i] * gfe_0 + 6.0 * ts_xxy_zz[i] * gfe_0 + 6.0 * ts_xxy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_x[i] * gfe_0 + 4.0 * ts_xxxy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xzz[i] * gfe_0 + ts_xxxy_xzz[i] * rgc2_0;

        gr_xxxy_yyy[i] = 6.0 * ts_xy_yyy[i] * gfe_0 + 6.0 * ts_xxy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_yy[i] * gfe_0 + 2.0 * ts_xxx_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxy_y[i] * gfe_0 + 6.0 * ts_xxxy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_yyy[i] * gfe_0 + ts_xxxy_yyy[i] * rgc2_0;

        gr_xxxy_yyz[i] = 6.0 * ts_xy_yyz[i] * gfe_0 + 6.0 * ts_xxy_yyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_yz[i] * gfe_0 + 2.0 * ts_xxx_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_z[i] * gfe_0 + 4.0 * ts_xxxy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_yyz[i] * gfe_0 + ts_xxxy_yyz[i] * rgc2_0;

        gr_xxxy_yzz[i] = 6.0 * ts_xy_yzz[i] * gfe_0 + 6.0 * ts_xxy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zz[i] * gfe_0 + 2.0 * ts_xxx_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_y[i] * gfe_0 + 4.0 * ts_xxxy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_yzz[i] * gfe_0 + ts_xxxy_yzz[i] * rgc2_0;

        gr_xxxy_zzz[i] = 6.0 * ts_xy_zzz[i] * gfe_0 + 6.0 * ts_xxy_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxxy_z[i] * gfe_0 + 6.0 * ts_xxxy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_zzz[i] * gfe_0 + ts_xxxy_zzz[i] * rgc2_0;
    }

    // Set up 20-30 components of targeted buffer : GF

    auto gr_xxxz_xxx = pbuffer.data(idx_g_gf + 20);

    auto gr_xxxz_xxy = pbuffer.data(idx_g_gf + 21);

    auto gr_xxxz_xxz = pbuffer.data(idx_g_gf + 22);

    auto gr_xxxz_xyy = pbuffer.data(idx_g_gf + 23);

    auto gr_xxxz_xyz = pbuffer.data(idx_g_gf + 24);

    auto gr_xxxz_xzz = pbuffer.data(idx_g_gf + 25);

    auto gr_xxxz_yyy = pbuffer.data(idx_g_gf + 26);

    auto gr_xxxz_yyz = pbuffer.data(idx_g_gf + 27);

    auto gr_xxxz_yzz = pbuffer.data(idx_g_gf + 28);

    auto gr_xxxz_zzz = pbuffer.data(idx_g_gf + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_xxx, gr_xxxz_xxy, gr_xxxz_xxz, gr_xxxz_xyy, gr_xxxz_xyz, gr_xxxz_xzz, gr_xxxz_yyy, gr_xxxz_yyz, gr_xxxz_yzz, gr_xxxz_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xxz, ts_xxxz_xy, ts_xxxz_xyy, ts_xxxz_xyz, ts_xxxz_xz, ts_xxxz_xzz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yyy, ts_xxxz_yyz, ts_xxxz_yz, ts_xxxz_yzz, ts_xxxz_z, ts_xxxz_zz, ts_xxxz_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxz_xxx[i] = 6.0 * ts_xz_xxx[i] * gfe_0 + 18.0 * ts_xxz_xx[i] * gfe_0 + 6.0 * ts_xxz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxxz_x[i] * gfe_0 + 6.0 * ts_xxxz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxz_xxx[i] * gfe_0 + ts_xxxz_xxx[i] * rgc2_0;

        gr_xxxz_xxy[i] = 6.0 * ts_xz_xxy[i] * gfe_0 + 12.0 * ts_xxz_xy[i] * gfe_0 + 6.0 * ts_xxz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_y[i] * gfe_0 + 4.0 * ts_xxxz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_xxy[i] * gfe_0 + ts_xxxz_xxy[i] * rgc2_0;

        gr_xxxz_xxz[i] = 6.0 * ts_xz_xxz[i] * gfe_0 + 12.0 * ts_xxz_xz[i] * gfe_0 + 6.0 * ts_xxz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe_0 + 2.0 * ts_xxx_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_z[i] * gfe_0 + 4.0 * ts_xxxz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xxz[i] * gfe_0 + ts_xxxz_xxz[i] * rgc2_0;

        gr_xxxz_xyy[i] = 6.0 * ts_xz_xyy[i] * gfe_0 + 6.0 * ts_xxz_yy[i] * gfe_0 + 6.0 * ts_xxz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_x[i] * gfe_0 + 4.0 * ts_xxxz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_xyy[i] * gfe_0 + ts_xxxz_xyy[i] * rgc2_0;

        gr_xxxz_xyz[i] = 6.0 * ts_xz_xyz[i] * gfe_0 + 6.0 * ts_xxz_yz[i] * gfe_0 + 6.0 * ts_xxz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xy[i] * gfe_0 + 2.0 * ts_xxx_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xyz[i] * gfe_0 + ts_xxxz_xyz[i] * rgc2_0;

        gr_xxxz_xzz[i] = 6.0 * ts_xz_xzz[i] * gfe_0 + 6.0 * ts_xxz_zz[i] * gfe_0 + 6.0 * ts_xxz_xzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_xz[i] * gfe_0 + 2.0 * ts_xxx_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_x[i] * gfe_0 + 4.0 * ts_xxxz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xzz[i] * gfe_0 + ts_xxxz_xzz[i] * rgc2_0;

        gr_xxxz_yyy[i] = 6.0 * ts_xz_yyy[i] * gfe_0 + 6.0 * ts_xxz_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxxz_y[i] * gfe_0 + 6.0 * ts_xxxz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_yyy[i] * gfe_0 + ts_xxxz_yyy[i] * rgc2_0;

        gr_xxxz_yyz[i] = 6.0 * ts_xz_yyz[i] * gfe_0 + 6.0 * ts_xxz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yy[i] * gfe_0 + 2.0 * ts_xxx_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_z[i] * gfe_0 + 4.0 * ts_xxxz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_yyz[i] * gfe_0 + ts_xxxz_yyz[i] * rgc2_0;

        gr_xxxz_yzz[i] = 6.0 * ts_xz_yzz[i] * gfe_0 + 6.0 * ts_xxz_yzz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_yz[i] * gfe_0 + 2.0 * ts_xxx_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_y[i] * gfe_0 + 4.0 * ts_xxxz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_yzz[i] * gfe_0 + ts_xxxz_yzz[i] * rgc2_0;

        gr_xxxz_zzz[i] = 6.0 * ts_xz_zzz[i] * gfe_0 + 6.0 * ts_xxz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxx_zz[i] * gfe_0 + 2.0 * ts_xxx_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxxz_z[i] * gfe_0 + 6.0 * ts_xxxz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_zzz[i] * gfe_0 + ts_xxxz_zzz[i] * rgc2_0;
    }

    // Set up 30-40 components of targeted buffer : GF

    auto gr_xxyy_xxx = pbuffer.data(idx_g_gf + 30);

    auto gr_xxyy_xxy = pbuffer.data(idx_g_gf + 31);

    auto gr_xxyy_xxz = pbuffer.data(idx_g_gf + 32);

    auto gr_xxyy_xyy = pbuffer.data(idx_g_gf + 33);

    auto gr_xxyy_xyz = pbuffer.data(idx_g_gf + 34);

    auto gr_xxyy_xzz = pbuffer.data(idx_g_gf + 35);

    auto gr_xxyy_yyy = pbuffer.data(idx_g_gf + 36);

    auto gr_xxyy_yyz = pbuffer.data(idx_g_gf + 37);

    auto gr_xxyy_yzz = pbuffer.data(idx_g_gf + 38);

    auto gr_xxyy_zzz = pbuffer.data(idx_g_gf + 39);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_xxx, gr_xxyy_xxy, gr_xxyy_xxz, gr_xxyy_xyy, gr_xxyy_xyz, gr_xxyy_xzz, gr_xxyy_yyy, gr_xxyy_yyz, gr_xxyy_yzz, gr_xxyy_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xxx, ts_xxyy_xxy, ts_xxyy_xxz, ts_xxyy_xy, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_xz, ts_xxyy_xzz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yz, ts_xxyy_yzz, ts_xxyy_z, ts_xxyy_zz, ts_xxyy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxyy_xxx[i] = 2.0 * ts_yy_xxx[i] * gfe_0 + 12.0 * ts_xyy_xx[i] * gfe_0 + 4.0 * ts_xyy_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe_0 + 4.0 * ts_xxy_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyy_x[i] * gfe_0 + 6.0 * ts_xxyy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyy_xxx[i] * gfe_0 + ts_xxyy_xxx[i] * rgc2_0;

        gr_xxyy_xxy[i] = 2.0 * ts_yy_xxy[i] * gfe_0 + 8.0 * ts_xyy_xy[i] * gfe_0 + 4.0 * ts_xyy_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxy[i] * gfe_0 + 4.0 * ts_xxy_xx[i] * gfe_0 + 4.0 * ts_xxy_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_y[i] * gfe_0 + 4.0 * ts_xxyy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_xxy[i] * gfe_0 + ts_xxyy_xxy[i] * rgc2_0;

        gr_xxyy_xxz[i] = 2.0 * ts_yy_xxz[i] * gfe_0 + 8.0 * ts_xyy_xz[i] * gfe_0 + 4.0 * ts_xyy_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxz[i] * gfe_0 + 4.0 * ts_xxy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_z[i] * gfe_0 + 4.0 * ts_xxyy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xxz[i] * gfe_0 + ts_xxyy_xxz[i] * rgc2_0;

        gr_xxyy_xyy[i] = 2.0 * ts_yy_xyy[i] * gfe_0 + 4.0 * ts_xyy_yy[i] * gfe_0 + 4.0 * ts_xyy_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyy[i] * gfe_0 + 8.0 * ts_xxy_xy[i] * gfe_0 + 4.0 * ts_xxy_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_x[i] * gfe_0 + 4.0 * ts_xxyy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_xyy[i] * gfe_0 + ts_xxyy_xyy[i] * rgc2_0;

        gr_xxyy_xyz[i] = 2.0 * ts_yy_xyz[i] * gfe_0 + 4.0 * ts_xyy_yz[i] * gfe_0 + 4.0 * ts_xyy_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyz[i] * gfe_0 + 4.0 * ts_xxy_xz[i] * gfe_0 + 4.0 * ts_xxy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xyz[i] * gfe_0 + ts_xxyy_xyz[i] * rgc2_0;

        gr_xxyy_xzz[i] = 2.0 * ts_yy_xzz[i] * gfe_0 + 4.0 * ts_xyy_zz[i] * gfe_0 + 4.0 * ts_xyy_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzz[i] * gfe_0 + 4.0 * ts_xxy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_x[i] * gfe_0 + 4.0 * ts_xxyy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xzz[i] * gfe_0 + ts_xxyy_xzz[i] * rgc2_0;

        gr_xxyy_yyy[i] = 2.0 * ts_yy_yyy[i] * gfe_0 + 4.0 * ts_xyy_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyy[i] * gfe_0 + 12.0 * ts_xxy_yy[i] * gfe_0 + 4.0 * ts_xxy_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyy_y[i] * gfe_0 + 6.0 * ts_xxyy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_yyy[i] * gfe_0 + ts_xxyy_yyy[i] * rgc2_0;

        gr_xxyy_yyz[i] = 2.0 * ts_yy_yyz[i] * gfe_0 + 4.0 * ts_xyy_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyz[i] * gfe_0 + 8.0 * ts_xxy_yz[i] * gfe_0 + 4.0 * ts_xxy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_z[i] * gfe_0 + 4.0 * ts_xxyy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_yyz[i] * gfe_0 + ts_xxyy_yyz[i] * rgc2_0;

        gr_xxyy_yzz[i] = 2.0 * ts_yy_yzz[i] * gfe_0 + 4.0 * ts_xyy_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yzz[i] * gfe_0 + 4.0 * ts_xxy_zz[i] * gfe_0 + 4.0 * ts_xxy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_y[i] * gfe_0 + 4.0 * ts_xxyy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_yzz[i] * gfe_0 + ts_xxyy_yzz[i] * rgc2_0;

        gr_xxyy_zzz[i] = 2.0 * ts_yy_zzz[i] * gfe_0 + 4.0 * ts_xyy_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzz[i] * gfe_0 + 4.0 * ts_xxy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxyy_z[i] * gfe_0 + 6.0 * ts_xxyy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_zzz[i] * gfe_0 + ts_xxyy_zzz[i] * rgc2_0;
    }

    // Set up 40-50 components of targeted buffer : GF

    auto gr_xxyz_xxx = pbuffer.data(idx_g_gf + 40);

    auto gr_xxyz_xxy = pbuffer.data(idx_g_gf + 41);

    auto gr_xxyz_xxz = pbuffer.data(idx_g_gf + 42);

    auto gr_xxyz_xyy = pbuffer.data(idx_g_gf + 43);

    auto gr_xxyz_xyz = pbuffer.data(idx_g_gf + 44);

    auto gr_xxyz_xzz = pbuffer.data(idx_g_gf + 45);

    auto gr_xxyz_yyy = pbuffer.data(idx_g_gf + 46);

    auto gr_xxyz_yyz = pbuffer.data(idx_g_gf + 47);

    auto gr_xxyz_yzz = pbuffer.data(idx_g_gf + 48);

    auto gr_xxyz_zzz = pbuffer.data(idx_g_gf + 49);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_xxx, gr_xxyz_xxy, gr_xxyz_xxz, gr_xxyz_xyy, gr_xxyz_xyz, gr_xxyz_xzz, gr_xxyz_yyy, gr_xxyz_yyz, gr_xxyz_yzz, gr_xxyz_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xxx, ts_xxyz_xxy, ts_xxyz_xxz, ts_xxyz_xy, ts_xxyz_xyy, ts_xxyz_xyz, ts_xxyz_xz, ts_xxyz_xzz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yyy, ts_xxyz_yyz, ts_xxyz_yz, ts_xxyz_yzz, ts_xxyz_z, ts_xxyz_zz, ts_xxyz_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxyz_xxx[i] = 2.0 * ts_yz_xxx[i] * gfe_0 + 12.0 * ts_xyz_xx[i] * gfe_0 + 4.0 * ts_xyz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxyz_x[i] * gfe_0 + 6.0 * ts_xxyz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyz_xxx[i] * gfe_0 + ts_xxyz_xxx[i] * rgc2_0;

        gr_xxyz_xxy[i] = 2.0 * ts_yz_xxy[i] * gfe_0 + 8.0 * ts_xyz_xy[i] * gfe_0 + 4.0 * ts_xyz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xx[i] * gfe_0 + 2.0 * ts_xxz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_y[i] * gfe_0 + 4.0 * ts_xxyz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_xxy[i] * gfe_0 + ts_xxyz_xxy[i] * rgc2_0;

        gr_xxyz_xxz[i] = 2.0 * ts_yz_xxz[i] * gfe_0 + 8.0 * ts_xyz_xz[i] * gfe_0 + 4.0 * ts_xyz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xx[i] * gfe_0 + 2.0 * ts_xxy_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_z[i] * gfe_0 + 4.0 * ts_xxyz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xxz[i] * gfe_0 + ts_xxyz_xxz[i] * rgc2_0;

        gr_xxyz_xyy[i] = 2.0 * ts_yz_xyy[i] * gfe_0 + 4.0 * ts_xyz_yy[i] * gfe_0 + 4.0 * ts_xyz_xyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxz_xy[i] * gfe_0 + 2.0 * ts_xxz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_x[i] * gfe_0 + 4.0 * ts_xxyz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_xyy[i] * gfe_0 + ts_xxyz_xyy[i] * rgc2_0;

        gr_xxyz_xyz[i] = 2.0 * ts_yz_xyz[i] * gfe_0 + 4.0 * ts_xyz_yz[i] * gfe_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xz[i] * gfe_0 + 2.0 * ts_xxz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xy[i] * gfe_0 + 2.0 * ts_xxy_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xyz[i] * gfe_0 + ts_xxyz_xyz[i] * rgc2_0;

        gr_xxyz_xzz[i] = 2.0 * ts_yz_xzz[i] * gfe_0 + 4.0 * ts_xyz_zz[i] * gfe_0 + 4.0 * ts_xyz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xxy_xz[i] * gfe_0 + 2.0 * ts_xxy_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_x[i] * gfe_0 + 4.0 * ts_xxyz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xzz[i] * gfe_0 + ts_xxyz_xzz[i] * rgc2_0;

        gr_xxyz_yyy[i] = 2.0 * ts_yz_yyy[i] * gfe_0 + 4.0 * ts_xyz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xxz_yy[i] * gfe_0 + 2.0 * ts_xxz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxyz_y[i] * gfe_0 + 6.0 * ts_xxyz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_yyy[i] * gfe_0 + ts_xxyz_yyy[i] * rgc2_0;

        gr_xxyz_yyz[i] = 2.0 * ts_yz_yyz[i] * gfe_0 + 4.0 * ts_xyz_yyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxz_yz[i] * gfe_0 + 2.0 * ts_xxz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yy[i] * gfe_0 + 2.0 * ts_xxy_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_z[i] * gfe_0 + 4.0 * ts_xxyz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_yyz[i] * gfe_0 + ts_xxyz_yyz[i] * rgc2_0;

        gr_xxyz_yzz[i] = 2.0 * ts_yz_yzz[i] * gfe_0 + 4.0 * ts_xyz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_zz[i] * gfe_0 + 2.0 * ts_xxz_yzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xxy_yz[i] * gfe_0 + 2.0 * ts_xxy_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_y[i] * gfe_0 + 4.0 * ts_xxyz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_yzz[i] * gfe_0 + ts_xxyz_yzz[i] * rgc2_0;

        gr_xxyz_zzz[i] = 2.0 * ts_yz_zzz[i] * gfe_0 + 4.0 * ts_xyz_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xxy_zz[i] * gfe_0 + 2.0 * ts_xxy_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxyz_z[i] * gfe_0 + 6.0 * ts_xxyz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_zzz[i] * gfe_0 + ts_xxyz_zzz[i] * rgc2_0;
    }

    // Set up 50-60 components of targeted buffer : GF

    auto gr_xxzz_xxx = pbuffer.data(idx_g_gf + 50);

    auto gr_xxzz_xxy = pbuffer.data(idx_g_gf + 51);

    auto gr_xxzz_xxz = pbuffer.data(idx_g_gf + 52);

    auto gr_xxzz_xyy = pbuffer.data(idx_g_gf + 53);

    auto gr_xxzz_xyz = pbuffer.data(idx_g_gf + 54);

    auto gr_xxzz_xzz = pbuffer.data(idx_g_gf + 55);

    auto gr_xxzz_yyy = pbuffer.data(idx_g_gf + 56);

    auto gr_xxzz_yyz = pbuffer.data(idx_g_gf + 57);

    auto gr_xxzz_yzz = pbuffer.data(idx_g_gf + 58);

    auto gr_xxzz_zzz = pbuffer.data(idx_g_gf + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_xxx, gr_xxzz_xxy, gr_xxzz_xxz, gr_xxzz_xyy, gr_xxzz_xyz, gr_xxzz_xzz, gr_xxzz_yyy, gr_xxzz_yyz, gr_xxzz_yzz, gr_xxzz_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xxz, ts_xxzz_xy, ts_xxzz_xyy, ts_xxzz_xyz, ts_xxzz_xz, ts_xxzz_xzz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yz, ts_xxzz_yzz, ts_xxzz_z, ts_xxzz_zz, ts_xxzz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxzz_xxx[i] = 2.0 * ts_zz_xxx[i] * gfe_0 + 12.0 * ts_xzz_xx[i] * gfe_0 + 4.0 * ts_xzz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxx[i] * gfe_0 + 4.0 * ts_xxz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxzz_x[i] * gfe_0 + 6.0 * ts_xxzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxzz_xxx[i] * gfe_0 + ts_xxzz_xxx[i] * rgc2_0;

        gr_xxzz_xxy[i] = 2.0 * ts_zz_xxy[i] * gfe_0 + 8.0 * ts_xzz_xy[i] * gfe_0 + 4.0 * ts_xzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxy[i] * gfe_0 + 4.0 * ts_xxz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_y[i] * gfe_0 + 4.0 * ts_xxzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_xxy[i] * gfe_0 + ts_xxzz_xxy[i] * rgc2_0;

        gr_xxzz_xxz[i] = 2.0 * ts_zz_xxz[i] * gfe_0 + 8.0 * ts_xzz_xz[i] * gfe_0 + 4.0 * ts_xzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xxz[i] * gfe_0 + 4.0 * ts_xxz_xx[i] * gfe_0 + 4.0 * ts_xxz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_z[i] * gfe_0 + 4.0 * ts_xxzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xxz[i] * gfe_0 + ts_xxzz_xxz[i] * rgc2_0;

        gr_xxzz_xyy[i] = 2.0 * ts_zz_xyy[i] * gfe_0 + 4.0 * ts_xzz_yy[i] * gfe_0 + 4.0 * ts_xzz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyy[i] * gfe_0 + 4.0 * ts_xxz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_x[i] * gfe_0 + 4.0 * ts_xxzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_xyy[i] * gfe_0 + ts_xxzz_xyy[i] * rgc2_0;

        gr_xxzz_xyz[i] = 2.0 * ts_zz_xyz[i] * gfe_0 + 4.0 * ts_xzz_yz[i] * gfe_0 + 4.0 * ts_xzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xyz[i] * gfe_0 + 4.0 * ts_xxz_xy[i] * gfe_0 + 4.0 * ts_xxz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xyz[i] * gfe_0 + ts_xxzz_xyz[i] * rgc2_0;

        gr_xxzz_xzz[i] = 2.0 * ts_zz_xzz[i] * gfe_0 + 4.0 * ts_xzz_zz[i] * gfe_0 + 4.0 * ts_xzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xzz[i] * gfe_0 + 8.0 * ts_xxz_xz[i] * gfe_0 + 4.0 * ts_xxz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_x[i] * gfe_0 + 4.0 * ts_xxzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xzz[i] * gfe_0 + ts_xxzz_xzz[i] * rgc2_0;

        gr_xxzz_yyy[i] = 2.0 * ts_zz_yyy[i] * gfe_0 + 4.0 * ts_xzz_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyy[i] * gfe_0 + 4.0 * ts_xxz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxzz_y[i] * gfe_0 + 6.0 * ts_xxzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_yyy[i] * gfe_0 + ts_xxzz_yyy[i] * rgc2_0;

        gr_xxzz_yyz[i] = 2.0 * ts_zz_yyz[i] * gfe_0 + 4.0 * ts_xzz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yyz[i] * gfe_0 + 4.0 * ts_xxz_yy[i] * gfe_0 + 4.0 * ts_xxz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_z[i] * gfe_0 + 4.0 * ts_xxzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_yyz[i] * gfe_0 + ts_xxzz_yyz[i] * rgc2_0;

        gr_xxzz_yzz[i] = 2.0 * ts_zz_yzz[i] * gfe_0 + 4.0 * ts_xzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yzz[i] * gfe_0 + 8.0 * ts_xxz_yz[i] * gfe_0 + 4.0 * ts_xxz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_y[i] * gfe_0 + 4.0 * ts_xxzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_yzz[i] * gfe_0 + ts_xxzz_yzz[i] * rgc2_0;

        gr_xxzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 + 4.0 * ts_xzz_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zzz[i] * gfe_0 + 12.0 * ts_xxz_zz[i] * gfe_0 + 4.0 * ts_xxz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xxzz_z[i] * gfe_0 + 6.0 * ts_xxzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_zzz[i] * gfe_0 + ts_xxzz_zzz[i] * rgc2_0;
    }

    // Set up 60-70 components of targeted buffer : GF

    auto gr_xyyy_xxx = pbuffer.data(idx_g_gf + 60);

    auto gr_xyyy_xxy = pbuffer.data(idx_g_gf + 61);

    auto gr_xyyy_xxz = pbuffer.data(idx_g_gf + 62);

    auto gr_xyyy_xyy = pbuffer.data(idx_g_gf + 63);

    auto gr_xyyy_xyz = pbuffer.data(idx_g_gf + 64);

    auto gr_xyyy_xzz = pbuffer.data(idx_g_gf + 65);

    auto gr_xyyy_yyy = pbuffer.data(idx_g_gf + 66);

    auto gr_xyyy_yyz = pbuffer.data(idx_g_gf + 67);

    auto gr_xyyy_yzz = pbuffer.data(idx_g_gf + 68);

    auto gr_xyyy_zzz = pbuffer.data(idx_g_gf + 69);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_xxx, gr_xyyy_xxy, gr_xyyy_xxz, gr_xyyy_xyy, gr_xyyy_xyz, gr_xyyy_xzz, gr_xyyy_yyy, gr_xyyy_yyz, gr_xyyy_yzz, gr_xyyy_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xxx, ts_xyyy_xxy, ts_xyyy_xxz, ts_xyyy_xy, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_xz, ts_xyyy_xzz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yz, ts_xyyy_yzz, ts_xyyy_z, ts_xyyy_zz, ts_xyyy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyyy_xxx[i] = 6.0 * ts_yyy_xx[i] * gfe_0 + 2.0 * ts_yyy_xxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxx[i] * gfe_0 + 6.0 * ts_xyy_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyy_x[i] * gfe_0 + 6.0 * ts_xyyy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyy_xxx[i] * gfe_0 + ts_xyyy_xxx[i] * rgc2_0;

        gr_xyyy_xxy[i] = 4.0 * ts_yyy_xy[i] * gfe_0 + 2.0 * ts_yyy_xxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxy[i] * gfe_0 + 6.0 * ts_xyy_xx[i] * gfe_0 + 6.0 * ts_xyy_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_y[i] * gfe_0 + 4.0 * ts_xyyy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_xxy[i] * gfe_0 + ts_xyyy_xxy[i] * rgc2_0;

        gr_xyyy_xxz[i] = 4.0 * ts_yyy_xz[i] * gfe_0 + 2.0 * ts_yyy_xxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xxz[i] * gfe_0 + 6.0 * ts_xyy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_z[i] * gfe_0 + 4.0 * ts_xyyy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xxz[i] * gfe_0 + ts_xyyy_xxz[i] * rgc2_0;

        gr_xyyy_xyy[i] = 2.0 * ts_yyy_yy[i] * gfe_0 + 2.0 * ts_yyy_xyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xyy[i] * gfe_0 + 12.0 * ts_xyy_xy[i] * gfe_0 + 6.0 * ts_xyy_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_x[i] * gfe_0 + 4.0 * ts_xyyy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_xyy[i] * gfe_0 + ts_xyyy_xyy[i] * rgc2_0;

        gr_xyyy_xyz[i] = 2.0 * ts_yyy_yz[i] * gfe_0 + 2.0 * ts_yyy_xyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xyz[i] * gfe_0 + 6.0 * ts_xyy_xz[i] * gfe_0 + 6.0 * ts_xyy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xyz[i] * gfe_0 + ts_xyyy_xyz[i] * rgc2_0;

        gr_xyyy_xzz[i] = 2.0 * ts_yyy_zz[i] * gfe_0 + 2.0 * ts_yyy_xzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xzz[i] * gfe_0 + 6.0 * ts_xyy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_x[i] * gfe_0 + 4.0 * ts_xyyy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xzz[i] * gfe_0 + ts_xyyy_xzz[i] * rgc2_0;

        gr_xyyy_yyy[i] = 2.0 * ts_yyy_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yyy[i] * gfe_0 + 18.0 * ts_xyy_yy[i] * gfe_0 + 6.0 * ts_xyy_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyy_y[i] * gfe_0 + 6.0 * ts_xyyy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_yyy[i] * gfe_0 + ts_xyyy_yyy[i] * rgc2_0;

        gr_xyyy_yyz[i] = 2.0 * ts_yyy_yyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yyz[i] * gfe_0 + 12.0 * ts_xyy_yz[i] * gfe_0 + 6.0 * ts_xyy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_z[i] * gfe_0 + 4.0 * ts_xyyy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_yyz[i] * gfe_0 + ts_xyyy_yyz[i] * rgc2_0;

        gr_xyyy_yzz[i] = 2.0 * ts_yyy_yzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yzz[i] * gfe_0 + 6.0 * ts_xyy_zz[i] * gfe_0 + 6.0 * ts_xyy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_y[i] * gfe_0 + 4.0 * ts_xyyy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_yzz[i] * gfe_0 + ts_xyyy_yzz[i] * rgc2_0;

        gr_xyyy_zzz[i] = 2.0 * ts_yyy_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_zzz[i] * gfe_0 + 6.0 * ts_xyy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyyy_z[i] * gfe_0 + 6.0 * ts_xyyy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_zzz[i] * gfe_0 + ts_xyyy_zzz[i] * rgc2_0;
    }

    // Set up 70-80 components of targeted buffer : GF

    auto gr_xyyz_xxx = pbuffer.data(idx_g_gf + 70);

    auto gr_xyyz_xxy = pbuffer.data(idx_g_gf + 71);

    auto gr_xyyz_xxz = pbuffer.data(idx_g_gf + 72);

    auto gr_xyyz_xyy = pbuffer.data(idx_g_gf + 73);

    auto gr_xyyz_xyz = pbuffer.data(idx_g_gf + 74);

    auto gr_xyyz_xzz = pbuffer.data(idx_g_gf + 75);

    auto gr_xyyz_yyy = pbuffer.data(idx_g_gf + 76);

    auto gr_xyyz_yyz = pbuffer.data(idx_g_gf + 77);

    auto gr_xyyz_yzz = pbuffer.data(idx_g_gf + 78);

    auto gr_xyyz_zzz = pbuffer.data(idx_g_gf + 79);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_xxx, gr_xyyz_xxy, gr_xyyz_xxz, gr_xyyz_xyy, gr_xyyz_xyz, gr_xyyz_xzz, gr_xyyz_yyy, gr_xyyz_yyz, gr_xyyz_yzz, gr_xyyz_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xxx, ts_xyyz_xxy, ts_xyyz_xxz, ts_xyyz_xy, ts_xyyz_xyy, ts_xyyz_xyz, ts_xyyz_xz, ts_xyyz_xzz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yyy, ts_xyyz_yyz, ts_xyyz_yz, ts_xyyz_yzz, ts_xyyz_z, ts_xyyz_zz, ts_xyyz_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyyz_xxx[i] = 6.0 * ts_yyz_xx[i] * gfe_0 + 2.0 * ts_yyz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxx[i] * gfe_0 + 4.0 * ts_xyz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyyz_x[i] * gfe_0 + 6.0 * ts_xyyz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyz_xxx[i] * gfe_0 + ts_xyyz_xxx[i] * rgc2_0;

        gr_xyyz_xxy[i] = 4.0 * ts_yyz_xy[i] * gfe_0 + 2.0 * ts_yyz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxy[i] * gfe_0 + 4.0 * ts_xyz_xx[i] * gfe_0 + 4.0 * ts_xyz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_y[i] * gfe_0 + 4.0 * ts_xyyz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_xxy[i] * gfe_0 + ts_xyyz_xxy[i] * rgc2_0;

        gr_xyyz_xxz[i] = 4.0 * ts_yyz_xz[i] * gfe_0 + 2.0 * ts_yyz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xxz[i] * gfe_0 + 4.0 * ts_xyz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xx[i] * gfe_0 + 2.0 * ts_xyy_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_z[i] * gfe_0 + 4.0 * ts_xyyz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xxz[i] * gfe_0 + ts_xyyz_xxz[i] * rgc2_0;

        gr_xyyz_xyy[i] = 2.0 * ts_yyz_yy[i] * gfe_0 + 2.0 * ts_yyz_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xyy[i] * gfe_0 + 8.0 * ts_xyz_xy[i] * gfe_0 + 4.0 * ts_xyz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_x[i] * gfe_0 + 4.0 * ts_xyyz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_xyy[i] * gfe_0 + ts_xyyz_xyy[i] * rgc2_0;

        gr_xyyz_xyz[i] = 2.0 * ts_yyz_yz[i] * gfe_0 + 2.0 * ts_yyz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xyz[i] * gfe_0 + 4.0 * ts_xyz_xz[i] * gfe_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xy[i] * gfe_0 + 2.0 * ts_xyy_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xyz[i] * gfe_0 + ts_xyyz_xyz[i] * rgc2_0;

        gr_xyyz_xzz[i] = 2.0 * ts_yyz_zz[i] * gfe_0 + 2.0 * ts_yyz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xzz[i] * gfe_0 + 4.0 * ts_xyz_xzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xyy_xz[i] * gfe_0 + 2.0 * ts_xyy_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_x[i] * gfe_0 + 4.0 * ts_xyyz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xzz[i] * gfe_0 + ts_xyyz_xzz[i] * rgc2_0;

        gr_xyyz_yyy[i] = 2.0 * ts_yyz_yyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yyy[i] * gfe_0 + 12.0 * ts_xyz_yy[i] * gfe_0 + 4.0 * ts_xyz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyyz_y[i] * gfe_0 + 6.0 * ts_xyyz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_yyy[i] * gfe_0 + ts_xyyz_yyy[i] * rgc2_0;

        gr_xyyz_yyz[i] = 2.0 * ts_yyz_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yyz[i] * gfe_0 + 8.0 * ts_xyz_yz[i] * gfe_0 + 4.0 * ts_xyz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yy[i] * gfe_0 + 2.0 * ts_xyy_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_z[i] * gfe_0 + 4.0 * ts_xyyz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_yyz[i] * gfe_0 + ts_xyyz_yyz[i] * rgc2_0;

        gr_xyyz_yzz[i] = 2.0 * ts_yyz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yzz[i] * gfe_0 + 4.0 * ts_xyz_zz[i] * gfe_0 + 4.0 * ts_xyz_yzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xyy_yz[i] * gfe_0 + 2.0 * ts_xyy_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_y[i] * gfe_0 + 4.0 * ts_xyyz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_yzz[i] * gfe_0 + ts_xyyz_yzz[i] * rgc2_0;

        gr_xyyz_zzz[i] = 2.0 * ts_yyz_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zzz[i] * gfe_0 + 4.0 * ts_xyz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_xyy_zz[i] * gfe_0 + 2.0 * ts_xyy_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyyz_z[i] * gfe_0 + 6.0 * ts_xyyz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_zzz[i] * gfe_0 + ts_xyyz_zzz[i] * rgc2_0;
    }

    // Set up 80-90 components of targeted buffer : GF

    auto gr_xyzz_xxx = pbuffer.data(idx_g_gf + 80);

    auto gr_xyzz_xxy = pbuffer.data(idx_g_gf + 81);

    auto gr_xyzz_xxz = pbuffer.data(idx_g_gf + 82);

    auto gr_xyzz_xyy = pbuffer.data(idx_g_gf + 83);

    auto gr_xyzz_xyz = pbuffer.data(idx_g_gf + 84);

    auto gr_xyzz_xzz = pbuffer.data(idx_g_gf + 85);

    auto gr_xyzz_yyy = pbuffer.data(idx_g_gf + 86);

    auto gr_xyzz_yyz = pbuffer.data(idx_g_gf + 87);

    auto gr_xyzz_yzz = pbuffer.data(idx_g_gf + 88);

    auto gr_xyzz_zzz = pbuffer.data(idx_g_gf + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_xxx, gr_xyzz_xxy, gr_xyzz_xxz, gr_xyzz_xyy, gr_xyzz_xyz, gr_xyzz_xzz, gr_xyzz_yyy, gr_xyzz_yyz, gr_xyzz_yzz, gr_xyzz_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xxx, ts_xyzz_xxy, ts_xyzz_xxz, ts_xyzz_xy, ts_xyzz_xyy, ts_xyzz_xyz, ts_xyzz_xz, ts_xyzz_xzz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yz, ts_xyzz_yzz, ts_xyzz_z, ts_xyzz_zz, ts_xyzz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyzz_xxx[i] = 6.0 * ts_yzz_xx[i] * gfe_0 + 2.0 * ts_yzz_xxx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxx[i] * gfe_0 + 4.0 * ts_xyz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyzz_x[i] * gfe_0 + 6.0 * ts_xyzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyzz_xxx[i] * gfe_0 + ts_xyzz_xxx[i] * rgc2_0;

        gr_xyzz_xxy[i] = 4.0 * ts_yzz_xy[i] * gfe_0 + 2.0 * ts_yzz_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xx[i] * gfe_0 + 2.0 * ts_xzz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxy[i] * gfe_0 + 4.0 * ts_xyz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_y[i] * gfe_0 + 4.0 * ts_xyzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_xxy[i] * gfe_0 + ts_xyzz_xxy[i] * rgc2_0;

        gr_xyzz_xxz[i] = 4.0 * ts_yzz_xz[i] * gfe_0 + 2.0 * ts_yzz_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xxz[i] * gfe_0 + 4.0 * ts_xyz_xx[i] * gfe_0 + 4.0 * ts_xyz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_z[i] * gfe_0 + 4.0 * ts_xyzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xxz[i] * gfe_0 + ts_xyzz_xxz[i] * rgc2_0;

        gr_xyzz_xyy[i] = 2.0 * ts_yzz_yy[i] * gfe_0 + 2.0 * ts_yzz_xyy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xzz_xy[i] * gfe_0 + 2.0 * ts_xzz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyy[i] * gfe_0 + 4.0 * ts_xyz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_x[i] * gfe_0 + 4.0 * ts_xyzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_xyy[i] * gfe_0 + ts_xyzz_xyy[i] * rgc2_0;

        gr_xyzz_xyz[i] = 2.0 * ts_yzz_yz[i] * gfe_0 + 2.0 * ts_yzz_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xz[i] * gfe_0 + 2.0 * ts_xzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xyz[i] * gfe_0 + 4.0 * ts_xyz_xy[i] * gfe_0 + 4.0 * ts_xyz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xyz[i] * gfe_0 + ts_xyzz_xyz[i] * rgc2_0;

        gr_xyzz_xzz[i] = 2.0 * ts_yzz_zz[i] * gfe_0 + 2.0 * ts_yzz_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xzz[i] * gfe_0 + 8.0 * ts_xyz_xz[i] * gfe_0 + 4.0 * ts_xyz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_x[i] * gfe_0 + 4.0 * ts_xyzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xzz[i] * gfe_0 + ts_xyzz_xzz[i] * rgc2_0;

        gr_xyzz_yyy[i] = 2.0 * ts_yzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xzz_yy[i] * gfe_0 + 2.0 * ts_xzz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyy[i] * gfe_0 + 4.0 * ts_xyz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyzz_y[i] * gfe_0 + 6.0 * ts_xyzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_yyy[i] * gfe_0 + ts_xyzz_yyy[i] * rgc2_0;

        gr_xyzz_yyz[i] = 2.0 * ts_yzz_yyz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xzz_yz[i] * gfe_0 + 2.0 * ts_xzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yyz[i] * gfe_0 + 4.0 * ts_xyz_yy[i] * gfe_0 + 4.0 * ts_xyz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_z[i] * gfe_0 + 4.0 * ts_xyzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_yyz[i] * gfe_0 + ts_xyzz_yyz[i] * rgc2_0;

        gr_xyzz_yzz[i] = 2.0 * ts_yzz_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_zz[i] * gfe_0 + 2.0 * ts_xzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yzz[i] * gfe_0 + 8.0 * ts_xyz_yz[i] * gfe_0 + 4.0 * ts_xyz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_y[i] * gfe_0 + 4.0 * ts_xyzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_yzz[i] * gfe_0 + ts_xyzz_yzz[i] * rgc2_0;

        gr_xyzz_zzz[i] = 2.0 * ts_yzz_zzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_zzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zzz[i] * gfe_0 + 12.0 * ts_xyz_zz[i] * gfe_0 + 4.0 * ts_xyz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xyzz_z[i] * gfe_0 + 6.0 * ts_xyzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_zzz[i] * gfe_0 + ts_xyzz_zzz[i] * rgc2_0;
    }

    // Set up 90-100 components of targeted buffer : GF

    auto gr_xzzz_xxx = pbuffer.data(idx_g_gf + 90);

    auto gr_xzzz_xxy = pbuffer.data(idx_g_gf + 91);

    auto gr_xzzz_xxz = pbuffer.data(idx_g_gf + 92);

    auto gr_xzzz_xyy = pbuffer.data(idx_g_gf + 93);

    auto gr_xzzz_xyz = pbuffer.data(idx_g_gf + 94);

    auto gr_xzzz_xzz = pbuffer.data(idx_g_gf + 95);

    auto gr_xzzz_yyy = pbuffer.data(idx_g_gf + 96);

    auto gr_xzzz_yyz = pbuffer.data(idx_g_gf + 97);

    auto gr_xzzz_yzz = pbuffer.data(idx_g_gf + 98);

    auto gr_xzzz_zzz = pbuffer.data(idx_g_gf + 99);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_xxx, gr_xzzz_xxy, gr_xzzz_xxz, gr_xzzz_xyy, gr_xzzz_xyz, gr_xzzz_xzz, gr_xzzz_yyy, gr_xzzz_yyz, gr_xzzz_yzz, gr_xzzz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xxx, ts_xzzz_xxy, ts_xzzz_xxz, ts_xzzz_xy, ts_xzzz_xyy, ts_xzzz_xyz, ts_xzzz_xz, ts_xzzz_xzz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yz, ts_xzzz_yzz, ts_xzzz_z, ts_xzzz_zz, ts_xzzz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xzzz_xxx[i] = 6.0 * ts_zzz_xx[i] * gfe_0 + 2.0 * ts_zzz_xxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxx[i] * gfe_0 + 6.0 * ts_xzz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzzz_x[i] * gfe_0 + 6.0 * ts_xzzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzzz_xxx[i] * gfe_0 + ts_xzzz_xxx[i] * rgc2_0;

        gr_xzzz_xxy[i] = 4.0 * ts_zzz_xy[i] * gfe_0 + 2.0 * ts_zzz_xxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxy[i] * gfe_0 + 6.0 * ts_xzz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_y[i] * gfe_0 + 4.0 * ts_xzzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_xxy[i] * gfe_0 + ts_xzzz_xxy[i] * rgc2_0;

        gr_xzzz_xxz[i] = 4.0 * ts_zzz_xz[i] * gfe_0 + 2.0 * ts_zzz_xxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xxz[i] * gfe_0 + 6.0 * ts_xzz_xx[i] * gfe_0 + 6.0 * ts_xzz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_z[i] * gfe_0 + 4.0 * ts_xzzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xxz[i] * gfe_0 + ts_xzzz_xxz[i] * rgc2_0;

        gr_xzzz_xyy[i] = 2.0 * ts_zzz_yy[i] * gfe_0 + 2.0 * ts_zzz_xyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xyy[i] * gfe_0 + 6.0 * ts_xzz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_x[i] * gfe_0 + 4.0 * ts_xzzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_xyy[i] * gfe_0 + ts_xzzz_xyy[i] * rgc2_0;

        gr_xzzz_xyz[i] = 2.0 * ts_zzz_yz[i] * gfe_0 + 2.0 * ts_zzz_xyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xyz[i] * gfe_0 + 6.0 * ts_xzz_xy[i] * gfe_0 + 6.0 * ts_xzz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xyz[i] * gfe_0 + ts_xzzz_xyz[i] * rgc2_0;

        gr_xzzz_xzz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 + 2.0 * ts_zzz_xzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xzz[i] * gfe_0 + 12.0 * ts_xzz_xz[i] * gfe_0 + 6.0 * ts_xzz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_x[i] * gfe_0 + 4.0 * ts_xzzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xzz[i] * gfe_0 + ts_xzzz_xzz[i] * rgc2_0;

        gr_xzzz_yyy[i] = 2.0 * ts_zzz_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yyy[i] * gfe_0 + 6.0 * ts_xzz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzzz_y[i] * gfe_0 + 6.0 * ts_xzzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_yyy[i] * gfe_0 + ts_xzzz_yyy[i] * rgc2_0;

        gr_xzzz_yyz[i] = 2.0 * ts_zzz_yyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yyz[i] * gfe_0 + 6.0 * ts_xzz_yy[i] * gfe_0 + 6.0 * ts_xzz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_z[i] * gfe_0 + 4.0 * ts_xzzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_yyz[i] * gfe_0 + ts_xzzz_yyz[i] * rgc2_0;

        gr_xzzz_yzz[i] = 2.0 * ts_zzz_yzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yzz[i] * gfe_0 + 12.0 * ts_xzz_yz[i] * gfe_0 + 6.0 * ts_xzz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_y[i] * gfe_0 + 4.0 * ts_xzzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_yzz[i] * gfe_0 + ts_xzzz_yzz[i] * rgc2_0;

        gr_xzzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_zzz[i] * gfe_0 + 18.0 * ts_xzz_zz[i] * gfe_0 + 6.0 * ts_xzz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_xzzz_z[i] * gfe_0 + 6.0 * ts_xzzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_zzz[i] * gfe_0 + ts_xzzz_zzz[i] * rgc2_0;
    }

    // Set up 100-110 components of targeted buffer : GF

    auto gr_yyyy_xxx = pbuffer.data(idx_g_gf + 100);

    auto gr_yyyy_xxy = pbuffer.data(idx_g_gf + 101);

    auto gr_yyyy_xxz = pbuffer.data(idx_g_gf + 102);

    auto gr_yyyy_xyy = pbuffer.data(idx_g_gf + 103);

    auto gr_yyyy_xyz = pbuffer.data(idx_g_gf + 104);

    auto gr_yyyy_xzz = pbuffer.data(idx_g_gf + 105);

    auto gr_yyyy_yyy = pbuffer.data(idx_g_gf + 106);

    auto gr_yyyy_yyz = pbuffer.data(idx_g_gf + 107);

    auto gr_yyyy_yzz = pbuffer.data(idx_g_gf + 108);

    auto gr_yyyy_zzz = pbuffer.data(idx_g_gf + 109);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_xxx, gr_yyyy_xxy, gr_yyyy_xxz, gr_yyyy_xyy, gr_yyyy_xyz, gr_yyyy_xzz, gr_yyyy_yyy, gr_yyyy_yyz, gr_yyyy_yzz, gr_yyyy_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xy, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xz, ts_yyyy_xzz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yz, ts_yyyy_yzz, ts_yyyy_z, ts_yyyy_zz, ts_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyyy_xxx[i] = 12.0 * ts_yy_xxx[i] * gfe_0 + 8.0 * ts_yyy_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyy_x[i] * gfe_0 + 6.0 * ts_yyyy_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyy_xxx[i] * gfe_0 + ts_yyyy_xxx[i] * rgc2_0;

        gr_yyyy_xxy[i] = 12.0 * ts_yy_xxy[i] * gfe_0 + 8.0 * ts_yyy_xx[i] * gfe_0 + 8.0 * ts_yyy_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_y[i] * gfe_0 + 4.0 * ts_yyyy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_xxy[i] * gfe_0 + ts_yyyy_xxy[i] * rgc2_0;

        gr_yyyy_xxz[i] = 12.0 * ts_yy_xxz[i] * gfe_0 + 8.0 * ts_yyy_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_z[i] * gfe_0 + 4.0 * ts_yyyy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xxz[i] * gfe_0 + ts_yyyy_xxz[i] * rgc2_0;

        gr_yyyy_xyy[i] = 12.0 * ts_yy_xyy[i] * gfe_0 + 16.0 * ts_yyy_xy[i] * gfe_0 + 8.0 * ts_yyy_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_x[i] * gfe_0 + 4.0 * ts_yyyy_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_xyy[i] * gfe_0 + ts_yyyy_xyy[i] * rgc2_0;

        gr_yyyy_xyz[i] = 12.0 * ts_yy_xyz[i] * gfe_0 + 8.0 * ts_yyy_xz[i] * gfe_0 + 8.0 * ts_yyy_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xyz[i] * gfe_0 + ts_yyyy_xyz[i] * rgc2_0;

        gr_yyyy_xzz[i] = 12.0 * ts_yy_xzz[i] * gfe_0 + 8.0 * ts_yyy_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_x[i] * gfe_0 + 4.0 * ts_yyyy_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xzz[i] * gfe_0 + ts_yyyy_xzz[i] * rgc2_0;

        gr_yyyy_yyy[i] = 12.0 * ts_yy_yyy[i] * gfe_0 + 24.0 * ts_yyy_yy[i] * gfe_0 + 8.0 * ts_yyy_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyy_y[i] * gfe_0 + 6.0 * ts_yyyy_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_yyy[i] * gfe_0 + ts_yyyy_yyy[i] * rgc2_0;

        gr_yyyy_yyz[i] = 12.0 * ts_yy_yyz[i] * gfe_0 + 16.0 * ts_yyy_yz[i] * gfe_0 + 8.0 * ts_yyy_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_z[i] * gfe_0 + 4.0 * ts_yyyy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_yyz[i] * gfe_0 + ts_yyyy_yyz[i] * rgc2_0;

        gr_yyyy_yzz[i] = 12.0 * ts_yy_yzz[i] * gfe_0 + 8.0 * ts_yyy_zz[i] * gfe_0 + 8.0 * ts_yyy_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_y[i] * gfe_0 + 4.0 * ts_yyyy_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_yzz[i] * gfe_0 + ts_yyyy_yzz[i] * rgc2_0;

        gr_yyyy_zzz[i] = 12.0 * ts_yy_zzz[i] * gfe_0 + 8.0 * ts_yyy_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyyy_z[i] * gfe_0 + 6.0 * ts_yyyy_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_zzz[i] * gfe_0 + ts_yyyy_zzz[i] * rgc2_0;
    }

    // Set up 110-120 components of targeted buffer : GF

    auto gr_yyyz_xxx = pbuffer.data(idx_g_gf + 110);

    auto gr_yyyz_xxy = pbuffer.data(idx_g_gf + 111);

    auto gr_yyyz_xxz = pbuffer.data(idx_g_gf + 112);

    auto gr_yyyz_xyy = pbuffer.data(idx_g_gf + 113);

    auto gr_yyyz_xyz = pbuffer.data(idx_g_gf + 114);

    auto gr_yyyz_xzz = pbuffer.data(idx_g_gf + 115);

    auto gr_yyyz_yyy = pbuffer.data(idx_g_gf + 116);

    auto gr_yyyz_yyz = pbuffer.data(idx_g_gf + 117);

    auto gr_yyyz_yzz = pbuffer.data(idx_g_gf + 118);

    auto gr_yyyz_zzz = pbuffer.data(idx_g_gf + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_xxx, gr_yyyz_xxy, gr_yyyz_xxz, gr_yyyz_xyy, gr_yyyz_xyz, gr_yyyz_xzz, gr_yyyz_yyy, gr_yyyz_yyz, gr_yyyz_yzz, gr_yyyz_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xxx, ts_yyyz_xxy, ts_yyyz_xxz, ts_yyyz_xy, ts_yyyz_xyy, ts_yyyz_xyz, ts_yyyz_xz, ts_yyyz_xzz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yyy, ts_yyyz_yyz, ts_yyyz_yz, ts_yyyz_yzz, ts_yyyz_z, ts_yyyz_zz, ts_yyyz_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyyz_xxx[i] = 6.0 * ts_yz_xxx[i] * gfe_0 + 6.0 * ts_yyz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyyz_x[i] * gfe_0 + 6.0 * ts_yyyz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyz_xxx[i] * gfe_0 + ts_yyyz_xxx[i] * rgc2_0;

        gr_yyyz_xxy[i] = 6.0 * ts_yz_xxy[i] * gfe_0 + 6.0 * ts_yyz_xx[i] * gfe_0 + 6.0 * ts_yyz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_y[i] * gfe_0 + 4.0 * ts_yyyz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_xxy[i] * gfe_0 + ts_yyyz_xxy[i] * rgc2_0;

        gr_yyyz_xxz[i] = 6.0 * ts_yz_xxz[i] * gfe_0 + 6.0 * ts_yyz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xx[i] * gfe_0 + 2.0 * ts_yyy_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_z[i] * gfe_0 + 4.0 * ts_yyyz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xxz[i] * gfe_0 + ts_yyyz_xxz[i] * rgc2_0;

        gr_yyyz_xyy[i] = 6.0 * ts_yz_xyy[i] * gfe_0 + 12.0 * ts_yyz_xy[i] * gfe_0 + 6.0 * ts_yyz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_x[i] * gfe_0 + 4.0 * ts_yyyz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_xyy[i] * gfe_0 + ts_yyyz_xyy[i] * rgc2_0;

        gr_yyyz_xyz[i] = 6.0 * ts_yz_xyz[i] * gfe_0 + 6.0 * ts_yyz_xz[i] * gfe_0 + 6.0 * ts_yyz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xy[i] * gfe_0 + 2.0 * ts_yyy_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xyz[i] * gfe_0 + ts_yyyz_xyz[i] * rgc2_0;

        gr_yyyz_xzz[i] = 6.0 * ts_yz_xzz[i] * gfe_0 + 6.0 * ts_yyz_xzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yyy_xz[i] * gfe_0 + 2.0 * ts_yyy_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_x[i] * gfe_0 + 4.0 * ts_yyyz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xzz[i] * gfe_0 + ts_yyyz_xzz[i] * rgc2_0;

        gr_yyyz_yyy[i] = 6.0 * ts_yz_yyy[i] * gfe_0 + 18.0 * ts_yyz_yy[i] * gfe_0 + 6.0 * ts_yyz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyyz_y[i] * gfe_0 + 6.0 * ts_yyyz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_yyy[i] * gfe_0 + ts_yyyz_yyy[i] * rgc2_0;

        gr_yyyz_yyz[i] = 6.0 * ts_yz_yyz[i] * gfe_0 + 12.0 * ts_yyz_yz[i] * gfe_0 + 6.0 * ts_yyz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yy[i] * gfe_0 + 2.0 * ts_yyy_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_z[i] * gfe_0 + 4.0 * ts_yyyz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_yyz[i] * gfe_0 + ts_yyyz_yyz[i] * rgc2_0;

        gr_yyyz_yzz[i] = 6.0 * ts_yz_yzz[i] * gfe_0 + 6.0 * ts_yyz_zz[i] * gfe_0 + 6.0 * ts_yyz_yzz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yyy_yz[i] * gfe_0 + 2.0 * ts_yyy_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_y[i] * gfe_0 + 4.0 * ts_yyyz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_yzz[i] * gfe_0 + ts_yyyz_yzz[i] * rgc2_0;

        gr_yyyz_zzz[i] = 6.0 * ts_yz_zzz[i] * gfe_0 + 6.0 * ts_yyz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yyy_zz[i] * gfe_0 + 2.0 * ts_yyy_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyyz_z[i] * gfe_0 + 6.0 * ts_yyyz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_zzz[i] * gfe_0 + ts_yyyz_zzz[i] * rgc2_0;
    }

    // Set up 120-130 components of targeted buffer : GF

    auto gr_yyzz_xxx = pbuffer.data(idx_g_gf + 120);

    auto gr_yyzz_xxy = pbuffer.data(idx_g_gf + 121);

    auto gr_yyzz_xxz = pbuffer.data(idx_g_gf + 122);

    auto gr_yyzz_xyy = pbuffer.data(idx_g_gf + 123);

    auto gr_yyzz_xyz = pbuffer.data(idx_g_gf + 124);

    auto gr_yyzz_xzz = pbuffer.data(idx_g_gf + 125);

    auto gr_yyzz_yyy = pbuffer.data(idx_g_gf + 126);

    auto gr_yyzz_yyz = pbuffer.data(idx_g_gf + 127);

    auto gr_yyzz_yzz = pbuffer.data(idx_g_gf + 128);

    auto gr_yyzz_zzz = pbuffer.data(idx_g_gf + 129);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_xxx, gr_yyzz_xxy, gr_yyzz_xxz, gr_yyzz_xyy, gr_yyzz_xyz, gr_yyzz_xzz, gr_yyzz_yyy, gr_yyzz_yyz, gr_yyzz_yzz, gr_yyzz_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xxx, ts_yyzz_xxy, ts_yyzz_xxz, ts_yyzz_xy, ts_yyzz_xyy, ts_yyzz_xyz, ts_yyzz_xz, ts_yyzz_xzz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yz, ts_yyzz_yzz, ts_yyzz_z, ts_yyzz_zz, ts_yyzz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyzz_xxx[i] = 2.0 * ts_zz_xxx[i] * gfe_0 + 4.0 * ts_yzz_xxx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxx[i] * gfe_0 + 4.0 * ts_yyz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyzz_x[i] * gfe_0 + 6.0 * ts_yyzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyzz_xxx[i] * gfe_0 + ts_yyzz_xxx[i] * rgc2_0;

        gr_yyzz_xxy[i] = 2.0 * ts_zz_xxy[i] * gfe_0 + 4.0 * ts_yzz_xx[i] * gfe_0 + 4.0 * ts_yzz_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxy[i] * gfe_0 + 4.0 * ts_yyz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_y[i] * gfe_0 + 4.0 * ts_yyzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_xxy[i] * gfe_0 + ts_yyzz_xxy[i] * rgc2_0;

        gr_yyzz_xxz[i] = 2.0 * ts_zz_xxz[i] * gfe_0 + 4.0 * ts_yzz_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xxz[i] * gfe_0 + 4.0 * ts_yyz_xx[i] * gfe_0 + 4.0 * ts_yyz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_z[i] * gfe_0 + 4.0 * ts_yyzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xxz[i] * gfe_0 + ts_yyzz_xxz[i] * rgc2_0;

        gr_yyzz_xyy[i] = 2.0 * ts_zz_xyy[i] * gfe_0 + 8.0 * ts_yzz_xy[i] * gfe_0 + 4.0 * ts_yzz_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyy[i] * gfe_0 + 4.0 * ts_yyz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_x[i] * gfe_0 + 4.0 * ts_yyzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_xyy[i] * gfe_0 + ts_yyzz_xyy[i] * rgc2_0;

        gr_yyzz_xyz[i] = 2.0 * ts_zz_xyz[i] * gfe_0 + 4.0 * ts_yzz_xz[i] * gfe_0 + 4.0 * ts_yzz_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xyz[i] * gfe_0 + 4.0 * ts_yyz_xy[i] * gfe_0 + 4.0 * ts_yyz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xyz[i] * gfe_0 + ts_yyzz_xyz[i] * rgc2_0;

        gr_yyzz_xzz[i] = 2.0 * ts_zz_xzz[i] * gfe_0 + 4.0 * ts_yzz_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xzz[i] * gfe_0 + 8.0 * ts_yyz_xz[i] * gfe_0 + 4.0 * ts_yyz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_x[i] * gfe_0 + 4.0 * ts_yyzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xzz[i] * gfe_0 + ts_yyzz_xzz[i] * rgc2_0;

        gr_yyzz_yyy[i] = 2.0 * ts_zz_yyy[i] * gfe_0 + 12.0 * ts_yzz_yy[i] * gfe_0 + 4.0 * ts_yzz_yyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyy[i] * gfe_0 + 4.0 * ts_yyz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyzz_y[i] * gfe_0 + 6.0 * ts_yyzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_yyy[i] * gfe_0 + ts_yyzz_yyy[i] * rgc2_0;

        gr_yyzz_yyz[i] = 2.0 * ts_zz_yyz[i] * gfe_0 + 8.0 * ts_yzz_yz[i] * gfe_0 + 4.0 * ts_yzz_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yyz[i] * gfe_0 + 4.0 * ts_yyz_yy[i] * gfe_0 + 4.0 * ts_yyz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_z[i] * gfe_0 + 4.0 * ts_yyzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_yyz[i] * gfe_0 + ts_yyzz_yyz[i] * rgc2_0;

        gr_yyzz_yzz[i] = 2.0 * ts_zz_yzz[i] * gfe_0 + 4.0 * ts_yzz_zz[i] * gfe_0 + 4.0 * ts_yzz_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yzz[i] * gfe_0 + 8.0 * ts_yyz_yz[i] * gfe_0 + 4.0 * ts_yyz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_y[i] * gfe_0 + 4.0 * ts_yyzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_yzz[i] * gfe_0 + ts_yyzz_yzz[i] * rgc2_0;

        gr_yyzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 + 4.0 * ts_yzz_zzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zzz[i] * gfe_0 + 12.0 * ts_yyz_zz[i] * gfe_0 + 4.0 * ts_yyz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yyzz_z[i] * gfe_0 + 6.0 * ts_yyzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_zzz[i] * gfe_0 + ts_yyzz_zzz[i] * rgc2_0;
    }

    // Set up 130-140 components of targeted buffer : GF

    auto gr_yzzz_xxx = pbuffer.data(idx_g_gf + 130);

    auto gr_yzzz_xxy = pbuffer.data(idx_g_gf + 131);

    auto gr_yzzz_xxz = pbuffer.data(idx_g_gf + 132);

    auto gr_yzzz_xyy = pbuffer.data(idx_g_gf + 133);

    auto gr_yzzz_xyz = pbuffer.data(idx_g_gf + 134);

    auto gr_yzzz_xzz = pbuffer.data(idx_g_gf + 135);

    auto gr_yzzz_yyy = pbuffer.data(idx_g_gf + 136);

    auto gr_yzzz_yyz = pbuffer.data(idx_g_gf + 137);

    auto gr_yzzz_yzz = pbuffer.data(idx_g_gf + 138);

    auto gr_yzzz_zzz = pbuffer.data(idx_g_gf + 139);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_xxx, gr_yzzz_xxy, gr_yzzz_xxz, gr_yzzz_xyy, gr_yzzz_xyz, gr_yzzz_xzz, gr_yzzz_yyy, gr_yzzz_yyz, gr_yzzz_yzz, gr_yzzz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xxx, ts_yzzz_xxy, ts_yzzz_xxz, ts_yzzz_xy, ts_yzzz_xyy, ts_yzzz_xyz, ts_yzzz_xz, ts_yzzz_xzz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yz, ts_yzzz_yzz, ts_yzzz_z, ts_yzzz_zz, ts_yzzz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yzzz_xxx[i] = 2.0 * ts_zzz_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxx[i] * gfe_0 + 6.0 * ts_yzz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzzz_x[i] * gfe_0 + 6.0 * ts_yzzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzzz_xxx[i] * gfe_0 + ts_yzzz_xxx[i] * rgc2_0;

        gr_yzzz_xxy[i] = 2.0 * ts_zzz_xx[i] * gfe_0 + 2.0 * ts_zzz_xxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxy[i] * gfe_0 + 6.0 * ts_yzz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_y[i] * gfe_0 + 4.0 * ts_yzzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_xxy[i] * gfe_0 + ts_yzzz_xxy[i] * rgc2_0;

        gr_yzzz_xxz[i] = 2.0 * ts_zzz_xxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xxz[i] * gfe_0 + 6.0 * ts_yzz_xx[i] * gfe_0 + 6.0 * ts_yzz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_z[i] * gfe_0 + 4.0 * ts_yzzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xxz[i] * gfe_0 + ts_yzzz_xxz[i] * rgc2_0;

        gr_yzzz_xyy[i] = 4.0 * ts_zzz_xy[i] * gfe_0 + 2.0 * ts_zzz_xyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xyy[i] * gfe_0 + 6.0 * ts_yzz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_x[i] * gfe_0 + 4.0 * ts_yzzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_xyy[i] * gfe_0 + ts_yzzz_xyy[i] * rgc2_0;

        gr_yzzz_xyz[i] = 2.0 * ts_zzz_xz[i] * gfe_0 + 2.0 * ts_zzz_xyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xyz[i] * gfe_0 + 6.0 * ts_yzz_xy[i] * gfe_0 + 6.0 * ts_yzz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xyz[i] * gfe_0 + ts_yzzz_xyz[i] * rgc2_0;

        gr_yzzz_xzz[i] = 2.0 * ts_zzz_xzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xzz[i] * gfe_0 + 12.0 * ts_yzz_xz[i] * gfe_0 + 6.0 * ts_yzz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_x[i] * gfe_0 + 4.0 * ts_yzzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xzz[i] * gfe_0 + ts_yzzz_xzz[i] * rgc2_0;

        gr_yzzz_yyy[i] = 6.0 * ts_zzz_yy[i] * gfe_0 + 2.0 * ts_zzz_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yyy[i] * gfe_0 + 6.0 * ts_yzz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzzz_y[i] * gfe_0 + 6.0 * ts_yzzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_yyy[i] * gfe_0 + ts_yzzz_yyy[i] * rgc2_0;

        gr_yzzz_yyz[i] = 4.0 * ts_zzz_yz[i] * gfe_0 + 2.0 * ts_zzz_yyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yyz[i] * gfe_0 + 6.0 * ts_yzz_yy[i] * gfe_0 + 6.0 * ts_yzz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_z[i] * gfe_0 + 4.0 * ts_yzzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_yyz[i] * gfe_0 + ts_yzzz_yyz[i] * rgc2_0;

        gr_yzzz_yzz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 + 2.0 * ts_zzz_yzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yzz[i] * gfe_0 + 12.0 * ts_yzz_yz[i] * gfe_0 + 6.0 * ts_yzz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_y[i] * gfe_0 + 4.0 * ts_yzzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_yzz[i] * gfe_0 + ts_yzzz_yzz[i] * rgc2_0;

        gr_yzzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_zzz[i] * gfe_0 + 18.0 * ts_yzz_zz[i] * gfe_0 + 6.0 * ts_yzz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_yzzz_z[i] * gfe_0 + 6.0 * ts_yzzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_zzz[i] * gfe_0 + ts_yzzz_zzz[i] * rgc2_0;
    }

    // Set up 140-150 components of targeted buffer : GF

    auto gr_zzzz_xxx = pbuffer.data(idx_g_gf + 140);

    auto gr_zzzz_xxy = pbuffer.data(idx_g_gf + 141);

    auto gr_zzzz_xxz = pbuffer.data(idx_g_gf + 142);

    auto gr_zzzz_xyy = pbuffer.data(idx_g_gf + 143);

    auto gr_zzzz_xyz = pbuffer.data(idx_g_gf + 144);

    auto gr_zzzz_xzz = pbuffer.data(idx_g_gf + 145);

    auto gr_zzzz_yyy = pbuffer.data(idx_g_gf + 146);

    auto gr_zzzz_yyz = pbuffer.data(idx_g_gf + 147);

    auto gr_zzzz_yzz = pbuffer.data(idx_g_gf + 148);

    auto gr_zzzz_zzz = pbuffer.data(idx_g_gf + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_xxx, gr_zzzz_xxy, gr_zzzz_xxz, gr_zzzz_xyy, gr_zzzz_xyz, gr_zzzz_xzz, gr_zzzz_yyy, gr_zzzz_yyz, gr_zzzz_yzz, gr_zzzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xy, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xz, ts_zzzz_xzz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yz, ts_zzzz_yzz, ts_zzzz_z, ts_zzzz_zz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zzzz_xxx[i] = 12.0 * ts_zz_xxx[i] * gfe_0 + 8.0 * ts_zzz_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzzz_x[i] * gfe_0 + 6.0 * ts_zzzz_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzzz_xxx[i] * gfe_0 + ts_zzzz_xxx[i] * rgc2_0;

        gr_zzzz_xxy[i] = 12.0 * ts_zz_xxy[i] * gfe_0 + 8.0 * ts_zzz_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_y[i] * gfe_0 + 4.0 * ts_zzzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_xxy[i] * gfe_0 + ts_zzzz_xxy[i] * rgc2_0;

        gr_zzzz_xxz[i] = 12.0 * ts_zz_xxz[i] * gfe_0 + 8.0 * ts_zzz_xx[i] * gfe_0 + 8.0 * ts_zzz_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_z[i] * gfe_0 + 4.0 * ts_zzzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xxz[i] * gfe_0 + ts_zzzz_xxz[i] * rgc2_0;

        gr_zzzz_xyy[i] = 12.0 * ts_zz_xyy[i] * gfe_0 + 8.0 * ts_zzz_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_x[i] * gfe_0 + 4.0 * ts_zzzz_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_xyy[i] * gfe_0 + ts_zzzz_xyy[i] * rgc2_0;

        gr_zzzz_xyz[i] = 12.0 * ts_zz_xyz[i] * gfe_0 + 8.0 * ts_zzz_xy[i] * gfe_0 + 8.0 * ts_zzz_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xyz[i] * gfe_0 + ts_zzzz_xyz[i] * rgc2_0;

        gr_zzzz_xzz[i] = 12.0 * ts_zz_xzz[i] * gfe_0 + 16.0 * ts_zzz_xz[i] * gfe_0 + 8.0 * ts_zzz_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_x[i] * gfe_0 + 4.0 * ts_zzzz_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xzz[i] * gfe_0 + ts_zzzz_xzz[i] * rgc2_0;

        gr_zzzz_yyy[i] = 12.0 * ts_zz_yyy[i] * gfe_0 + 8.0 * ts_zzz_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzzz_y[i] * gfe_0 + 6.0 * ts_zzzz_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_yyy[i] * gfe_0 + ts_zzzz_yyy[i] * rgc2_0;

        gr_zzzz_yyz[i] = 12.0 * ts_zz_yyz[i] * gfe_0 + 8.0 * ts_zzz_yy[i] * gfe_0 + 8.0 * ts_zzz_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_z[i] * gfe_0 + 4.0 * ts_zzzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_yyz[i] * gfe_0 + ts_zzzz_yyz[i] * rgc2_0;

        gr_zzzz_yzz[i] = 12.0 * ts_zz_yzz[i] * gfe_0 + 16.0 * ts_zzz_yz[i] * gfe_0 + 8.0 * ts_zzz_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_y[i] * gfe_0 + 4.0 * ts_zzzz_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_yzz[i] * gfe_0 + ts_zzzz_yzz[i] * rgc2_0;

        gr_zzzz_zzz[i] = 12.0 * ts_zz_zzz[i] * gfe_0 + 24.0 * ts_zzz_zz[i] * gfe_0 + 8.0 * ts_zzz_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_zzzz_z[i] * gfe_0 + 6.0 * ts_zzzz_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_zzz[i] * gfe_0 + ts_zzzz_zzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

