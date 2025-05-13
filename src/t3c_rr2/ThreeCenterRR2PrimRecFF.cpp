#include "ThreeCenterRR2PrimRecFF.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_ff(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_ff,
                  const size_t idx_df,
                  const size_t idx_g_df,
                  const size_t idx_fd,
                  const size_t idx_g_fd,
                  const size_t idx_ff,
                  const size_t idx_g_ff,
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

    // Set up components of auxiliary buffer : DF

    auto gr_xx_xxx = pbuffer.data(idx_g_df);

    auto gr_xx_xxy = pbuffer.data(idx_g_df + 1);

    auto gr_xx_xxz = pbuffer.data(idx_g_df + 2);

    auto gr_xx_xyy = pbuffer.data(idx_g_df + 3);

    auto gr_xx_xyz = pbuffer.data(idx_g_df + 4);

    auto gr_xx_xzz = pbuffer.data(idx_g_df + 5);

    auto gr_xx_yyy = pbuffer.data(idx_g_df + 6);

    auto gr_xx_yyz = pbuffer.data(idx_g_df + 7);

    auto gr_xx_yzz = pbuffer.data(idx_g_df + 8);

    auto gr_xx_zzz = pbuffer.data(idx_g_df + 9);

    auto gr_xy_xxx = pbuffer.data(idx_g_df + 10);

    auto gr_xy_xxy = pbuffer.data(idx_g_df + 11);

    auto gr_xy_xxz = pbuffer.data(idx_g_df + 12);

    auto gr_xy_xyy = pbuffer.data(idx_g_df + 13);

    auto gr_xy_xyz = pbuffer.data(idx_g_df + 14);

    auto gr_xy_xzz = pbuffer.data(idx_g_df + 15);

    auto gr_xy_yyy = pbuffer.data(idx_g_df + 16);

    auto gr_xy_yyz = pbuffer.data(idx_g_df + 17);

    auto gr_xy_yzz = pbuffer.data(idx_g_df + 18);

    auto gr_xy_zzz = pbuffer.data(idx_g_df + 19);

    auto gr_xz_xxx = pbuffer.data(idx_g_df + 20);

    auto gr_xz_xxy = pbuffer.data(idx_g_df + 21);

    auto gr_xz_xxz = pbuffer.data(idx_g_df + 22);

    auto gr_xz_xyy = pbuffer.data(idx_g_df + 23);

    auto gr_xz_xyz = pbuffer.data(idx_g_df + 24);

    auto gr_xz_xzz = pbuffer.data(idx_g_df + 25);

    auto gr_xz_yyy = pbuffer.data(idx_g_df + 26);

    auto gr_xz_yyz = pbuffer.data(idx_g_df + 27);

    auto gr_xz_yzz = pbuffer.data(idx_g_df + 28);

    auto gr_xz_zzz = pbuffer.data(idx_g_df + 29);

    auto gr_yy_xxx = pbuffer.data(idx_g_df + 30);

    auto gr_yy_xxy = pbuffer.data(idx_g_df + 31);

    auto gr_yy_xxz = pbuffer.data(idx_g_df + 32);

    auto gr_yy_xyy = pbuffer.data(idx_g_df + 33);

    auto gr_yy_xyz = pbuffer.data(idx_g_df + 34);

    auto gr_yy_xzz = pbuffer.data(idx_g_df + 35);

    auto gr_yy_yyy = pbuffer.data(idx_g_df + 36);

    auto gr_yy_yyz = pbuffer.data(idx_g_df + 37);

    auto gr_yy_yzz = pbuffer.data(idx_g_df + 38);

    auto gr_yy_zzz = pbuffer.data(idx_g_df + 39);

    auto gr_yz_xxx = pbuffer.data(idx_g_df + 40);

    auto gr_yz_xxy = pbuffer.data(idx_g_df + 41);

    auto gr_yz_xxz = pbuffer.data(idx_g_df + 42);

    auto gr_yz_xyy = pbuffer.data(idx_g_df + 43);

    auto gr_yz_xyz = pbuffer.data(idx_g_df + 44);

    auto gr_yz_xzz = pbuffer.data(idx_g_df + 45);

    auto gr_yz_yyy = pbuffer.data(idx_g_df + 46);

    auto gr_yz_yyz = pbuffer.data(idx_g_df + 47);

    auto gr_yz_yzz = pbuffer.data(idx_g_df + 48);

    auto gr_yz_zzz = pbuffer.data(idx_g_df + 49);

    auto gr_zz_xxx = pbuffer.data(idx_g_df + 50);

    auto gr_zz_xxy = pbuffer.data(idx_g_df + 51);

    auto gr_zz_xxz = pbuffer.data(idx_g_df + 52);

    auto gr_zz_xyy = pbuffer.data(idx_g_df + 53);

    auto gr_zz_xyz = pbuffer.data(idx_g_df + 54);

    auto gr_zz_xzz = pbuffer.data(idx_g_df + 55);

    auto gr_zz_yyy = pbuffer.data(idx_g_df + 56);

    auto gr_zz_yyz = pbuffer.data(idx_g_df + 57);

    auto gr_zz_yzz = pbuffer.data(idx_g_df + 58);

    auto gr_zz_zzz = pbuffer.data(idx_g_df + 59);

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

    // Set up components of auxiliary buffer : FD

    auto gr_xxx_xx = pbuffer.data(idx_g_fd);

    auto gr_xxx_xy = pbuffer.data(idx_g_fd + 1);

    auto gr_xxx_xz = pbuffer.data(idx_g_fd + 2);

    auto gr_xxx_yy = pbuffer.data(idx_g_fd + 3);

    auto gr_xxx_yz = pbuffer.data(idx_g_fd + 4);

    auto gr_xxx_zz = pbuffer.data(idx_g_fd + 5);

    auto gr_xxy_xx = pbuffer.data(idx_g_fd + 6);

    auto gr_xxy_xy = pbuffer.data(idx_g_fd + 7);

    auto gr_xxy_xz = pbuffer.data(idx_g_fd + 8);

    auto gr_xxy_yy = pbuffer.data(idx_g_fd + 9);

    auto gr_xxy_yz = pbuffer.data(idx_g_fd + 10);

    auto gr_xxy_zz = pbuffer.data(idx_g_fd + 11);

    auto gr_xxz_xx = pbuffer.data(idx_g_fd + 12);

    auto gr_xxz_xy = pbuffer.data(idx_g_fd + 13);

    auto gr_xxz_xz = pbuffer.data(idx_g_fd + 14);

    auto gr_xxz_yy = pbuffer.data(idx_g_fd + 15);

    auto gr_xxz_yz = pbuffer.data(idx_g_fd + 16);

    auto gr_xxz_zz = pbuffer.data(idx_g_fd + 17);

    auto gr_xyy_xx = pbuffer.data(idx_g_fd + 18);

    auto gr_xyy_xy = pbuffer.data(idx_g_fd + 19);

    auto gr_xyy_xz = pbuffer.data(idx_g_fd + 20);

    auto gr_xyy_yy = pbuffer.data(idx_g_fd + 21);

    auto gr_xyy_yz = pbuffer.data(idx_g_fd + 22);

    auto gr_xyy_zz = pbuffer.data(idx_g_fd + 23);

    auto gr_xyz_xx = pbuffer.data(idx_g_fd + 24);

    auto gr_xyz_xy = pbuffer.data(idx_g_fd + 25);

    auto gr_xyz_xz = pbuffer.data(idx_g_fd + 26);

    auto gr_xyz_yy = pbuffer.data(idx_g_fd + 27);

    auto gr_xyz_yz = pbuffer.data(idx_g_fd + 28);

    auto gr_xyz_zz = pbuffer.data(idx_g_fd + 29);

    auto gr_xzz_xx = pbuffer.data(idx_g_fd + 30);

    auto gr_xzz_xy = pbuffer.data(idx_g_fd + 31);

    auto gr_xzz_xz = pbuffer.data(idx_g_fd + 32);

    auto gr_xzz_yy = pbuffer.data(idx_g_fd + 33);

    auto gr_xzz_yz = pbuffer.data(idx_g_fd + 34);

    auto gr_xzz_zz = pbuffer.data(idx_g_fd + 35);

    auto gr_yyy_xx = pbuffer.data(idx_g_fd + 36);

    auto gr_yyy_xy = pbuffer.data(idx_g_fd + 37);

    auto gr_yyy_xz = pbuffer.data(idx_g_fd + 38);

    auto gr_yyy_yy = pbuffer.data(idx_g_fd + 39);

    auto gr_yyy_yz = pbuffer.data(idx_g_fd + 40);

    auto gr_yyy_zz = pbuffer.data(idx_g_fd + 41);

    auto gr_yyz_xx = pbuffer.data(idx_g_fd + 42);

    auto gr_yyz_xy = pbuffer.data(idx_g_fd + 43);

    auto gr_yyz_xz = pbuffer.data(idx_g_fd + 44);

    auto gr_yyz_yy = pbuffer.data(idx_g_fd + 45);

    auto gr_yyz_yz = pbuffer.data(idx_g_fd + 46);

    auto gr_yyz_zz = pbuffer.data(idx_g_fd + 47);

    auto gr_yzz_xx = pbuffer.data(idx_g_fd + 48);

    auto gr_yzz_xy = pbuffer.data(idx_g_fd + 49);

    auto gr_yzz_xz = pbuffer.data(idx_g_fd + 50);

    auto gr_yzz_yy = pbuffer.data(idx_g_fd + 51);

    auto gr_yzz_yz = pbuffer.data(idx_g_fd + 52);

    auto gr_yzz_zz = pbuffer.data(idx_g_fd + 53);

    auto gr_zzz_xx = pbuffer.data(idx_g_fd + 54);

    auto gr_zzz_xy = pbuffer.data(idx_g_fd + 55);

    auto gr_zzz_xz = pbuffer.data(idx_g_fd + 56);

    auto gr_zzz_yy = pbuffer.data(idx_g_fd + 57);

    auto gr_zzz_yz = pbuffer.data(idx_g_fd + 58);

    auto gr_zzz_zz = pbuffer.data(idx_g_fd + 59);

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

    // Set up components of auxiliary buffer : FF

    auto gr_xxx_xxx = pbuffer.data(idx_g_ff);

    auto gr_xxx_xxy = pbuffer.data(idx_g_ff + 1);

    auto gr_xxx_xxz = pbuffer.data(idx_g_ff + 2);

    auto gr_xxx_xyy = pbuffer.data(idx_g_ff + 3);

    auto gr_xxx_xyz = pbuffer.data(idx_g_ff + 4);

    auto gr_xxx_xzz = pbuffer.data(idx_g_ff + 5);

    auto gr_xxx_yyy = pbuffer.data(idx_g_ff + 6);

    auto gr_xxx_yyz = pbuffer.data(idx_g_ff + 7);

    auto gr_xxx_yzz = pbuffer.data(idx_g_ff + 8);

    auto gr_xxx_zzz = pbuffer.data(idx_g_ff + 9);

    auto gr_xxy_xxx = pbuffer.data(idx_g_ff + 10);

    auto gr_xxy_xxy = pbuffer.data(idx_g_ff + 11);

    auto gr_xxy_xxz = pbuffer.data(idx_g_ff + 12);

    auto gr_xxy_xyy = pbuffer.data(idx_g_ff + 13);

    auto gr_xxy_xyz = pbuffer.data(idx_g_ff + 14);

    auto gr_xxy_xzz = pbuffer.data(idx_g_ff + 15);

    auto gr_xxy_yyy = pbuffer.data(idx_g_ff + 16);

    auto gr_xxy_yyz = pbuffer.data(idx_g_ff + 17);

    auto gr_xxy_yzz = pbuffer.data(idx_g_ff + 18);

    auto gr_xxy_zzz = pbuffer.data(idx_g_ff + 19);

    auto gr_xxz_xxx = pbuffer.data(idx_g_ff + 20);

    auto gr_xxz_xxy = pbuffer.data(idx_g_ff + 21);

    auto gr_xxz_xxz = pbuffer.data(idx_g_ff + 22);

    auto gr_xxz_xyy = pbuffer.data(idx_g_ff + 23);

    auto gr_xxz_xyz = pbuffer.data(idx_g_ff + 24);

    auto gr_xxz_xzz = pbuffer.data(idx_g_ff + 25);

    auto gr_xxz_yyy = pbuffer.data(idx_g_ff + 26);

    auto gr_xxz_yyz = pbuffer.data(idx_g_ff + 27);

    auto gr_xxz_yzz = pbuffer.data(idx_g_ff + 28);

    auto gr_xxz_zzz = pbuffer.data(idx_g_ff + 29);

    auto gr_xyy_xxx = pbuffer.data(idx_g_ff + 30);

    auto gr_xyy_xxy = pbuffer.data(idx_g_ff + 31);

    auto gr_xyy_xxz = pbuffer.data(idx_g_ff + 32);

    auto gr_xyy_xyy = pbuffer.data(idx_g_ff + 33);

    auto gr_xyy_xyz = pbuffer.data(idx_g_ff + 34);

    auto gr_xyy_xzz = pbuffer.data(idx_g_ff + 35);

    auto gr_xyy_yyy = pbuffer.data(idx_g_ff + 36);

    auto gr_xyy_yyz = pbuffer.data(idx_g_ff + 37);

    auto gr_xyy_yzz = pbuffer.data(idx_g_ff + 38);

    auto gr_xyy_zzz = pbuffer.data(idx_g_ff + 39);

    auto gr_xyz_xxx = pbuffer.data(idx_g_ff + 40);

    auto gr_xyz_xxy = pbuffer.data(idx_g_ff + 41);

    auto gr_xyz_xxz = pbuffer.data(idx_g_ff + 42);

    auto gr_xyz_xyy = pbuffer.data(idx_g_ff + 43);

    auto gr_xyz_xyz = pbuffer.data(idx_g_ff + 44);

    auto gr_xyz_xzz = pbuffer.data(idx_g_ff + 45);

    auto gr_xyz_yyy = pbuffer.data(idx_g_ff + 46);

    auto gr_xyz_yyz = pbuffer.data(idx_g_ff + 47);

    auto gr_xyz_yzz = pbuffer.data(idx_g_ff + 48);

    auto gr_xyz_zzz = pbuffer.data(idx_g_ff + 49);

    auto gr_xzz_xxx = pbuffer.data(idx_g_ff + 50);

    auto gr_xzz_xxy = pbuffer.data(idx_g_ff + 51);

    auto gr_xzz_xxz = pbuffer.data(idx_g_ff + 52);

    auto gr_xzz_xyy = pbuffer.data(idx_g_ff + 53);

    auto gr_xzz_xyz = pbuffer.data(idx_g_ff + 54);

    auto gr_xzz_xzz = pbuffer.data(idx_g_ff + 55);

    auto gr_xzz_yyy = pbuffer.data(idx_g_ff + 56);

    auto gr_xzz_yyz = pbuffer.data(idx_g_ff + 57);

    auto gr_xzz_yzz = pbuffer.data(idx_g_ff + 58);

    auto gr_xzz_zzz = pbuffer.data(idx_g_ff + 59);

    auto gr_yyy_xxx = pbuffer.data(idx_g_ff + 60);

    auto gr_yyy_xxy = pbuffer.data(idx_g_ff + 61);

    auto gr_yyy_xxz = pbuffer.data(idx_g_ff + 62);

    auto gr_yyy_xyy = pbuffer.data(idx_g_ff + 63);

    auto gr_yyy_xyz = pbuffer.data(idx_g_ff + 64);

    auto gr_yyy_xzz = pbuffer.data(idx_g_ff + 65);

    auto gr_yyy_yyy = pbuffer.data(idx_g_ff + 66);

    auto gr_yyy_yyz = pbuffer.data(idx_g_ff + 67);

    auto gr_yyy_yzz = pbuffer.data(idx_g_ff + 68);

    auto gr_yyy_zzz = pbuffer.data(idx_g_ff + 69);

    auto gr_yyz_xxx = pbuffer.data(idx_g_ff + 70);

    auto gr_yyz_xxy = pbuffer.data(idx_g_ff + 71);

    auto gr_yyz_xxz = pbuffer.data(idx_g_ff + 72);

    auto gr_yyz_xyy = pbuffer.data(idx_g_ff + 73);

    auto gr_yyz_xyz = pbuffer.data(idx_g_ff + 74);

    auto gr_yyz_xzz = pbuffer.data(idx_g_ff + 75);

    auto gr_yyz_yyy = pbuffer.data(idx_g_ff + 76);

    auto gr_yyz_yyz = pbuffer.data(idx_g_ff + 77);

    auto gr_yyz_yzz = pbuffer.data(idx_g_ff + 78);

    auto gr_yyz_zzz = pbuffer.data(idx_g_ff + 79);

    auto gr_yzz_xxx = pbuffer.data(idx_g_ff + 80);

    auto gr_yzz_xxy = pbuffer.data(idx_g_ff + 81);

    auto gr_yzz_xxz = pbuffer.data(idx_g_ff + 82);

    auto gr_yzz_xyy = pbuffer.data(idx_g_ff + 83);

    auto gr_yzz_xyz = pbuffer.data(idx_g_ff + 84);

    auto gr_yzz_xzz = pbuffer.data(idx_g_ff + 85);

    auto gr_yzz_yyy = pbuffer.data(idx_g_ff + 86);

    auto gr_yzz_yyz = pbuffer.data(idx_g_ff + 87);

    auto gr_yzz_yzz = pbuffer.data(idx_g_ff + 88);

    auto gr_yzz_zzz = pbuffer.data(idx_g_ff + 89);

    auto gr_zzz_xxx = pbuffer.data(idx_g_ff + 90);

    auto gr_zzz_xxy = pbuffer.data(idx_g_ff + 91);

    auto gr_zzz_xxz = pbuffer.data(idx_g_ff + 92);

    auto gr_zzz_xyy = pbuffer.data(idx_g_ff + 93);

    auto gr_zzz_xyz = pbuffer.data(idx_g_ff + 94);

    auto gr_zzz_xzz = pbuffer.data(idx_g_ff + 95);

    auto gr_zzz_yyy = pbuffer.data(idx_g_ff + 96);

    auto gr_zzz_yyz = pbuffer.data(idx_g_ff + 97);

    auto gr_zzz_yzz = pbuffer.data(idx_g_ff + 98);

    auto gr_zzz_zzz = pbuffer.data(idx_g_ff + 99);

    // Set up 0-10 components of targeted buffer : FF

    auto grr_x_xxx_xxx = pbuffer.data(idx_gr_ff);

    auto grr_x_xxx_xxy = pbuffer.data(idx_gr_ff + 1);

    auto grr_x_xxx_xxz = pbuffer.data(idx_gr_ff + 2);

    auto grr_x_xxx_xyy = pbuffer.data(idx_gr_ff + 3);

    auto grr_x_xxx_xyz = pbuffer.data(idx_gr_ff + 4);

    auto grr_x_xxx_xzz = pbuffer.data(idx_gr_ff + 5);

    auto grr_x_xxx_yyy = pbuffer.data(idx_gr_ff + 6);

    auto grr_x_xxx_yyz = pbuffer.data(idx_gr_ff + 7);

    auto grr_x_xxx_yzz = pbuffer.data(idx_gr_ff + 8);

    auto grr_x_xxx_zzz = pbuffer.data(idx_gr_ff + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxx, gr_xx_xxy, gr_xx_xxz, gr_xx_xyy, gr_xx_xyz, gr_xx_xzz, gr_xx_yyy, gr_xx_yyz, gr_xx_yzz, gr_xx_zzz, gr_xxx_xx, gr_xxx_xxx, gr_xxx_xxy, gr_xxx_xxz, gr_xxx_xy, gr_xxx_xyy, gr_xxx_xyz, gr_xxx_xz, gr_xxx_xzz, gr_xxx_yy, gr_xxx_yyy, gr_xxx_yyz, gr_xxx_yz, gr_xxx_yzz, gr_xxx_zz, gr_xxx_zzz, grr_x_xxx_xxx, grr_x_xxx_xxy, grr_x_xxx_xxz, grr_x_xxx_xyy, grr_x_xxx_xyz, grr_x_xxx_xzz, grr_x_xxx_yyy, grr_x_xxx_yyz, grr_x_xxx_yzz, grr_x_xxx_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxx_xxx[i] = 6.0 * ts_xx_xxx[i] * gfe2_0 + 3.0 * gr_xx_xxx[i] * gfe_0 + 6.0 * ts_xxx_xx[i] * gfe2_0 + 3.0 * gr_xxx_xx[i] * gfe_0 + 2.0 * ts_xxx_xxx[i] * gfe_0 * gc_x[i] + gr_xxx_xxx[i] * gc_x[i];

        grr_x_xxx_xxy[i] = 6.0 * ts_xx_xxy[i] * gfe2_0 + 3.0 * gr_xx_xxy[i] * gfe_0 + 4.0 * ts_xxx_xy[i] * gfe2_0 + 2.0 * gr_xxx_xy[i] * gfe_0 + 2.0 * ts_xxx_xxy[i] * gfe_0 * gc_x[i] + gr_xxx_xxy[i] * gc_x[i];

        grr_x_xxx_xxz[i] = 6.0 * ts_xx_xxz[i] * gfe2_0 + 3.0 * gr_xx_xxz[i] * gfe_0 + 4.0 * ts_xxx_xz[i] * gfe2_0 + 2.0 * gr_xxx_xz[i] * gfe_0 + 2.0 * ts_xxx_xxz[i] * gfe_0 * gc_x[i] + gr_xxx_xxz[i] * gc_x[i];

        grr_x_xxx_xyy[i] = 6.0 * ts_xx_xyy[i] * gfe2_0 + 3.0 * gr_xx_xyy[i] * gfe_0 + 2.0 * ts_xxx_yy[i] * gfe2_0 + gr_xxx_yy[i] * gfe_0 + 2.0 * ts_xxx_xyy[i] * gfe_0 * gc_x[i] + gr_xxx_xyy[i] * gc_x[i];

        grr_x_xxx_xyz[i] = 6.0 * ts_xx_xyz[i] * gfe2_0 + 3.0 * gr_xx_xyz[i] * gfe_0 + 2.0 * ts_xxx_yz[i] * gfe2_0 + gr_xxx_yz[i] * gfe_0 + 2.0 * ts_xxx_xyz[i] * gfe_0 * gc_x[i] + gr_xxx_xyz[i] * gc_x[i];

        grr_x_xxx_xzz[i] = 6.0 * ts_xx_xzz[i] * gfe2_0 + 3.0 * gr_xx_xzz[i] * gfe_0 + 2.0 * ts_xxx_zz[i] * gfe2_0 + gr_xxx_zz[i] * gfe_0 + 2.0 * ts_xxx_xzz[i] * gfe_0 * gc_x[i] + gr_xxx_xzz[i] * gc_x[i];

        grr_x_xxx_yyy[i] = 6.0 * ts_xx_yyy[i] * gfe2_0 + 3.0 * gr_xx_yyy[i] * gfe_0 + 2.0 * ts_xxx_yyy[i] * gfe_0 * gc_x[i] + gr_xxx_yyy[i] * gc_x[i];

        grr_x_xxx_yyz[i] = 6.0 * ts_xx_yyz[i] * gfe2_0 + 3.0 * gr_xx_yyz[i] * gfe_0 + 2.0 * ts_xxx_yyz[i] * gfe_0 * gc_x[i] + gr_xxx_yyz[i] * gc_x[i];

        grr_x_xxx_yzz[i] = 6.0 * ts_xx_yzz[i] * gfe2_0 + 3.0 * gr_xx_yzz[i] * gfe_0 + 2.0 * ts_xxx_yzz[i] * gfe_0 * gc_x[i] + gr_xxx_yzz[i] * gc_x[i];

        grr_x_xxx_zzz[i] = 6.0 * ts_xx_zzz[i] * gfe2_0 + 3.0 * gr_xx_zzz[i] * gfe_0 + 2.0 * ts_xxx_zzz[i] * gfe_0 * gc_x[i] + gr_xxx_zzz[i] * gc_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto grr_x_xxy_xxx = pbuffer.data(idx_gr_ff + 10);

    auto grr_x_xxy_xxy = pbuffer.data(idx_gr_ff + 11);

    auto grr_x_xxy_xxz = pbuffer.data(idx_gr_ff + 12);

    auto grr_x_xxy_xyy = pbuffer.data(idx_gr_ff + 13);

    auto grr_x_xxy_xyz = pbuffer.data(idx_gr_ff + 14);

    auto grr_x_xxy_xzz = pbuffer.data(idx_gr_ff + 15);

    auto grr_x_xxy_yyy = pbuffer.data(idx_gr_ff + 16);

    auto grr_x_xxy_yyz = pbuffer.data(idx_gr_ff + 17);

    auto grr_x_xxy_yzz = pbuffer.data(idx_gr_ff + 18);

    auto grr_x_xxy_zzz = pbuffer.data(idx_gr_ff + 19);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xx, gr_xxy_xxx, gr_xxy_xxy, gr_xxy_xxz, gr_xxy_xy, gr_xxy_xyy, gr_xxy_xyz, gr_xxy_xz, gr_xxy_xzz, gr_xxy_yy, gr_xxy_yyy, gr_xxy_yyz, gr_xxy_yz, gr_xxy_yzz, gr_xxy_zz, gr_xxy_zzz, gr_xy_xxx, gr_xy_xxy, gr_xy_xxz, gr_xy_xyy, gr_xy_xyz, gr_xy_xzz, gr_xy_yyy, gr_xy_yyz, gr_xy_yzz, gr_xy_zzz, grr_x_xxy_xxx, grr_x_xxy_xxy, grr_x_xxy_xxz, grr_x_xxy_xyy, grr_x_xxy_xyz, grr_x_xxy_xzz, grr_x_xxy_yyy, grr_x_xxy_yyz, grr_x_xxy_yzz, grr_x_xxy_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxy_xxx[i] = 4.0 * ts_xy_xxx[i] * gfe2_0 + 2.0 * gr_xy_xxx[i] * gfe_0 + 6.0 * ts_xxy_xx[i] * gfe2_0 + 3.0 * gr_xxy_xx[i] * gfe_0 + 2.0 * ts_xxy_xxx[i] * gfe_0 * gc_x[i] + gr_xxy_xxx[i] * gc_x[i];

        grr_x_xxy_xxy[i] = 4.0 * ts_xy_xxy[i] * gfe2_0 + 2.0 * gr_xy_xxy[i] * gfe_0 + 4.0 * ts_xxy_xy[i] * gfe2_0 + 2.0 * gr_xxy_xy[i] * gfe_0 + 2.0 * ts_xxy_xxy[i] * gfe_0 * gc_x[i] + gr_xxy_xxy[i] * gc_x[i];

        grr_x_xxy_xxz[i] = 4.0 * ts_xy_xxz[i] * gfe2_0 + 2.0 * gr_xy_xxz[i] * gfe_0 + 4.0 * ts_xxy_xz[i] * gfe2_0 + 2.0 * gr_xxy_xz[i] * gfe_0 + 2.0 * ts_xxy_xxz[i] * gfe_0 * gc_x[i] + gr_xxy_xxz[i] * gc_x[i];

        grr_x_xxy_xyy[i] = 4.0 * ts_xy_xyy[i] * gfe2_0 + 2.0 * gr_xy_xyy[i] * gfe_0 + 2.0 * ts_xxy_yy[i] * gfe2_0 + gr_xxy_yy[i] * gfe_0 + 2.0 * ts_xxy_xyy[i] * gfe_0 * gc_x[i] + gr_xxy_xyy[i] * gc_x[i];

        grr_x_xxy_xyz[i] = 4.0 * ts_xy_xyz[i] * gfe2_0 + 2.0 * gr_xy_xyz[i] * gfe_0 + 2.0 * ts_xxy_yz[i] * gfe2_0 + gr_xxy_yz[i] * gfe_0 + 2.0 * ts_xxy_xyz[i] * gfe_0 * gc_x[i] + gr_xxy_xyz[i] * gc_x[i];

        grr_x_xxy_xzz[i] = 4.0 * ts_xy_xzz[i] * gfe2_0 + 2.0 * gr_xy_xzz[i] * gfe_0 + 2.0 * ts_xxy_zz[i] * gfe2_0 + gr_xxy_zz[i] * gfe_0 + 2.0 * ts_xxy_xzz[i] * gfe_0 * gc_x[i] + gr_xxy_xzz[i] * gc_x[i];

        grr_x_xxy_yyy[i] = 4.0 * ts_xy_yyy[i] * gfe2_0 + 2.0 * gr_xy_yyy[i] * gfe_0 + 2.0 * ts_xxy_yyy[i] * gfe_0 * gc_x[i] + gr_xxy_yyy[i] * gc_x[i];

        grr_x_xxy_yyz[i] = 4.0 * ts_xy_yyz[i] * gfe2_0 + 2.0 * gr_xy_yyz[i] * gfe_0 + 2.0 * ts_xxy_yyz[i] * gfe_0 * gc_x[i] + gr_xxy_yyz[i] * gc_x[i];

        grr_x_xxy_yzz[i] = 4.0 * ts_xy_yzz[i] * gfe2_0 + 2.0 * gr_xy_yzz[i] * gfe_0 + 2.0 * ts_xxy_yzz[i] * gfe_0 * gc_x[i] + gr_xxy_yzz[i] * gc_x[i];

        grr_x_xxy_zzz[i] = 4.0 * ts_xy_zzz[i] * gfe2_0 + 2.0 * gr_xy_zzz[i] * gfe_0 + 2.0 * ts_xxy_zzz[i] * gfe_0 * gc_x[i] + gr_xxy_zzz[i] * gc_x[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto grr_x_xxz_xxx = pbuffer.data(idx_gr_ff + 20);

    auto grr_x_xxz_xxy = pbuffer.data(idx_gr_ff + 21);

    auto grr_x_xxz_xxz = pbuffer.data(idx_gr_ff + 22);

    auto grr_x_xxz_xyy = pbuffer.data(idx_gr_ff + 23);

    auto grr_x_xxz_xyz = pbuffer.data(idx_gr_ff + 24);

    auto grr_x_xxz_xzz = pbuffer.data(idx_gr_ff + 25);

    auto grr_x_xxz_yyy = pbuffer.data(idx_gr_ff + 26);

    auto grr_x_xxz_yyz = pbuffer.data(idx_gr_ff + 27);

    auto grr_x_xxz_yzz = pbuffer.data(idx_gr_ff + 28);

    auto grr_x_xxz_zzz = pbuffer.data(idx_gr_ff + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xx, gr_xxz_xxx, gr_xxz_xxy, gr_xxz_xxz, gr_xxz_xy, gr_xxz_xyy, gr_xxz_xyz, gr_xxz_xz, gr_xxz_xzz, gr_xxz_yy, gr_xxz_yyy, gr_xxz_yyz, gr_xxz_yz, gr_xxz_yzz, gr_xxz_zz, gr_xxz_zzz, gr_xz_xxx, gr_xz_xxy, gr_xz_xxz, gr_xz_xyy, gr_xz_xyz, gr_xz_xzz, gr_xz_yyy, gr_xz_yyz, gr_xz_yzz, gr_xz_zzz, grr_x_xxz_xxx, grr_x_xxz_xxy, grr_x_xxz_xxz, grr_x_xxz_xyy, grr_x_xxz_xyz, grr_x_xxz_xzz, grr_x_xxz_yyy, grr_x_xxz_yyz, grr_x_xxz_yzz, grr_x_xxz_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxz_xxx[i] = 4.0 * ts_xz_xxx[i] * gfe2_0 + 2.0 * gr_xz_xxx[i] * gfe_0 + 6.0 * ts_xxz_xx[i] * gfe2_0 + 3.0 * gr_xxz_xx[i] * gfe_0 + 2.0 * ts_xxz_xxx[i] * gfe_0 * gc_x[i] + gr_xxz_xxx[i] * gc_x[i];

        grr_x_xxz_xxy[i] = 4.0 * ts_xz_xxy[i] * gfe2_0 + 2.0 * gr_xz_xxy[i] * gfe_0 + 4.0 * ts_xxz_xy[i] * gfe2_0 + 2.0 * gr_xxz_xy[i] * gfe_0 + 2.0 * ts_xxz_xxy[i] * gfe_0 * gc_x[i] + gr_xxz_xxy[i] * gc_x[i];

        grr_x_xxz_xxz[i] = 4.0 * ts_xz_xxz[i] * gfe2_0 + 2.0 * gr_xz_xxz[i] * gfe_0 + 4.0 * ts_xxz_xz[i] * gfe2_0 + 2.0 * gr_xxz_xz[i] * gfe_0 + 2.0 * ts_xxz_xxz[i] * gfe_0 * gc_x[i] + gr_xxz_xxz[i] * gc_x[i];

        grr_x_xxz_xyy[i] = 4.0 * ts_xz_xyy[i] * gfe2_0 + 2.0 * gr_xz_xyy[i] * gfe_0 + 2.0 * ts_xxz_yy[i] * gfe2_0 + gr_xxz_yy[i] * gfe_0 + 2.0 * ts_xxz_xyy[i] * gfe_0 * gc_x[i] + gr_xxz_xyy[i] * gc_x[i];

        grr_x_xxz_xyz[i] = 4.0 * ts_xz_xyz[i] * gfe2_0 + 2.0 * gr_xz_xyz[i] * gfe_0 + 2.0 * ts_xxz_yz[i] * gfe2_0 + gr_xxz_yz[i] * gfe_0 + 2.0 * ts_xxz_xyz[i] * gfe_0 * gc_x[i] + gr_xxz_xyz[i] * gc_x[i];

        grr_x_xxz_xzz[i] = 4.0 * ts_xz_xzz[i] * gfe2_0 + 2.0 * gr_xz_xzz[i] * gfe_0 + 2.0 * ts_xxz_zz[i] * gfe2_0 + gr_xxz_zz[i] * gfe_0 + 2.0 * ts_xxz_xzz[i] * gfe_0 * gc_x[i] + gr_xxz_xzz[i] * gc_x[i];

        grr_x_xxz_yyy[i] = 4.0 * ts_xz_yyy[i] * gfe2_0 + 2.0 * gr_xz_yyy[i] * gfe_0 + 2.0 * ts_xxz_yyy[i] * gfe_0 * gc_x[i] + gr_xxz_yyy[i] * gc_x[i];

        grr_x_xxz_yyz[i] = 4.0 * ts_xz_yyz[i] * gfe2_0 + 2.0 * gr_xz_yyz[i] * gfe_0 + 2.0 * ts_xxz_yyz[i] * gfe_0 * gc_x[i] + gr_xxz_yyz[i] * gc_x[i];

        grr_x_xxz_yzz[i] = 4.0 * ts_xz_yzz[i] * gfe2_0 + 2.0 * gr_xz_yzz[i] * gfe_0 + 2.0 * ts_xxz_yzz[i] * gfe_0 * gc_x[i] + gr_xxz_yzz[i] * gc_x[i];

        grr_x_xxz_zzz[i] = 4.0 * ts_xz_zzz[i] * gfe2_0 + 2.0 * gr_xz_zzz[i] * gfe_0 + 2.0 * ts_xxz_zzz[i] * gfe_0 * gc_x[i] + gr_xxz_zzz[i] * gc_x[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto grr_x_xyy_xxx = pbuffer.data(idx_gr_ff + 30);

    auto grr_x_xyy_xxy = pbuffer.data(idx_gr_ff + 31);

    auto grr_x_xyy_xxz = pbuffer.data(idx_gr_ff + 32);

    auto grr_x_xyy_xyy = pbuffer.data(idx_gr_ff + 33);

    auto grr_x_xyy_xyz = pbuffer.data(idx_gr_ff + 34);

    auto grr_x_xyy_xzz = pbuffer.data(idx_gr_ff + 35);

    auto grr_x_xyy_yyy = pbuffer.data(idx_gr_ff + 36);

    auto grr_x_xyy_yyz = pbuffer.data(idx_gr_ff + 37);

    auto grr_x_xyy_yzz = pbuffer.data(idx_gr_ff + 38);

    auto grr_x_xyy_zzz = pbuffer.data(idx_gr_ff + 39);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xx, gr_xyy_xxx, gr_xyy_xxy, gr_xyy_xxz, gr_xyy_xy, gr_xyy_xyy, gr_xyy_xyz, gr_xyy_xz, gr_xyy_xzz, gr_xyy_yy, gr_xyy_yyy, gr_xyy_yyz, gr_xyy_yz, gr_xyy_yzz, gr_xyy_zz, gr_xyy_zzz, gr_yy_xxx, gr_yy_xxy, gr_yy_xxz, gr_yy_xyy, gr_yy_xyz, gr_yy_xzz, gr_yy_yyy, gr_yy_yyz, gr_yy_yzz, gr_yy_zzz, grr_x_xyy_xxx, grr_x_xyy_xxy, grr_x_xyy_xxz, grr_x_xyy_xyy, grr_x_xyy_xyz, grr_x_xyy_xzz, grr_x_xyy_yyy, grr_x_xyy_yyz, grr_x_xyy_yzz, grr_x_xyy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyy_xxx[i] = 2.0 * ts_yy_xxx[i] * gfe2_0 + gr_yy_xxx[i] * gfe_0 + 6.0 * ts_xyy_xx[i] * gfe2_0 + 3.0 * gr_xyy_xx[i] * gfe_0 + 2.0 * ts_xyy_xxx[i] * gfe_0 * gc_x[i] + gr_xyy_xxx[i] * gc_x[i];

        grr_x_xyy_xxy[i] = 2.0 * ts_yy_xxy[i] * gfe2_0 + gr_yy_xxy[i] * gfe_0 + 4.0 * ts_xyy_xy[i] * gfe2_0 + 2.0 * gr_xyy_xy[i] * gfe_0 + 2.0 * ts_xyy_xxy[i] * gfe_0 * gc_x[i] + gr_xyy_xxy[i] * gc_x[i];

        grr_x_xyy_xxz[i] = 2.0 * ts_yy_xxz[i] * gfe2_0 + gr_yy_xxz[i] * gfe_0 + 4.0 * ts_xyy_xz[i] * gfe2_0 + 2.0 * gr_xyy_xz[i] * gfe_0 + 2.0 * ts_xyy_xxz[i] * gfe_0 * gc_x[i] + gr_xyy_xxz[i] * gc_x[i];

        grr_x_xyy_xyy[i] = 2.0 * ts_yy_xyy[i] * gfe2_0 + gr_yy_xyy[i] * gfe_0 + 2.0 * ts_xyy_yy[i] * gfe2_0 + gr_xyy_yy[i] * gfe_0 + 2.0 * ts_xyy_xyy[i] * gfe_0 * gc_x[i] + gr_xyy_xyy[i] * gc_x[i];

        grr_x_xyy_xyz[i] = 2.0 * ts_yy_xyz[i] * gfe2_0 + gr_yy_xyz[i] * gfe_0 + 2.0 * ts_xyy_yz[i] * gfe2_0 + gr_xyy_yz[i] * gfe_0 + 2.0 * ts_xyy_xyz[i] * gfe_0 * gc_x[i] + gr_xyy_xyz[i] * gc_x[i];

        grr_x_xyy_xzz[i] = 2.0 * ts_yy_xzz[i] * gfe2_0 + gr_yy_xzz[i] * gfe_0 + 2.0 * ts_xyy_zz[i] * gfe2_0 + gr_xyy_zz[i] * gfe_0 + 2.0 * ts_xyy_xzz[i] * gfe_0 * gc_x[i] + gr_xyy_xzz[i] * gc_x[i];

        grr_x_xyy_yyy[i] = 2.0 * ts_yy_yyy[i] * gfe2_0 + gr_yy_yyy[i] * gfe_0 + 2.0 * ts_xyy_yyy[i] * gfe_0 * gc_x[i] + gr_xyy_yyy[i] * gc_x[i];

        grr_x_xyy_yyz[i] = 2.0 * ts_yy_yyz[i] * gfe2_0 + gr_yy_yyz[i] * gfe_0 + 2.0 * ts_xyy_yyz[i] * gfe_0 * gc_x[i] + gr_xyy_yyz[i] * gc_x[i];

        grr_x_xyy_yzz[i] = 2.0 * ts_yy_yzz[i] * gfe2_0 + gr_yy_yzz[i] * gfe_0 + 2.0 * ts_xyy_yzz[i] * gfe_0 * gc_x[i] + gr_xyy_yzz[i] * gc_x[i];

        grr_x_xyy_zzz[i] = 2.0 * ts_yy_zzz[i] * gfe2_0 + gr_yy_zzz[i] * gfe_0 + 2.0 * ts_xyy_zzz[i] * gfe_0 * gc_x[i] + gr_xyy_zzz[i] * gc_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto grr_x_xyz_xxx = pbuffer.data(idx_gr_ff + 40);

    auto grr_x_xyz_xxy = pbuffer.data(idx_gr_ff + 41);

    auto grr_x_xyz_xxz = pbuffer.data(idx_gr_ff + 42);

    auto grr_x_xyz_xyy = pbuffer.data(idx_gr_ff + 43);

    auto grr_x_xyz_xyz = pbuffer.data(idx_gr_ff + 44);

    auto grr_x_xyz_xzz = pbuffer.data(idx_gr_ff + 45);

    auto grr_x_xyz_yyy = pbuffer.data(idx_gr_ff + 46);

    auto grr_x_xyz_yyz = pbuffer.data(idx_gr_ff + 47);

    auto grr_x_xyz_yzz = pbuffer.data(idx_gr_ff + 48);

    auto grr_x_xyz_zzz = pbuffer.data(idx_gr_ff + 49);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xx, gr_xyz_xxx, gr_xyz_xxy, gr_xyz_xxz, gr_xyz_xy, gr_xyz_xyy, gr_xyz_xyz, gr_xyz_xz, gr_xyz_xzz, gr_xyz_yy, gr_xyz_yyy, gr_xyz_yyz, gr_xyz_yz, gr_xyz_yzz, gr_xyz_zz, gr_xyz_zzz, gr_yz_xxx, gr_yz_xxy, gr_yz_xxz, gr_yz_xyy, gr_yz_xyz, gr_yz_xzz, gr_yz_yyy, gr_yz_yyz, gr_yz_yzz, gr_yz_zzz, grr_x_xyz_xxx, grr_x_xyz_xxy, grr_x_xyz_xxz, grr_x_xyz_xyy, grr_x_xyz_xyz, grr_x_xyz_xzz, grr_x_xyz_yyy, grr_x_xyz_yyz, grr_x_xyz_yzz, grr_x_xyz_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyz_xxx[i] = 2.0 * ts_yz_xxx[i] * gfe2_0 + gr_yz_xxx[i] * gfe_0 + 6.0 * ts_xyz_xx[i] * gfe2_0 + 3.0 * gr_xyz_xx[i] * gfe_0 + 2.0 * ts_xyz_xxx[i] * gfe_0 * gc_x[i] + gr_xyz_xxx[i] * gc_x[i];

        grr_x_xyz_xxy[i] = 2.0 * ts_yz_xxy[i] * gfe2_0 + gr_yz_xxy[i] * gfe_0 + 4.0 * ts_xyz_xy[i] * gfe2_0 + 2.0 * gr_xyz_xy[i] * gfe_0 + 2.0 * ts_xyz_xxy[i] * gfe_0 * gc_x[i] + gr_xyz_xxy[i] * gc_x[i];

        grr_x_xyz_xxz[i] = 2.0 * ts_yz_xxz[i] * gfe2_0 + gr_yz_xxz[i] * gfe_0 + 4.0 * ts_xyz_xz[i] * gfe2_0 + 2.0 * gr_xyz_xz[i] * gfe_0 + 2.0 * ts_xyz_xxz[i] * gfe_0 * gc_x[i] + gr_xyz_xxz[i] * gc_x[i];

        grr_x_xyz_xyy[i] = 2.0 * ts_yz_xyy[i] * gfe2_0 + gr_yz_xyy[i] * gfe_0 + 2.0 * ts_xyz_yy[i] * gfe2_0 + gr_xyz_yy[i] * gfe_0 + 2.0 * ts_xyz_xyy[i] * gfe_0 * gc_x[i] + gr_xyz_xyy[i] * gc_x[i];

        grr_x_xyz_xyz[i] = 2.0 * ts_yz_xyz[i] * gfe2_0 + gr_yz_xyz[i] * gfe_0 + 2.0 * ts_xyz_yz[i] * gfe2_0 + gr_xyz_yz[i] * gfe_0 + 2.0 * ts_xyz_xyz[i] * gfe_0 * gc_x[i] + gr_xyz_xyz[i] * gc_x[i];

        grr_x_xyz_xzz[i] = 2.0 * ts_yz_xzz[i] * gfe2_0 + gr_yz_xzz[i] * gfe_0 + 2.0 * ts_xyz_zz[i] * gfe2_0 + gr_xyz_zz[i] * gfe_0 + 2.0 * ts_xyz_xzz[i] * gfe_0 * gc_x[i] + gr_xyz_xzz[i] * gc_x[i];

        grr_x_xyz_yyy[i] = 2.0 * ts_yz_yyy[i] * gfe2_0 + gr_yz_yyy[i] * gfe_0 + 2.0 * ts_xyz_yyy[i] * gfe_0 * gc_x[i] + gr_xyz_yyy[i] * gc_x[i];

        grr_x_xyz_yyz[i] = 2.0 * ts_yz_yyz[i] * gfe2_0 + gr_yz_yyz[i] * gfe_0 + 2.0 * ts_xyz_yyz[i] * gfe_0 * gc_x[i] + gr_xyz_yyz[i] * gc_x[i];

        grr_x_xyz_yzz[i] = 2.0 * ts_yz_yzz[i] * gfe2_0 + gr_yz_yzz[i] * gfe_0 + 2.0 * ts_xyz_yzz[i] * gfe_0 * gc_x[i] + gr_xyz_yzz[i] * gc_x[i];

        grr_x_xyz_zzz[i] = 2.0 * ts_yz_zzz[i] * gfe2_0 + gr_yz_zzz[i] * gfe_0 + 2.0 * ts_xyz_zzz[i] * gfe_0 * gc_x[i] + gr_xyz_zzz[i] * gc_x[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto grr_x_xzz_xxx = pbuffer.data(idx_gr_ff + 50);

    auto grr_x_xzz_xxy = pbuffer.data(idx_gr_ff + 51);

    auto grr_x_xzz_xxz = pbuffer.data(idx_gr_ff + 52);

    auto grr_x_xzz_xyy = pbuffer.data(idx_gr_ff + 53);

    auto grr_x_xzz_xyz = pbuffer.data(idx_gr_ff + 54);

    auto grr_x_xzz_xzz = pbuffer.data(idx_gr_ff + 55);

    auto grr_x_xzz_yyy = pbuffer.data(idx_gr_ff + 56);

    auto grr_x_xzz_yyz = pbuffer.data(idx_gr_ff + 57);

    auto grr_x_xzz_yzz = pbuffer.data(idx_gr_ff + 58);

    auto grr_x_xzz_zzz = pbuffer.data(idx_gr_ff + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xx, gr_xzz_xxx, gr_xzz_xxy, gr_xzz_xxz, gr_xzz_xy, gr_xzz_xyy, gr_xzz_xyz, gr_xzz_xz, gr_xzz_xzz, gr_xzz_yy, gr_xzz_yyy, gr_xzz_yyz, gr_xzz_yz, gr_xzz_yzz, gr_xzz_zz, gr_xzz_zzz, gr_zz_xxx, gr_zz_xxy, gr_zz_xxz, gr_zz_xyy, gr_zz_xyz, gr_zz_xzz, gr_zz_yyy, gr_zz_yyz, gr_zz_yzz, gr_zz_zzz, grr_x_xzz_xxx, grr_x_xzz_xxy, grr_x_xzz_xxz, grr_x_xzz_xyy, grr_x_xzz_xyz, grr_x_xzz_xzz, grr_x_xzz_yyy, grr_x_xzz_yyz, grr_x_xzz_yzz, grr_x_xzz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xzz_xxx[i] = 2.0 * ts_zz_xxx[i] * gfe2_0 + gr_zz_xxx[i] * gfe_0 + 6.0 * ts_xzz_xx[i] * gfe2_0 + 3.0 * gr_xzz_xx[i] * gfe_0 + 2.0 * ts_xzz_xxx[i] * gfe_0 * gc_x[i] + gr_xzz_xxx[i] * gc_x[i];

        grr_x_xzz_xxy[i] = 2.0 * ts_zz_xxy[i] * gfe2_0 + gr_zz_xxy[i] * gfe_0 + 4.0 * ts_xzz_xy[i] * gfe2_0 + 2.0 * gr_xzz_xy[i] * gfe_0 + 2.0 * ts_xzz_xxy[i] * gfe_0 * gc_x[i] + gr_xzz_xxy[i] * gc_x[i];

        grr_x_xzz_xxz[i] = 2.0 * ts_zz_xxz[i] * gfe2_0 + gr_zz_xxz[i] * gfe_0 + 4.0 * ts_xzz_xz[i] * gfe2_0 + 2.0 * gr_xzz_xz[i] * gfe_0 + 2.0 * ts_xzz_xxz[i] * gfe_0 * gc_x[i] + gr_xzz_xxz[i] * gc_x[i];

        grr_x_xzz_xyy[i] = 2.0 * ts_zz_xyy[i] * gfe2_0 + gr_zz_xyy[i] * gfe_0 + 2.0 * ts_xzz_yy[i] * gfe2_0 + gr_xzz_yy[i] * gfe_0 + 2.0 * ts_xzz_xyy[i] * gfe_0 * gc_x[i] + gr_xzz_xyy[i] * gc_x[i];

        grr_x_xzz_xyz[i] = 2.0 * ts_zz_xyz[i] * gfe2_0 + gr_zz_xyz[i] * gfe_0 + 2.0 * ts_xzz_yz[i] * gfe2_0 + gr_xzz_yz[i] * gfe_0 + 2.0 * ts_xzz_xyz[i] * gfe_0 * gc_x[i] + gr_xzz_xyz[i] * gc_x[i];

        grr_x_xzz_xzz[i] = 2.0 * ts_zz_xzz[i] * gfe2_0 + gr_zz_xzz[i] * gfe_0 + 2.0 * ts_xzz_zz[i] * gfe2_0 + gr_xzz_zz[i] * gfe_0 + 2.0 * ts_xzz_xzz[i] * gfe_0 * gc_x[i] + gr_xzz_xzz[i] * gc_x[i];

        grr_x_xzz_yyy[i] = 2.0 * ts_zz_yyy[i] * gfe2_0 + gr_zz_yyy[i] * gfe_0 + 2.0 * ts_xzz_yyy[i] * gfe_0 * gc_x[i] + gr_xzz_yyy[i] * gc_x[i];

        grr_x_xzz_yyz[i] = 2.0 * ts_zz_yyz[i] * gfe2_0 + gr_zz_yyz[i] * gfe_0 + 2.0 * ts_xzz_yyz[i] * gfe_0 * gc_x[i] + gr_xzz_yyz[i] * gc_x[i];

        grr_x_xzz_yzz[i] = 2.0 * ts_zz_yzz[i] * gfe2_0 + gr_zz_yzz[i] * gfe_0 + 2.0 * ts_xzz_yzz[i] * gfe_0 * gc_x[i] + gr_xzz_yzz[i] * gc_x[i];

        grr_x_xzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe2_0 + gr_zz_zzz[i] * gfe_0 + 2.0 * ts_xzz_zzz[i] * gfe_0 * gc_x[i] + gr_xzz_zzz[i] * gc_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto grr_x_yyy_xxx = pbuffer.data(idx_gr_ff + 60);

    auto grr_x_yyy_xxy = pbuffer.data(idx_gr_ff + 61);

    auto grr_x_yyy_xxz = pbuffer.data(idx_gr_ff + 62);

    auto grr_x_yyy_xyy = pbuffer.data(idx_gr_ff + 63);

    auto grr_x_yyy_xyz = pbuffer.data(idx_gr_ff + 64);

    auto grr_x_yyy_xzz = pbuffer.data(idx_gr_ff + 65);

    auto grr_x_yyy_yyy = pbuffer.data(idx_gr_ff + 66);

    auto grr_x_yyy_yyz = pbuffer.data(idx_gr_ff + 67);

    auto grr_x_yyy_yzz = pbuffer.data(idx_gr_ff + 68);

    auto grr_x_yyy_zzz = pbuffer.data(idx_gr_ff + 69);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xx, gr_yyy_xxx, gr_yyy_xxy, gr_yyy_xxz, gr_yyy_xy, gr_yyy_xyy, gr_yyy_xyz, gr_yyy_xz, gr_yyy_xzz, gr_yyy_yy, gr_yyy_yyy, gr_yyy_yyz, gr_yyy_yz, gr_yyy_yzz, gr_yyy_zz, gr_yyy_zzz, grr_x_yyy_xxx, grr_x_yyy_xxy, grr_x_yyy_xxz, grr_x_yyy_xyy, grr_x_yyy_xyz, grr_x_yyy_xzz, grr_x_yyy_yyy, grr_x_yyy_yyz, grr_x_yyy_yzz, grr_x_yyy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyy_xxx[i] = 6.0 * ts_yyy_xx[i] * gfe2_0 + 3.0 * gr_yyy_xx[i] * gfe_0 + 2.0 * ts_yyy_xxx[i] * gfe_0 * gc_x[i] + gr_yyy_xxx[i] * gc_x[i];

        grr_x_yyy_xxy[i] = 4.0 * ts_yyy_xy[i] * gfe2_0 + 2.0 * gr_yyy_xy[i] * gfe_0 + 2.0 * ts_yyy_xxy[i] * gfe_0 * gc_x[i] + gr_yyy_xxy[i] * gc_x[i];

        grr_x_yyy_xxz[i] = 4.0 * ts_yyy_xz[i] * gfe2_0 + 2.0 * gr_yyy_xz[i] * gfe_0 + 2.0 * ts_yyy_xxz[i] * gfe_0 * gc_x[i] + gr_yyy_xxz[i] * gc_x[i];

        grr_x_yyy_xyy[i] = 2.0 * ts_yyy_yy[i] * gfe2_0 + gr_yyy_yy[i] * gfe_0 + 2.0 * ts_yyy_xyy[i] * gfe_0 * gc_x[i] + gr_yyy_xyy[i] * gc_x[i];

        grr_x_yyy_xyz[i] = 2.0 * ts_yyy_yz[i] * gfe2_0 + gr_yyy_yz[i] * gfe_0 + 2.0 * ts_yyy_xyz[i] * gfe_0 * gc_x[i] + gr_yyy_xyz[i] * gc_x[i];

        grr_x_yyy_xzz[i] = 2.0 * ts_yyy_zz[i] * gfe2_0 + gr_yyy_zz[i] * gfe_0 + 2.0 * ts_yyy_xzz[i] * gfe_0 * gc_x[i] + gr_yyy_xzz[i] * gc_x[i];

        grr_x_yyy_yyy[i] = 2.0 * ts_yyy_yyy[i] * gfe_0 * gc_x[i] + gr_yyy_yyy[i] * gc_x[i];

        grr_x_yyy_yyz[i] = 2.0 * ts_yyy_yyz[i] * gfe_0 * gc_x[i] + gr_yyy_yyz[i] * gc_x[i];

        grr_x_yyy_yzz[i] = 2.0 * ts_yyy_yzz[i] * gfe_0 * gc_x[i] + gr_yyy_yzz[i] * gc_x[i];

        grr_x_yyy_zzz[i] = 2.0 * ts_yyy_zzz[i] * gfe_0 * gc_x[i] + gr_yyy_zzz[i] * gc_x[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto grr_x_yyz_xxx = pbuffer.data(idx_gr_ff + 70);

    auto grr_x_yyz_xxy = pbuffer.data(idx_gr_ff + 71);

    auto grr_x_yyz_xxz = pbuffer.data(idx_gr_ff + 72);

    auto grr_x_yyz_xyy = pbuffer.data(idx_gr_ff + 73);

    auto grr_x_yyz_xyz = pbuffer.data(idx_gr_ff + 74);

    auto grr_x_yyz_xzz = pbuffer.data(idx_gr_ff + 75);

    auto grr_x_yyz_yyy = pbuffer.data(idx_gr_ff + 76);

    auto grr_x_yyz_yyz = pbuffer.data(idx_gr_ff + 77);

    auto grr_x_yyz_yzz = pbuffer.data(idx_gr_ff + 78);

    auto grr_x_yyz_zzz = pbuffer.data(idx_gr_ff + 79);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xx, gr_yyz_xxx, gr_yyz_xxy, gr_yyz_xxz, gr_yyz_xy, gr_yyz_xyy, gr_yyz_xyz, gr_yyz_xz, gr_yyz_xzz, gr_yyz_yy, gr_yyz_yyy, gr_yyz_yyz, gr_yyz_yz, gr_yyz_yzz, gr_yyz_zz, gr_yyz_zzz, grr_x_yyz_xxx, grr_x_yyz_xxy, grr_x_yyz_xxz, grr_x_yyz_xyy, grr_x_yyz_xyz, grr_x_yyz_xzz, grr_x_yyz_yyy, grr_x_yyz_yyz, grr_x_yyz_yzz, grr_x_yyz_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyz_xxx[i] = 6.0 * ts_yyz_xx[i] * gfe2_0 + 3.0 * gr_yyz_xx[i] * gfe_0 + 2.0 * ts_yyz_xxx[i] * gfe_0 * gc_x[i] + gr_yyz_xxx[i] * gc_x[i];

        grr_x_yyz_xxy[i] = 4.0 * ts_yyz_xy[i] * gfe2_0 + 2.0 * gr_yyz_xy[i] * gfe_0 + 2.0 * ts_yyz_xxy[i] * gfe_0 * gc_x[i] + gr_yyz_xxy[i] * gc_x[i];

        grr_x_yyz_xxz[i] = 4.0 * ts_yyz_xz[i] * gfe2_0 + 2.0 * gr_yyz_xz[i] * gfe_0 + 2.0 * ts_yyz_xxz[i] * gfe_0 * gc_x[i] + gr_yyz_xxz[i] * gc_x[i];

        grr_x_yyz_xyy[i] = 2.0 * ts_yyz_yy[i] * gfe2_0 + gr_yyz_yy[i] * gfe_0 + 2.0 * ts_yyz_xyy[i] * gfe_0 * gc_x[i] + gr_yyz_xyy[i] * gc_x[i];

        grr_x_yyz_xyz[i] = 2.0 * ts_yyz_yz[i] * gfe2_0 + gr_yyz_yz[i] * gfe_0 + 2.0 * ts_yyz_xyz[i] * gfe_0 * gc_x[i] + gr_yyz_xyz[i] * gc_x[i];

        grr_x_yyz_xzz[i] = 2.0 * ts_yyz_zz[i] * gfe2_0 + gr_yyz_zz[i] * gfe_0 + 2.0 * ts_yyz_xzz[i] * gfe_0 * gc_x[i] + gr_yyz_xzz[i] * gc_x[i];

        grr_x_yyz_yyy[i] = 2.0 * ts_yyz_yyy[i] * gfe_0 * gc_x[i] + gr_yyz_yyy[i] * gc_x[i];

        grr_x_yyz_yyz[i] = 2.0 * ts_yyz_yyz[i] * gfe_0 * gc_x[i] + gr_yyz_yyz[i] * gc_x[i];

        grr_x_yyz_yzz[i] = 2.0 * ts_yyz_yzz[i] * gfe_0 * gc_x[i] + gr_yyz_yzz[i] * gc_x[i];

        grr_x_yyz_zzz[i] = 2.0 * ts_yyz_zzz[i] * gfe_0 * gc_x[i] + gr_yyz_zzz[i] * gc_x[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto grr_x_yzz_xxx = pbuffer.data(idx_gr_ff + 80);

    auto grr_x_yzz_xxy = pbuffer.data(idx_gr_ff + 81);

    auto grr_x_yzz_xxz = pbuffer.data(idx_gr_ff + 82);

    auto grr_x_yzz_xyy = pbuffer.data(idx_gr_ff + 83);

    auto grr_x_yzz_xyz = pbuffer.data(idx_gr_ff + 84);

    auto grr_x_yzz_xzz = pbuffer.data(idx_gr_ff + 85);

    auto grr_x_yzz_yyy = pbuffer.data(idx_gr_ff + 86);

    auto grr_x_yzz_yyz = pbuffer.data(idx_gr_ff + 87);

    auto grr_x_yzz_yzz = pbuffer.data(idx_gr_ff + 88);

    auto grr_x_yzz_zzz = pbuffer.data(idx_gr_ff + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xx, gr_yzz_xxx, gr_yzz_xxy, gr_yzz_xxz, gr_yzz_xy, gr_yzz_xyy, gr_yzz_xyz, gr_yzz_xz, gr_yzz_xzz, gr_yzz_yy, gr_yzz_yyy, gr_yzz_yyz, gr_yzz_yz, gr_yzz_yzz, gr_yzz_zz, gr_yzz_zzz, grr_x_yzz_xxx, grr_x_yzz_xxy, grr_x_yzz_xxz, grr_x_yzz_xyy, grr_x_yzz_xyz, grr_x_yzz_xzz, grr_x_yzz_yyy, grr_x_yzz_yyz, grr_x_yzz_yzz, grr_x_yzz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yzz_xxx[i] = 6.0 * ts_yzz_xx[i] * gfe2_0 + 3.0 * gr_yzz_xx[i] * gfe_0 + 2.0 * ts_yzz_xxx[i] * gfe_0 * gc_x[i] + gr_yzz_xxx[i] * gc_x[i];

        grr_x_yzz_xxy[i] = 4.0 * ts_yzz_xy[i] * gfe2_0 + 2.0 * gr_yzz_xy[i] * gfe_0 + 2.0 * ts_yzz_xxy[i] * gfe_0 * gc_x[i] + gr_yzz_xxy[i] * gc_x[i];

        grr_x_yzz_xxz[i] = 4.0 * ts_yzz_xz[i] * gfe2_0 + 2.0 * gr_yzz_xz[i] * gfe_0 + 2.0 * ts_yzz_xxz[i] * gfe_0 * gc_x[i] + gr_yzz_xxz[i] * gc_x[i];

        grr_x_yzz_xyy[i] = 2.0 * ts_yzz_yy[i] * gfe2_0 + gr_yzz_yy[i] * gfe_0 + 2.0 * ts_yzz_xyy[i] * gfe_0 * gc_x[i] + gr_yzz_xyy[i] * gc_x[i];

        grr_x_yzz_xyz[i] = 2.0 * ts_yzz_yz[i] * gfe2_0 + gr_yzz_yz[i] * gfe_0 + 2.0 * ts_yzz_xyz[i] * gfe_0 * gc_x[i] + gr_yzz_xyz[i] * gc_x[i];

        grr_x_yzz_xzz[i] = 2.0 * ts_yzz_zz[i] * gfe2_0 + gr_yzz_zz[i] * gfe_0 + 2.0 * ts_yzz_xzz[i] * gfe_0 * gc_x[i] + gr_yzz_xzz[i] * gc_x[i];

        grr_x_yzz_yyy[i] = 2.0 * ts_yzz_yyy[i] * gfe_0 * gc_x[i] + gr_yzz_yyy[i] * gc_x[i];

        grr_x_yzz_yyz[i] = 2.0 * ts_yzz_yyz[i] * gfe_0 * gc_x[i] + gr_yzz_yyz[i] * gc_x[i];

        grr_x_yzz_yzz[i] = 2.0 * ts_yzz_yzz[i] * gfe_0 * gc_x[i] + gr_yzz_yzz[i] * gc_x[i];

        grr_x_yzz_zzz[i] = 2.0 * ts_yzz_zzz[i] * gfe_0 * gc_x[i] + gr_yzz_zzz[i] * gc_x[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto grr_x_zzz_xxx = pbuffer.data(idx_gr_ff + 90);

    auto grr_x_zzz_xxy = pbuffer.data(idx_gr_ff + 91);

    auto grr_x_zzz_xxz = pbuffer.data(idx_gr_ff + 92);

    auto grr_x_zzz_xyy = pbuffer.data(idx_gr_ff + 93);

    auto grr_x_zzz_xyz = pbuffer.data(idx_gr_ff + 94);

    auto grr_x_zzz_xzz = pbuffer.data(idx_gr_ff + 95);

    auto grr_x_zzz_yyy = pbuffer.data(idx_gr_ff + 96);

    auto grr_x_zzz_yyz = pbuffer.data(idx_gr_ff + 97);

    auto grr_x_zzz_yzz = pbuffer.data(idx_gr_ff + 98);

    auto grr_x_zzz_zzz = pbuffer.data(idx_gr_ff + 99);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xx, gr_zzz_xxx, gr_zzz_xxy, gr_zzz_xxz, gr_zzz_xy, gr_zzz_xyy, gr_zzz_xyz, gr_zzz_xz, gr_zzz_xzz, gr_zzz_yy, gr_zzz_yyy, gr_zzz_yyz, gr_zzz_yz, gr_zzz_yzz, gr_zzz_zz, gr_zzz_zzz, grr_x_zzz_xxx, grr_x_zzz_xxy, grr_x_zzz_xxz, grr_x_zzz_xyy, grr_x_zzz_xyz, grr_x_zzz_xzz, grr_x_zzz_yyy, grr_x_zzz_yyz, grr_x_zzz_yzz, grr_x_zzz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zzz_xxx[i] = 6.0 * ts_zzz_xx[i] * gfe2_0 + 3.0 * gr_zzz_xx[i] * gfe_0 + 2.0 * ts_zzz_xxx[i] * gfe_0 * gc_x[i] + gr_zzz_xxx[i] * gc_x[i];

        grr_x_zzz_xxy[i] = 4.0 * ts_zzz_xy[i] * gfe2_0 + 2.0 * gr_zzz_xy[i] * gfe_0 + 2.0 * ts_zzz_xxy[i] * gfe_0 * gc_x[i] + gr_zzz_xxy[i] * gc_x[i];

        grr_x_zzz_xxz[i] = 4.0 * ts_zzz_xz[i] * gfe2_0 + 2.0 * gr_zzz_xz[i] * gfe_0 + 2.0 * ts_zzz_xxz[i] * gfe_0 * gc_x[i] + gr_zzz_xxz[i] * gc_x[i];

        grr_x_zzz_xyy[i] = 2.0 * ts_zzz_yy[i] * gfe2_0 + gr_zzz_yy[i] * gfe_0 + 2.0 * ts_zzz_xyy[i] * gfe_0 * gc_x[i] + gr_zzz_xyy[i] * gc_x[i];

        grr_x_zzz_xyz[i] = 2.0 * ts_zzz_yz[i] * gfe2_0 + gr_zzz_yz[i] * gfe_0 + 2.0 * ts_zzz_xyz[i] * gfe_0 * gc_x[i] + gr_zzz_xyz[i] * gc_x[i];

        grr_x_zzz_xzz[i] = 2.0 * ts_zzz_zz[i] * gfe2_0 + gr_zzz_zz[i] * gfe_0 + 2.0 * ts_zzz_xzz[i] * gfe_0 * gc_x[i] + gr_zzz_xzz[i] * gc_x[i];

        grr_x_zzz_yyy[i] = 2.0 * ts_zzz_yyy[i] * gfe_0 * gc_x[i] + gr_zzz_yyy[i] * gc_x[i];

        grr_x_zzz_yyz[i] = 2.0 * ts_zzz_yyz[i] * gfe_0 * gc_x[i] + gr_zzz_yyz[i] * gc_x[i];

        grr_x_zzz_yzz[i] = 2.0 * ts_zzz_yzz[i] * gfe_0 * gc_x[i] + gr_zzz_yzz[i] * gc_x[i];

        grr_x_zzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * gc_x[i] + gr_zzz_zzz[i] * gc_x[i];
    }

    // Set up 100-110 components of targeted buffer : FF

    auto grr_y_xxx_xxx = pbuffer.data(idx_gr_ff + 100);

    auto grr_y_xxx_xxy = pbuffer.data(idx_gr_ff + 101);

    auto grr_y_xxx_xxz = pbuffer.data(idx_gr_ff + 102);

    auto grr_y_xxx_xyy = pbuffer.data(idx_gr_ff + 103);

    auto grr_y_xxx_xyz = pbuffer.data(idx_gr_ff + 104);

    auto grr_y_xxx_xzz = pbuffer.data(idx_gr_ff + 105);

    auto grr_y_xxx_yyy = pbuffer.data(idx_gr_ff + 106);

    auto grr_y_xxx_yyz = pbuffer.data(idx_gr_ff + 107);

    auto grr_y_xxx_yzz = pbuffer.data(idx_gr_ff + 108);

    auto grr_y_xxx_zzz = pbuffer.data(idx_gr_ff + 109);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xx, gr_xxx_xxx, gr_xxx_xxy, gr_xxx_xxz, gr_xxx_xy, gr_xxx_xyy, gr_xxx_xyz, gr_xxx_xz, gr_xxx_xzz, gr_xxx_yy, gr_xxx_yyy, gr_xxx_yyz, gr_xxx_yz, gr_xxx_yzz, gr_xxx_zz, gr_xxx_zzz, grr_y_xxx_xxx, grr_y_xxx_xxy, grr_y_xxx_xxz, grr_y_xxx_xyy, grr_y_xxx_xyz, grr_y_xxx_xzz, grr_y_xxx_yyy, grr_y_xxx_yyz, grr_y_xxx_yzz, grr_y_xxx_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxx_xxx[i] = 2.0 * ts_xxx_xxx[i] * gfe_0 * gc_y[i] + gr_xxx_xxx[i] * gc_y[i];

        grr_y_xxx_xxy[i] = 2.0 * ts_xxx_xx[i] * gfe2_0 + gr_xxx_xx[i] * gfe_0 + 2.0 * ts_xxx_xxy[i] * gfe_0 * gc_y[i] + gr_xxx_xxy[i] * gc_y[i];

        grr_y_xxx_xxz[i] = 2.0 * ts_xxx_xxz[i] * gfe_0 * gc_y[i] + gr_xxx_xxz[i] * gc_y[i];

        grr_y_xxx_xyy[i] = 4.0 * ts_xxx_xy[i] * gfe2_0 + 2.0 * gr_xxx_xy[i] * gfe_0 + 2.0 * ts_xxx_xyy[i] * gfe_0 * gc_y[i] + gr_xxx_xyy[i] * gc_y[i];

        grr_y_xxx_xyz[i] = 2.0 * ts_xxx_xz[i] * gfe2_0 + gr_xxx_xz[i] * gfe_0 + 2.0 * ts_xxx_xyz[i] * gfe_0 * gc_y[i] + gr_xxx_xyz[i] * gc_y[i];

        grr_y_xxx_xzz[i] = 2.0 * ts_xxx_xzz[i] * gfe_0 * gc_y[i] + gr_xxx_xzz[i] * gc_y[i];

        grr_y_xxx_yyy[i] = 6.0 * ts_xxx_yy[i] * gfe2_0 + 3.0 * gr_xxx_yy[i] * gfe_0 + 2.0 * ts_xxx_yyy[i] * gfe_0 * gc_y[i] + gr_xxx_yyy[i] * gc_y[i];

        grr_y_xxx_yyz[i] = 4.0 * ts_xxx_yz[i] * gfe2_0 + 2.0 * gr_xxx_yz[i] * gfe_0 + 2.0 * ts_xxx_yyz[i] * gfe_0 * gc_y[i] + gr_xxx_yyz[i] * gc_y[i];

        grr_y_xxx_yzz[i] = 2.0 * ts_xxx_zz[i] * gfe2_0 + gr_xxx_zz[i] * gfe_0 + 2.0 * ts_xxx_yzz[i] * gfe_0 * gc_y[i] + gr_xxx_yzz[i] * gc_y[i];

        grr_y_xxx_zzz[i] = 2.0 * ts_xxx_zzz[i] * gfe_0 * gc_y[i] + gr_xxx_zzz[i] * gc_y[i];
    }

    // Set up 110-120 components of targeted buffer : FF

    auto grr_y_xxy_xxx = pbuffer.data(idx_gr_ff + 110);

    auto grr_y_xxy_xxy = pbuffer.data(idx_gr_ff + 111);

    auto grr_y_xxy_xxz = pbuffer.data(idx_gr_ff + 112);

    auto grr_y_xxy_xyy = pbuffer.data(idx_gr_ff + 113);

    auto grr_y_xxy_xyz = pbuffer.data(idx_gr_ff + 114);

    auto grr_y_xxy_xzz = pbuffer.data(idx_gr_ff + 115);

    auto grr_y_xxy_yyy = pbuffer.data(idx_gr_ff + 116);

    auto grr_y_xxy_yyz = pbuffer.data(idx_gr_ff + 117);

    auto grr_y_xxy_yzz = pbuffer.data(idx_gr_ff + 118);

    auto grr_y_xxy_zzz = pbuffer.data(idx_gr_ff + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxx, gr_xx_xxy, gr_xx_xxz, gr_xx_xyy, gr_xx_xyz, gr_xx_xzz, gr_xx_yyy, gr_xx_yyz, gr_xx_yzz, gr_xx_zzz, gr_xxy_xx, gr_xxy_xxx, gr_xxy_xxy, gr_xxy_xxz, gr_xxy_xy, gr_xxy_xyy, gr_xxy_xyz, gr_xxy_xz, gr_xxy_xzz, gr_xxy_yy, gr_xxy_yyy, gr_xxy_yyz, gr_xxy_yz, gr_xxy_yzz, gr_xxy_zz, gr_xxy_zzz, grr_y_xxy_xxx, grr_y_xxy_xxy, grr_y_xxy_xxz, grr_y_xxy_xyy, grr_y_xxy_xyz, grr_y_xxy_xzz, grr_y_xxy_yyy, grr_y_xxy_yyz, grr_y_xxy_yzz, grr_y_xxy_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxy_xxx[i] = 2.0 * ts_xx_xxx[i] * gfe2_0 + gr_xx_xxx[i] * gfe_0 + 2.0 * ts_xxy_xxx[i] * gfe_0 * gc_y[i] + gr_xxy_xxx[i] * gc_y[i];

        grr_y_xxy_xxy[i] = 2.0 * ts_xx_xxy[i] * gfe2_0 + gr_xx_xxy[i] * gfe_0 + 2.0 * ts_xxy_xx[i] * gfe2_0 + gr_xxy_xx[i] * gfe_0 + 2.0 * ts_xxy_xxy[i] * gfe_0 * gc_y[i] + gr_xxy_xxy[i] * gc_y[i];

        grr_y_xxy_xxz[i] = 2.0 * ts_xx_xxz[i] * gfe2_0 + gr_xx_xxz[i] * gfe_0 + 2.0 * ts_xxy_xxz[i] * gfe_0 * gc_y[i] + gr_xxy_xxz[i] * gc_y[i];

        grr_y_xxy_xyy[i] = 2.0 * ts_xx_xyy[i] * gfe2_0 + gr_xx_xyy[i] * gfe_0 + 4.0 * ts_xxy_xy[i] * gfe2_0 + 2.0 * gr_xxy_xy[i] * gfe_0 + 2.0 * ts_xxy_xyy[i] * gfe_0 * gc_y[i] + gr_xxy_xyy[i] * gc_y[i];

        grr_y_xxy_xyz[i] = 2.0 * ts_xx_xyz[i] * gfe2_0 + gr_xx_xyz[i] * gfe_0 + 2.0 * ts_xxy_xz[i] * gfe2_0 + gr_xxy_xz[i] * gfe_0 + 2.0 * ts_xxy_xyz[i] * gfe_0 * gc_y[i] + gr_xxy_xyz[i] * gc_y[i];

        grr_y_xxy_xzz[i] = 2.0 * ts_xx_xzz[i] * gfe2_0 + gr_xx_xzz[i] * gfe_0 + 2.0 * ts_xxy_xzz[i] * gfe_0 * gc_y[i] + gr_xxy_xzz[i] * gc_y[i];

        grr_y_xxy_yyy[i] = 2.0 * ts_xx_yyy[i] * gfe2_0 + gr_xx_yyy[i] * gfe_0 + 6.0 * ts_xxy_yy[i] * gfe2_0 + 3.0 * gr_xxy_yy[i] * gfe_0 + 2.0 * ts_xxy_yyy[i] * gfe_0 * gc_y[i] + gr_xxy_yyy[i] * gc_y[i];

        grr_y_xxy_yyz[i] = 2.0 * ts_xx_yyz[i] * gfe2_0 + gr_xx_yyz[i] * gfe_0 + 4.0 * ts_xxy_yz[i] * gfe2_0 + 2.0 * gr_xxy_yz[i] * gfe_0 + 2.0 * ts_xxy_yyz[i] * gfe_0 * gc_y[i] + gr_xxy_yyz[i] * gc_y[i];

        grr_y_xxy_yzz[i] = 2.0 * ts_xx_yzz[i] * gfe2_0 + gr_xx_yzz[i] * gfe_0 + 2.0 * ts_xxy_zz[i] * gfe2_0 + gr_xxy_zz[i] * gfe_0 + 2.0 * ts_xxy_yzz[i] * gfe_0 * gc_y[i] + gr_xxy_yzz[i] * gc_y[i];

        grr_y_xxy_zzz[i] = 2.0 * ts_xx_zzz[i] * gfe2_0 + gr_xx_zzz[i] * gfe_0 + 2.0 * ts_xxy_zzz[i] * gfe_0 * gc_y[i] + gr_xxy_zzz[i] * gc_y[i];
    }

    // Set up 120-130 components of targeted buffer : FF

    auto grr_y_xxz_xxx = pbuffer.data(idx_gr_ff + 120);

    auto grr_y_xxz_xxy = pbuffer.data(idx_gr_ff + 121);

    auto grr_y_xxz_xxz = pbuffer.data(idx_gr_ff + 122);

    auto grr_y_xxz_xyy = pbuffer.data(idx_gr_ff + 123);

    auto grr_y_xxz_xyz = pbuffer.data(idx_gr_ff + 124);

    auto grr_y_xxz_xzz = pbuffer.data(idx_gr_ff + 125);

    auto grr_y_xxz_yyy = pbuffer.data(idx_gr_ff + 126);

    auto grr_y_xxz_yyz = pbuffer.data(idx_gr_ff + 127);

    auto grr_y_xxz_yzz = pbuffer.data(idx_gr_ff + 128);

    auto grr_y_xxz_zzz = pbuffer.data(idx_gr_ff + 129);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xx, gr_xxz_xxx, gr_xxz_xxy, gr_xxz_xxz, gr_xxz_xy, gr_xxz_xyy, gr_xxz_xyz, gr_xxz_xz, gr_xxz_xzz, gr_xxz_yy, gr_xxz_yyy, gr_xxz_yyz, gr_xxz_yz, gr_xxz_yzz, gr_xxz_zz, gr_xxz_zzz, grr_y_xxz_xxx, grr_y_xxz_xxy, grr_y_xxz_xxz, grr_y_xxz_xyy, grr_y_xxz_xyz, grr_y_xxz_xzz, grr_y_xxz_yyy, grr_y_xxz_yyz, grr_y_xxz_yzz, grr_y_xxz_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxz_xxx[i] = 2.0 * ts_xxz_xxx[i] * gfe_0 * gc_y[i] + gr_xxz_xxx[i] * gc_y[i];

        grr_y_xxz_xxy[i] = 2.0 * ts_xxz_xx[i] * gfe2_0 + gr_xxz_xx[i] * gfe_0 + 2.0 * ts_xxz_xxy[i] * gfe_0 * gc_y[i] + gr_xxz_xxy[i] * gc_y[i];

        grr_y_xxz_xxz[i] = 2.0 * ts_xxz_xxz[i] * gfe_0 * gc_y[i] + gr_xxz_xxz[i] * gc_y[i];

        grr_y_xxz_xyy[i] = 4.0 * ts_xxz_xy[i] * gfe2_0 + 2.0 * gr_xxz_xy[i] * gfe_0 + 2.0 * ts_xxz_xyy[i] * gfe_0 * gc_y[i] + gr_xxz_xyy[i] * gc_y[i];

        grr_y_xxz_xyz[i] = 2.0 * ts_xxz_xz[i] * gfe2_0 + gr_xxz_xz[i] * gfe_0 + 2.0 * ts_xxz_xyz[i] * gfe_0 * gc_y[i] + gr_xxz_xyz[i] * gc_y[i];

        grr_y_xxz_xzz[i] = 2.0 * ts_xxz_xzz[i] * gfe_0 * gc_y[i] + gr_xxz_xzz[i] * gc_y[i];

        grr_y_xxz_yyy[i] = 6.0 * ts_xxz_yy[i] * gfe2_0 + 3.0 * gr_xxz_yy[i] * gfe_0 + 2.0 * ts_xxz_yyy[i] * gfe_0 * gc_y[i] + gr_xxz_yyy[i] * gc_y[i];

        grr_y_xxz_yyz[i] = 4.0 * ts_xxz_yz[i] * gfe2_0 + 2.0 * gr_xxz_yz[i] * gfe_0 + 2.0 * ts_xxz_yyz[i] * gfe_0 * gc_y[i] + gr_xxz_yyz[i] * gc_y[i];

        grr_y_xxz_yzz[i] = 2.0 * ts_xxz_zz[i] * gfe2_0 + gr_xxz_zz[i] * gfe_0 + 2.0 * ts_xxz_yzz[i] * gfe_0 * gc_y[i] + gr_xxz_yzz[i] * gc_y[i];

        grr_y_xxz_zzz[i] = 2.0 * ts_xxz_zzz[i] * gfe_0 * gc_y[i] + gr_xxz_zzz[i] * gc_y[i];
    }

    // Set up 130-140 components of targeted buffer : FF

    auto grr_y_xyy_xxx = pbuffer.data(idx_gr_ff + 130);

    auto grr_y_xyy_xxy = pbuffer.data(idx_gr_ff + 131);

    auto grr_y_xyy_xxz = pbuffer.data(idx_gr_ff + 132);

    auto grr_y_xyy_xyy = pbuffer.data(idx_gr_ff + 133);

    auto grr_y_xyy_xyz = pbuffer.data(idx_gr_ff + 134);

    auto grr_y_xyy_xzz = pbuffer.data(idx_gr_ff + 135);

    auto grr_y_xyy_yyy = pbuffer.data(idx_gr_ff + 136);

    auto grr_y_xyy_yyz = pbuffer.data(idx_gr_ff + 137);

    auto grr_y_xyy_yzz = pbuffer.data(idx_gr_ff + 138);

    auto grr_y_xyy_zzz = pbuffer.data(idx_gr_ff + 139);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxx, gr_xy_xxy, gr_xy_xxz, gr_xy_xyy, gr_xy_xyz, gr_xy_xzz, gr_xy_yyy, gr_xy_yyz, gr_xy_yzz, gr_xy_zzz, gr_xyy_xx, gr_xyy_xxx, gr_xyy_xxy, gr_xyy_xxz, gr_xyy_xy, gr_xyy_xyy, gr_xyy_xyz, gr_xyy_xz, gr_xyy_xzz, gr_xyy_yy, gr_xyy_yyy, gr_xyy_yyz, gr_xyy_yz, gr_xyy_yzz, gr_xyy_zz, gr_xyy_zzz, grr_y_xyy_xxx, grr_y_xyy_xxy, grr_y_xyy_xxz, grr_y_xyy_xyy, grr_y_xyy_xyz, grr_y_xyy_xzz, grr_y_xyy_yyy, grr_y_xyy_yyz, grr_y_xyy_yzz, grr_y_xyy_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyy_xxx[i] = 4.0 * ts_xy_xxx[i] * gfe2_0 + 2.0 * gr_xy_xxx[i] * gfe_0 + 2.0 * ts_xyy_xxx[i] * gfe_0 * gc_y[i] + gr_xyy_xxx[i] * gc_y[i];

        grr_y_xyy_xxy[i] = 4.0 * ts_xy_xxy[i] * gfe2_0 + 2.0 * gr_xy_xxy[i] * gfe_0 + 2.0 * ts_xyy_xx[i] * gfe2_0 + gr_xyy_xx[i] * gfe_0 + 2.0 * ts_xyy_xxy[i] * gfe_0 * gc_y[i] + gr_xyy_xxy[i] * gc_y[i];

        grr_y_xyy_xxz[i] = 4.0 * ts_xy_xxz[i] * gfe2_0 + 2.0 * gr_xy_xxz[i] * gfe_0 + 2.0 * ts_xyy_xxz[i] * gfe_0 * gc_y[i] + gr_xyy_xxz[i] * gc_y[i];

        grr_y_xyy_xyy[i] = 4.0 * ts_xy_xyy[i] * gfe2_0 + 2.0 * gr_xy_xyy[i] * gfe_0 + 4.0 * ts_xyy_xy[i] * gfe2_0 + 2.0 * gr_xyy_xy[i] * gfe_0 + 2.0 * ts_xyy_xyy[i] * gfe_0 * gc_y[i] + gr_xyy_xyy[i] * gc_y[i];

        grr_y_xyy_xyz[i] = 4.0 * ts_xy_xyz[i] * gfe2_0 + 2.0 * gr_xy_xyz[i] * gfe_0 + 2.0 * ts_xyy_xz[i] * gfe2_0 + gr_xyy_xz[i] * gfe_0 + 2.0 * ts_xyy_xyz[i] * gfe_0 * gc_y[i] + gr_xyy_xyz[i] * gc_y[i];

        grr_y_xyy_xzz[i] = 4.0 * ts_xy_xzz[i] * gfe2_0 + 2.0 * gr_xy_xzz[i] * gfe_0 + 2.0 * ts_xyy_xzz[i] * gfe_0 * gc_y[i] + gr_xyy_xzz[i] * gc_y[i];

        grr_y_xyy_yyy[i] = 4.0 * ts_xy_yyy[i] * gfe2_0 + 2.0 * gr_xy_yyy[i] * gfe_0 + 6.0 * ts_xyy_yy[i] * gfe2_0 + 3.0 * gr_xyy_yy[i] * gfe_0 + 2.0 * ts_xyy_yyy[i] * gfe_0 * gc_y[i] + gr_xyy_yyy[i] * gc_y[i];

        grr_y_xyy_yyz[i] = 4.0 * ts_xy_yyz[i] * gfe2_0 + 2.0 * gr_xy_yyz[i] * gfe_0 + 4.0 * ts_xyy_yz[i] * gfe2_0 + 2.0 * gr_xyy_yz[i] * gfe_0 + 2.0 * ts_xyy_yyz[i] * gfe_0 * gc_y[i] + gr_xyy_yyz[i] * gc_y[i];

        grr_y_xyy_yzz[i] = 4.0 * ts_xy_yzz[i] * gfe2_0 + 2.0 * gr_xy_yzz[i] * gfe_0 + 2.0 * ts_xyy_zz[i] * gfe2_0 + gr_xyy_zz[i] * gfe_0 + 2.0 * ts_xyy_yzz[i] * gfe_0 * gc_y[i] + gr_xyy_yzz[i] * gc_y[i];

        grr_y_xyy_zzz[i] = 4.0 * ts_xy_zzz[i] * gfe2_0 + 2.0 * gr_xy_zzz[i] * gfe_0 + 2.0 * ts_xyy_zzz[i] * gfe_0 * gc_y[i] + gr_xyy_zzz[i] * gc_y[i];
    }

    // Set up 140-150 components of targeted buffer : FF

    auto grr_y_xyz_xxx = pbuffer.data(idx_gr_ff + 140);

    auto grr_y_xyz_xxy = pbuffer.data(idx_gr_ff + 141);

    auto grr_y_xyz_xxz = pbuffer.data(idx_gr_ff + 142);

    auto grr_y_xyz_xyy = pbuffer.data(idx_gr_ff + 143);

    auto grr_y_xyz_xyz = pbuffer.data(idx_gr_ff + 144);

    auto grr_y_xyz_xzz = pbuffer.data(idx_gr_ff + 145);

    auto grr_y_xyz_yyy = pbuffer.data(idx_gr_ff + 146);

    auto grr_y_xyz_yyz = pbuffer.data(idx_gr_ff + 147);

    auto grr_y_xyz_yzz = pbuffer.data(idx_gr_ff + 148);

    auto grr_y_xyz_zzz = pbuffer.data(idx_gr_ff + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xx, gr_xyz_xxx, gr_xyz_xxy, gr_xyz_xxz, gr_xyz_xy, gr_xyz_xyy, gr_xyz_xyz, gr_xyz_xz, gr_xyz_xzz, gr_xyz_yy, gr_xyz_yyy, gr_xyz_yyz, gr_xyz_yz, gr_xyz_yzz, gr_xyz_zz, gr_xyz_zzz, gr_xz_xxx, gr_xz_xxy, gr_xz_xxz, gr_xz_xyy, gr_xz_xyz, gr_xz_xzz, gr_xz_yyy, gr_xz_yyz, gr_xz_yzz, gr_xz_zzz, grr_y_xyz_xxx, grr_y_xyz_xxy, grr_y_xyz_xxz, grr_y_xyz_xyy, grr_y_xyz_xyz, grr_y_xyz_xzz, grr_y_xyz_yyy, grr_y_xyz_yyz, grr_y_xyz_yzz, grr_y_xyz_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyz_xxx[i] = 2.0 * ts_xz_xxx[i] * gfe2_0 + gr_xz_xxx[i] * gfe_0 + 2.0 * ts_xyz_xxx[i] * gfe_0 * gc_y[i] + gr_xyz_xxx[i] * gc_y[i];

        grr_y_xyz_xxy[i] = 2.0 * ts_xz_xxy[i] * gfe2_0 + gr_xz_xxy[i] * gfe_0 + 2.0 * ts_xyz_xx[i] * gfe2_0 + gr_xyz_xx[i] * gfe_0 + 2.0 * ts_xyz_xxy[i] * gfe_0 * gc_y[i] + gr_xyz_xxy[i] * gc_y[i];

        grr_y_xyz_xxz[i] = 2.0 * ts_xz_xxz[i] * gfe2_0 + gr_xz_xxz[i] * gfe_0 + 2.0 * ts_xyz_xxz[i] * gfe_0 * gc_y[i] + gr_xyz_xxz[i] * gc_y[i];

        grr_y_xyz_xyy[i] = 2.0 * ts_xz_xyy[i] * gfe2_0 + gr_xz_xyy[i] * gfe_0 + 4.0 * ts_xyz_xy[i] * gfe2_0 + 2.0 * gr_xyz_xy[i] * gfe_0 + 2.0 * ts_xyz_xyy[i] * gfe_0 * gc_y[i] + gr_xyz_xyy[i] * gc_y[i];

        grr_y_xyz_xyz[i] = 2.0 * ts_xz_xyz[i] * gfe2_0 + gr_xz_xyz[i] * gfe_0 + 2.0 * ts_xyz_xz[i] * gfe2_0 + gr_xyz_xz[i] * gfe_0 + 2.0 * ts_xyz_xyz[i] * gfe_0 * gc_y[i] + gr_xyz_xyz[i] * gc_y[i];

        grr_y_xyz_xzz[i] = 2.0 * ts_xz_xzz[i] * gfe2_0 + gr_xz_xzz[i] * gfe_0 + 2.0 * ts_xyz_xzz[i] * gfe_0 * gc_y[i] + gr_xyz_xzz[i] * gc_y[i];

        grr_y_xyz_yyy[i] = 2.0 * ts_xz_yyy[i] * gfe2_0 + gr_xz_yyy[i] * gfe_0 + 6.0 * ts_xyz_yy[i] * gfe2_0 + 3.0 * gr_xyz_yy[i] * gfe_0 + 2.0 * ts_xyz_yyy[i] * gfe_0 * gc_y[i] + gr_xyz_yyy[i] * gc_y[i];

        grr_y_xyz_yyz[i] = 2.0 * ts_xz_yyz[i] * gfe2_0 + gr_xz_yyz[i] * gfe_0 + 4.0 * ts_xyz_yz[i] * gfe2_0 + 2.0 * gr_xyz_yz[i] * gfe_0 + 2.0 * ts_xyz_yyz[i] * gfe_0 * gc_y[i] + gr_xyz_yyz[i] * gc_y[i];

        grr_y_xyz_yzz[i] = 2.0 * ts_xz_yzz[i] * gfe2_0 + gr_xz_yzz[i] * gfe_0 + 2.0 * ts_xyz_zz[i] * gfe2_0 + gr_xyz_zz[i] * gfe_0 + 2.0 * ts_xyz_yzz[i] * gfe_0 * gc_y[i] + gr_xyz_yzz[i] * gc_y[i];

        grr_y_xyz_zzz[i] = 2.0 * ts_xz_zzz[i] * gfe2_0 + gr_xz_zzz[i] * gfe_0 + 2.0 * ts_xyz_zzz[i] * gfe_0 * gc_y[i] + gr_xyz_zzz[i] * gc_y[i];
    }

    // Set up 150-160 components of targeted buffer : FF

    auto grr_y_xzz_xxx = pbuffer.data(idx_gr_ff + 150);

    auto grr_y_xzz_xxy = pbuffer.data(idx_gr_ff + 151);

    auto grr_y_xzz_xxz = pbuffer.data(idx_gr_ff + 152);

    auto grr_y_xzz_xyy = pbuffer.data(idx_gr_ff + 153);

    auto grr_y_xzz_xyz = pbuffer.data(idx_gr_ff + 154);

    auto grr_y_xzz_xzz = pbuffer.data(idx_gr_ff + 155);

    auto grr_y_xzz_yyy = pbuffer.data(idx_gr_ff + 156);

    auto grr_y_xzz_yyz = pbuffer.data(idx_gr_ff + 157);

    auto grr_y_xzz_yzz = pbuffer.data(idx_gr_ff + 158);

    auto grr_y_xzz_zzz = pbuffer.data(idx_gr_ff + 159);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xx, gr_xzz_xxx, gr_xzz_xxy, gr_xzz_xxz, gr_xzz_xy, gr_xzz_xyy, gr_xzz_xyz, gr_xzz_xz, gr_xzz_xzz, gr_xzz_yy, gr_xzz_yyy, gr_xzz_yyz, gr_xzz_yz, gr_xzz_yzz, gr_xzz_zz, gr_xzz_zzz, grr_y_xzz_xxx, grr_y_xzz_xxy, grr_y_xzz_xxz, grr_y_xzz_xyy, grr_y_xzz_xyz, grr_y_xzz_xzz, grr_y_xzz_yyy, grr_y_xzz_yyz, grr_y_xzz_yzz, grr_y_xzz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xzz_xxx[i] = 2.0 * ts_xzz_xxx[i] * gfe_0 * gc_y[i] + gr_xzz_xxx[i] * gc_y[i];

        grr_y_xzz_xxy[i] = 2.0 * ts_xzz_xx[i] * gfe2_0 + gr_xzz_xx[i] * gfe_0 + 2.0 * ts_xzz_xxy[i] * gfe_0 * gc_y[i] + gr_xzz_xxy[i] * gc_y[i];

        grr_y_xzz_xxz[i] = 2.0 * ts_xzz_xxz[i] * gfe_0 * gc_y[i] + gr_xzz_xxz[i] * gc_y[i];

        grr_y_xzz_xyy[i] = 4.0 * ts_xzz_xy[i] * gfe2_0 + 2.0 * gr_xzz_xy[i] * gfe_0 + 2.0 * ts_xzz_xyy[i] * gfe_0 * gc_y[i] + gr_xzz_xyy[i] * gc_y[i];

        grr_y_xzz_xyz[i] = 2.0 * ts_xzz_xz[i] * gfe2_0 + gr_xzz_xz[i] * gfe_0 + 2.0 * ts_xzz_xyz[i] * gfe_0 * gc_y[i] + gr_xzz_xyz[i] * gc_y[i];

        grr_y_xzz_xzz[i] = 2.0 * ts_xzz_xzz[i] * gfe_0 * gc_y[i] + gr_xzz_xzz[i] * gc_y[i];

        grr_y_xzz_yyy[i] = 6.0 * ts_xzz_yy[i] * gfe2_0 + 3.0 * gr_xzz_yy[i] * gfe_0 + 2.0 * ts_xzz_yyy[i] * gfe_0 * gc_y[i] + gr_xzz_yyy[i] * gc_y[i];

        grr_y_xzz_yyz[i] = 4.0 * ts_xzz_yz[i] * gfe2_0 + 2.0 * gr_xzz_yz[i] * gfe_0 + 2.0 * ts_xzz_yyz[i] * gfe_0 * gc_y[i] + gr_xzz_yyz[i] * gc_y[i];

        grr_y_xzz_yzz[i] = 2.0 * ts_xzz_zz[i] * gfe2_0 + gr_xzz_zz[i] * gfe_0 + 2.0 * ts_xzz_yzz[i] * gfe_0 * gc_y[i] + gr_xzz_yzz[i] * gc_y[i];

        grr_y_xzz_zzz[i] = 2.0 * ts_xzz_zzz[i] * gfe_0 * gc_y[i] + gr_xzz_zzz[i] * gc_y[i];
    }

    // Set up 160-170 components of targeted buffer : FF

    auto grr_y_yyy_xxx = pbuffer.data(idx_gr_ff + 160);

    auto grr_y_yyy_xxy = pbuffer.data(idx_gr_ff + 161);

    auto grr_y_yyy_xxz = pbuffer.data(idx_gr_ff + 162);

    auto grr_y_yyy_xyy = pbuffer.data(idx_gr_ff + 163);

    auto grr_y_yyy_xyz = pbuffer.data(idx_gr_ff + 164);

    auto grr_y_yyy_xzz = pbuffer.data(idx_gr_ff + 165);

    auto grr_y_yyy_yyy = pbuffer.data(idx_gr_ff + 166);

    auto grr_y_yyy_yyz = pbuffer.data(idx_gr_ff + 167);

    auto grr_y_yyy_yzz = pbuffer.data(idx_gr_ff + 168);

    auto grr_y_yyy_zzz = pbuffer.data(idx_gr_ff + 169);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxx, gr_yy_xxy, gr_yy_xxz, gr_yy_xyy, gr_yy_xyz, gr_yy_xzz, gr_yy_yyy, gr_yy_yyz, gr_yy_yzz, gr_yy_zzz, gr_yyy_xx, gr_yyy_xxx, gr_yyy_xxy, gr_yyy_xxz, gr_yyy_xy, gr_yyy_xyy, gr_yyy_xyz, gr_yyy_xz, gr_yyy_xzz, gr_yyy_yy, gr_yyy_yyy, gr_yyy_yyz, gr_yyy_yz, gr_yyy_yzz, gr_yyy_zz, gr_yyy_zzz, grr_y_yyy_xxx, grr_y_yyy_xxy, grr_y_yyy_xxz, grr_y_yyy_xyy, grr_y_yyy_xyz, grr_y_yyy_xzz, grr_y_yyy_yyy, grr_y_yyy_yyz, grr_y_yyy_yzz, grr_y_yyy_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyy_xxx[i] = 6.0 * ts_yy_xxx[i] * gfe2_0 + 3.0 * gr_yy_xxx[i] * gfe_0 + 2.0 * ts_yyy_xxx[i] * gfe_0 * gc_y[i] + gr_yyy_xxx[i] * gc_y[i];

        grr_y_yyy_xxy[i] = 6.0 * ts_yy_xxy[i] * gfe2_0 + 3.0 * gr_yy_xxy[i] * gfe_0 + 2.0 * ts_yyy_xx[i] * gfe2_0 + gr_yyy_xx[i] * gfe_0 + 2.0 * ts_yyy_xxy[i] * gfe_0 * gc_y[i] + gr_yyy_xxy[i] * gc_y[i];

        grr_y_yyy_xxz[i] = 6.0 * ts_yy_xxz[i] * gfe2_0 + 3.0 * gr_yy_xxz[i] * gfe_0 + 2.0 * ts_yyy_xxz[i] * gfe_0 * gc_y[i] + gr_yyy_xxz[i] * gc_y[i];

        grr_y_yyy_xyy[i] = 6.0 * ts_yy_xyy[i] * gfe2_0 + 3.0 * gr_yy_xyy[i] * gfe_0 + 4.0 * ts_yyy_xy[i] * gfe2_0 + 2.0 * gr_yyy_xy[i] * gfe_0 + 2.0 * ts_yyy_xyy[i] * gfe_0 * gc_y[i] + gr_yyy_xyy[i] * gc_y[i];

        grr_y_yyy_xyz[i] = 6.0 * ts_yy_xyz[i] * gfe2_0 + 3.0 * gr_yy_xyz[i] * gfe_0 + 2.0 * ts_yyy_xz[i] * gfe2_0 + gr_yyy_xz[i] * gfe_0 + 2.0 * ts_yyy_xyz[i] * gfe_0 * gc_y[i] + gr_yyy_xyz[i] * gc_y[i];

        grr_y_yyy_xzz[i] = 6.0 * ts_yy_xzz[i] * gfe2_0 + 3.0 * gr_yy_xzz[i] * gfe_0 + 2.0 * ts_yyy_xzz[i] * gfe_0 * gc_y[i] + gr_yyy_xzz[i] * gc_y[i];

        grr_y_yyy_yyy[i] = 6.0 * ts_yy_yyy[i] * gfe2_0 + 3.0 * gr_yy_yyy[i] * gfe_0 + 6.0 * ts_yyy_yy[i] * gfe2_0 + 3.0 * gr_yyy_yy[i] * gfe_0 + 2.0 * ts_yyy_yyy[i] * gfe_0 * gc_y[i] + gr_yyy_yyy[i] * gc_y[i];

        grr_y_yyy_yyz[i] = 6.0 * ts_yy_yyz[i] * gfe2_0 + 3.0 * gr_yy_yyz[i] * gfe_0 + 4.0 * ts_yyy_yz[i] * gfe2_0 + 2.0 * gr_yyy_yz[i] * gfe_0 + 2.0 * ts_yyy_yyz[i] * gfe_0 * gc_y[i] + gr_yyy_yyz[i] * gc_y[i];

        grr_y_yyy_yzz[i] = 6.0 * ts_yy_yzz[i] * gfe2_0 + 3.0 * gr_yy_yzz[i] * gfe_0 + 2.0 * ts_yyy_zz[i] * gfe2_0 + gr_yyy_zz[i] * gfe_0 + 2.0 * ts_yyy_yzz[i] * gfe_0 * gc_y[i] + gr_yyy_yzz[i] * gc_y[i];

        grr_y_yyy_zzz[i] = 6.0 * ts_yy_zzz[i] * gfe2_0 + 3.0 * gr_yy_zzz[i] * gfe_0 + 2.0 * ts_yyy_zzz[i] * gfe_0 * gc_y[i] + gr_yyy_zzz[i] * gc_y[i];
    }

    // Set up 170-180 components of targeted buffer : FF

    auto grr_y_yyz_xxx = pbuffer.data(idx_gr_ff + 170);

    auto grr_y_yyz_xxy = pbuffer.data(idx_gr_ff + 171);

    auto grr_y_yyz_xxz = pbuffer.data(idx_gr_ff + 172);

    auto grr_y_yyz_xyy = pbuffer.data(idx_gr_ff + 173);

    auto grr_y_yyz_xyz = pbuffer.data(idx_gr_ff + 174);

    auto grr_y_yyz_xzz = pbuffer.data(idx_gr_ff + 175);

    auto grr_y_yyz_yyy = pbuffer.data(idx_gr_ff + 176);

    auto grr_y_yyz_yyz = pbuffer.data(idx_gr_ff + 177);

    auto grr_y_yyz_yzz = pbuffer.data(idx_gr_ff + 178);

    auto grr_y_yyz_zzz = pbuffer.data(idx_gr_ff + 179);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xx, gr_yyz_xxx, gr_yyz_xxy, gr_yyz_xxz, gr_yyz_xy, gr_yyz_xyy, gr_yyz_xyz, gr_yyz_xz, gr_yyz_xzz, gr_yyz_yy, gr_yyz_yyy, gr_yyz_yyz, gr_yyz_yz, gr_yyz_yzz, gr_yyz_zz, gr_yyz_zzz, gr_yz_xxx, gr_yz_xxy, gr_yz_xxz, gr_yz_xyy, gr_yz_xyz, gr_yz_xzz, gr_yz_yyy, gr_yz_yyz, gr_yz_yzz, gr_yz_zzz, grr_y_yyz_xxx, grr_y_yyz_xxy, grr_y_yyz_xxz, grr_y_yyz_xyy, grr_y_yyz_xyz, grr_y_yyz_xzz, grr_y_yyz_yyy, grr_y_yyz_yyz, grr_y_yyz_yzz, grr_y_yyz_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyz_xxx[i] = 4.0 * ts_yz_xxx[i] * gfe2_0 + 2.0 * gr_yz_xxx[i] * gfe_0 + 2.0 * ts_yyz_xxx[i] * gfe_0 * gc_y[i] + gr_yyz_xxx[i] * gc_y[i];

        grr_y_yyz_xxy[i] = 4.0 * ts_yz_xxy[i] * gfe2_0 + 2.0 * gr_yz_xxy[i] * gfe_0 + 2.0 * ts_yyz_xx[i] * gfe2_0 + gr_yyz_xx[i] * gfe_0 + 2.0 * ts_yyz_xxy[i] * gfe_0 * gc_y[i] + gr_yyz_xxy[i] * gc_y[i];

        grr_y_yyz_xxz[i] = 4.0 * ts_yz_xxz[i] * gfe2_0 + 2.0 * gr_yz_xxz[i] * gfe_0 + 2.0 * ts_yyz_xxz[i] * gfe_0 * gc_y[i] + gr_yyz_xxz[i] * gc_y[i];

        grr_y_yyz_xyy[i] = 4.0 * ts_yz_xyy[i] * gfe2_0 + 2.0 * gr_yz_xyy[i] * gfe_0 + 4.0 * ts_yyz_xy[i] * gfe2_0 + 2.0 * gr_yyz_xy[i] * gfe_0 + 2.0 * ts_yyz_xyy[i] * gfe_0 * gc_y[i] + gr_yyz_xyy[i] * gc_y[i];

        grr_y_yyz_xyz[i] = 4.0 * ts_yz_xyz[i] * gfe2_0 + 2.0 * gr_yz_xyz[i] * gfe_0 + 2.0 * ts_yyz_xz[i] * gfe2_0 + gr_yyz_xz[i] * gfe_0 + 2.0 * ts_yyz_xyz[i] * gfe_0 * gc_y[i] + gr_yyz_xyz[i] * gc_y[i];

        grr_y_yyz_xzz[i] = 4.0 * ts_yz_xzz[i] * gfe2_0 + 2.0 * gr_yz_xzz[i] * gfe_0 + 2.0 * ts_yyz_xzz[i] * gfe_0 * gc_y[i] + gr_yyz_xzz[i] * gc_y[i];

        grr_y_yyz_yyy[i] = 4.0 * ts_yz_yyy[i] * gfe2_0 + 2.0 * gr_yz_yyy[i] * gfe_0 + 6.0 * ts_yyz_yy[i] * gfe2_0 + 3.0 * gr_yyz_yy[i] * gfe_0 + 2.0 * ts_yyz_yyy[i] * gfe_0 * gc_y[i] + gr_yyz_yyy[i] * gc_y[i];

        grr_y_yyz_yyz[i] = 4.0 * ts_yz_yyz[i] * gfe2_0 + 2.0 * gr_yz_yyz[i] * gfe_0 + 4.0 * ts_yyz_yz[i] * gfe2_0 + 2.0 * gr_yyz_yz[i] * gfe_0 + 2.0 * ts_yyz_yyz[i] * gfe_0 * gc_y[i] + gr_yyz_yyz[i] * gc_y[i];

        grr_y_yyz_yzz[i] = 4.0 * ts_yz_yzz[i] * gfe2_0 + 2.0 * gr_yz_yzz[i] * gfe_0 + 2.0 * ts_yyz_zz[i] * gfe2_0 + gr_yyz_zz[i] * gfe_0 + 2.0 * ts_yyz_yzz[i] * gfe_0 * gc_y[i] + gr_yyz_yzz[i] * gc_y[i];

        grr_y_yyz_zzz[i] = 4.0 * ts_yz_zzz[i] * gfe2_0 + 2.0 * gr_yz_zzz[i] * gfe_0 + 2.0 * ts_yyz_zzz[i] * gfe_0 * gc_y[i] + gr_yyz_zzz[i] * gc_y[i];
    }

    // Set up 180-190 components of targeted buffer : FF

    auto grr_y_yzz_xxx = pbuffer.data(idx_gr_ff + 180);

    auto grr_y_yzz_xxy = pbuffer.data(idx_gr_ff + 181);

    auto grr_y_yzz_xxz = pbuffer.data(idx_gr_ff + 182);

    auto grr_y_yzz_xyy = pbuffer.data(idx_gr_ff + 183);

    auto grr_y_yzz_xyz = pbuffer.data(idx_gr_ff + 184);

    auto grr_y_yzz_xzz = pbuffer.data(idx_gr_ff + 185);

    auto grr_y_yzz_yyy = pbuffer.data(idx_gr_ff + 186);

    auto grr_y_yzz_yyz = pbuffer.data(idx_gr_ff + 187);

    auto grr_y_yzz_yzz = pbuffer.data(idx_gr_ff + 188);

    auto grr_y_yzz_zzz = pbuffer.data(idx_gr_ff + 189);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xx, gr_yzz_xxx, gr_yzz_xxy, gr_yzz_xxz, gr_yzz_xy, gr_yzz_xyy, gr_yzz_xyz, gr_yzz_xz, gr_yzz_xzz, gr_yzz_yy, gr_yzz_yyy, gr_yzz_yyz, gr_yzz_yz, gr_yzz_yzz, gr_yzz_zz, gr_yzz_zzz, gr_zz_xxx, gr_zz_xxy, gr_zz_xxz, gr_zz_xyy, gr_zz_xyz, gr_zz_xzz, gr_zz_yyy, gr_zz_yyz, gr_zz_yzz, gr_zz_zzz, grr_y_yzz_xxx, grr_y_yzz_xxy, grr_y_yzz_xxz, grr_y_yzz_xyy, grr_y_yzz_xyz, grr_y_yzz_xzz, grr_y_yzz_yyy, grr_y_yzz_yyz, grr_y_yzz_yzz, grr_y_yzz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yzz_xxx[i] = 2.0 * ts_zz_xxx[i] * gfe2_0 + gr_zz_xxx[i] * gfe_0 + 2.0 * ts_yzz_xxx[i] * gfe_0 * gc_y[i] + gr_yzz_xxx[i] * gc_y[i];

        grr_y_yzz_xxy[i] = 2.0 * ts_zz_xxy[i] * gfe2_0 + gr_zz_xxy[i] * gfe_0 + 2.0 * ts_yzz_xx[i] * gfe2_0 + gr_yzz_xx[i] * gfe_0 + 2.0 * ts_yzz_xxy[i] * gfe_0 * gc_y[i] + gr_yzz_xxy[i] * gc_y[i];

        grr_y_yzz_xxz[i] = 2.0 * ts_zz_xxz[i] * gfe2_0 + gr_zz_xxz[i] * gfe_0 + 2.0 * ts_yzz_xxz[i] * gfe_0 * gc_y[i] + gr_yzz_xxz[i] * gc_y[i];

        grr_y_yzz_xyy[i] = 2.0 * ts_zz_xyy[i] * gfe2_0 + gr_zz_xyy[i] * gfe_0 + 4.0 * ts_yzz_xy[i] * gfe2_0 + 2.0 * gr_yzz_xy[i] * gfe_0 + 2.0 * ts_yzz_xyy[i] * gfe_0 * gc_y[i] + gr_yzz_xyy[i] * gc_y[i];

        grr_y_yzz_xyz[i] = 2.0 * ts_zz_xyz[i] * gfe2_0 + gr_zz_xyz[i] * gfe_0 + 2.0 * ts_yzz_xz[i] * gfe2_0 + gr_yzz_xz[i] * gfe_0 + 2.0 * ts_yzz_xyz[i] * gfe_0 * gc_y[i] + gr_yzz_xyz[i] * gc_y[i];

        grr_y_yzz_xzz[i] = 2.0 * ts_zz_xzz[i] * gfe2_0 + gr_zz_xzz[i] * gfe_0 + 2.0 * ts_yzz_xzz[i] * gfe_0 * gc_y[i] + gr_yzz_xzz[i] * gc_y[i];

        grr_y_yzz_yyy[i] = 2.0 * ts_zz_yyy[i] * gfe2_0 + gr_zz_yyy[i] * gfe_0 + 6.0 * ts_yzz_yy[i] * gfe2_0 + 3.0 * gr_yzz_yy[i] * gfe_0 + 2.0 * ts_yzz_yyy[i] * gfe_0 * gc_y[i] + gr_yzz_yyy[i] * gc_y[i];

        grr_y_yzz_yyz[i] = 2.0 * ts_zz_yyz[i] * gfe2_0 + gr_zz_yyz[i] * gfe_0 + 4.0 * ts_yzz_yz[i] * gfe2_0 + 2.0 * gr_yzz_yz[i] * gfe_0 + 2.0 * ts_yzz_yyz[i] * gfe_0 * gc_y[i] + gr_yzz_yyz[i] * gc_y[i];

        grr_y_yzz_yzz[i] = 2.0 * ts_zz_yzz[i] * gfe2_0 + gr_zz_yzz[i] * gfe_0 + 2.0 * ts_yzz_zz[i] * gfe2_0 + gr_yzz_zz[i] * gfe_0 + 2.0 * ts_yzz_yzz[i] * gfe_0 * gc_y[i] + gr_yzz_yzz[i] * gc_y[i];

        grr_y_yzz_zzz[i] = 2.0 * ts_zz_zzz[i] * gfe2_0 + gr_zz_zzz[i] * gfe_0 + 2.0 * ts_yzz_zzz[i] * gfe_0 * gc_y[i] + gr_yzz_zzz[i] * gc_y[i];
    }

    // Set up 190-200 components of targeted buffer : FF

    auto grr_y_zzz_xxx = pbuffer.data(idx_gr_ff + 190);

    auto grr_y_zzz_xxy = pbuffer.data(idx_gr_ff + 191);

    auto grr_y_zzz_xxz = pbuffer.data(idx_gr_ff + 192);

    auto grr_y_zzz_xyy = pbuffer.data(idx_gr_ff + 193);

    auto grr_y_zzz_xyz = pbuffer.data(idx_gr_ff + 194);

    auto grr_y_zzz_xzz = pbuffer.data(idx_gr_ff + 195);

    auto grr_y_zzz_yyy = pbuffer.data(idx_gr_ff + 196);

    auto grr_y_zzz_yyz = pbuffer.data(idx_gr_ff + 197);

    auto grr_y_zzz_yzz = pbuffer.data(idx_gr_ff + 198);

    auto grr_y_zzz_zzz = pbuffer.data(idx_gr_ff + 199);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xx, gr_zzz_xxx, gr_zzz_xxy, gr_zzz_xxz, gr_zzz_xy, gr_zzz_xyy, gr_zzz_xyz, gr_zzz_xz, gr_zzz_xzz, gr_zzz_yy, gr_zzz_yyy, gr_zzz_yyz, gr_zzz_yz, gr_zzz_yzz, gr_zzz_zz, gr_zzz_zzz, grr_y_zzz_xxx, grr_y_zzz_xxy, grr_y_zzz_xxz, grr_y_zzz_xyy, grr_y_zzz_xyz, grr_y_zzz_xzz, grr_y_zzz_yyy, grr_y_zzz_yyz, grr_y_zzz_yzz, grr_y_zzz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zzz_xxx[i] = 2.0 * ts_zzz_xxx[i] * gfe_0 * gc_y[i] + gr_zzz_xxx[i] * gc_y[i];

        grr_y_zzz_xxy[i] = 2.0 * ts_zzz_xx[i] * gfe2_0 + gr_zzz_xx[i] * gfe_0 + 2.0 * ts_zzz_xxy[i] * gfe_0 * gc_y[i] + gr_zzz_xxy[i] * gc_y[i];

        grr_y_zzz_xxz[i] = 2.0 * ts_zzz_xxz[i] * gfe_0 * gc_y[i] + gr_zzz_xxz[i] * gc_y[i];

        grr_y_zzz_xyy[i] = 4.0 * ts_zzz_xy[i] * gfe2_0 + 2.0 * gr_zzz_xy[i] * gfe_0 + 2.0 * ts_zzz_xyy[i] * gfe_0 * gc_y[i] + gr_zzz_xyy[i] * gc_y[i];

        grr_y_zzz_xyz[i] = 2.0 * ts_zzz_xz[i] * gfe2_0 + gr_zzz_xz[i] * gfe_0 + 2.0 * ts_zzz_xyz[i] * gfe_0 * gc_y[i] + gr_zzz_xyz[i] * gc_y[i];

        grr_y_zzz_xzz[i] = 2.0 * ts_zzz_xzz[i] * gfe_0 * gc_y[i] + gr_zzz_xzz[i] * gc_y[i];

        grr_y_zzz_yyy[i] = 6.0 * ts_zzz_yy[i] * gfe2_0 + 3.0 * gr_zzz_yy[i] * gfe_0 + 2.0 * ts_zzz_yyy[i] * gfe_0 * gc_y[i] + gr_zzz_yyy[i] * gc_y[i];

        grr_y_zzz_yyz[i] = 4.0 * ts_zzz_yz[i] * gfe2_0 + 2.0 * gr_zzz_yz[i] * gfe_0 + 2.0 * ts_zzz_yyz[i] * gfe_0 * gc_y[i] + gr_zzz_yyz[i] * gc_y[i];

        grr_y_zzz_yzz[i] = 2.0 * ts_zzz_zz[i] * gfe2_0 + gr_zzz_zz[i] * gfe_0 + 2.0 * ts_zzz_yzz[i] * gfe_0 * gc_y[i] + gr_zzz_yzz[i] * gc_y[i];

        grr_y_zzz_zzz[i] = 2.0 * ts_zzz_zzz[i] * gfe_0 * gc_y[i] + gr_zzz_zzz[i] * gc_y[i];
    }

    // Set up 200-210 components of targeted buffer : FF

    auto grr_z_xxx_xxx = pbuffer.data(idx_gr_ff + 200);

    auto grr_z_xxx_xxy = pbuffer.data(idx_gr_ff + 201);

    auto grr_z_xxx_xxz = pbuffer.data(idx_gr_ff + 202);

    auto grr_z_xxx_xyy = pbuffer.data(idx_gr_ff + 203);

    auto grr_z_xxx_xyz = pbuffer.data(idx_gr_ff + 204);

    auto grr_z_xxx_xzz = pbuffer.data(idx_gr_ff + 205);

    auto grr_z_xxx_yyy = pbuffer.data(idx_gr_ff + 206);

    auto grr_z_xxx_yyz = pbuffer.data(idx_gr_ff + 207);

    auto grr_z_xxx_yzz = pbuffer.data(idx_gr_ff + 208);

    auto grr_z_xxx_zzz = pbuffer.data(idx_gr_ff + 209);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xx, gr_xxx_xxx, gr_xxx_xxy, gr_xxx_xxz, gr_xxx_xy, gr_xxx_xyy, gr_xxx_xyz, gr_xxx_xz, gr_xxx_xzz, gr_xxx_yy, gr_xxx_yyy, gr_xxx_yyz, gr_xxx_yz, gr_xxx_yzz, gr_xxx_zz, gr_xxx_zzz, grr_z_xxx_xxx, grr_z_xxx_xxy, grr_z_xxx_xxz, grr_z_xxx_xyy, grr_z_xxx_xyz, grr_z_xxx_xzz, grr_z_xxx_yyy, grr_z_xxx_yyz, grr_z_xxx_yzz, grr_z_xxx_zzz, ts_xxx_xx, ts_xxx_xxx, ts_xxx_xxy, ts_xxx_xxz, ts_xxx_xy, ts_xxx_xyy, ts_xxx_xyz, ts_xxx_xz, ts_xxx_xzz, ts_xxx_yy, ts_xxx_yyy, ts_xxx_yyz, ts_xxx_yz, ts_xxx_yzz, ts_xxx_zz, ts_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxx_xxx[i] = 2.0 * ts_xxx_xxx[i] * gfe_0 * gc_z[i] + gr_xxx_xxx[i] * gc_z[i];

        grr_z_xxx_xxy[i] = 2.0 * ts_xxx_xxy[i] * gfe_0 * gc_z[i] + gr_xxx_xxy[i] * gc_z[i];

        grr_z_xxx_xxz[i] = 2.0 * ts_xxx_xx[i] * gfe2_0 + gr_xxx_xx[i] * gfe_0 + 2.0 * ts_xxx_xxz[i] * gfe_0 * gc_z[i] + gr_xxx_xxz[i] * gc_z[i];

        grr_z_xxx_xyy[i] = 2.0 * ts_xxx_xyy[i] * gfe_0 * gc_z[i] + gr_xxx_xyy[i] * gc_z[i];

        grr_z_xxx_xyz[i] = 2.0 * ts_xxx_xy[i] * gfe2_0 + gr_xxx_xy[i] * gfe_0 + 2.0 * ts_xxx_xyz[i] * gfe_0 * gc_z[i] + gr_xxx_xyz[i] * gc_z[i];

        grr_z_xxx_xzz[i] = 4.0 * ts_xxx_xz[i] * gfe2_0 + 2.0 * gr_xxx_xz[i] * gfe_0 + 2.0 * ts_xxx_xzz[i] * gfe_0 * gc_z[i] + gr_xxx_xzz[i] * gc_z[i];

        grr_z_xxx_yyy[i] = 2.0 * ts_xxx_yyy[i] * gfe_0 * gc_z[i] + gr_xxx_yyy[i] * gc_z[i];

        grr_z_xxx_yyz[i] = 2.0 * ts_xxx_yy[i] * gfe2_0 + gr_xxx_yy[i] * gfe_0 + 2.0 * ts_xxx_yyz[i] * gfe_0 * gc_z[i] + gr_xxx_yyz[i] * gc_z[i];

        grr_z_xxx_yzz[i] = 4.0 * ts_xxx_yz[i] * gfe2_0 + 2.0 * gr_xxx_yz[i] * gfe_0 + 2.0 * ts_xxx_yzz[i] * gfe_0 * gc_z[i] + gr_xxx_yzz[i] * gc_z[i];

        grr_z_xxx_zzz[i] = 6.0 * ts_xxx_zz[i] * gfe2_0 + 3.0 * gr_xxx_zz[i] * gfe_0 + 2.0 * ts_xxx_zzz[i] * gfe_0 * gc_z[i] + gr_xxx_zzz[i] * gc_z[i];
    }

    // Set up 210-220 components of targeted buffer : FF

    auto grr_z_xxy_xxx = pbuffer.data(idx_gr_ff + 210);

    auto grr_z_xxy_xxy = pbuffer.data(idx_gr_ff + 211);

    auto grr_z_xxy_xxz = pbuffer.data(idx_gr_ff + 212);

    auto grr_z_xxy_xyy = pbuffer.data(idx_gr_ff + 213);

    auto grr_z_xxy_xyz = pbuffer.data(idx_gr_ff + 214);

    auto grr_z_xxy_xzz = pbuffer.data(idx_gr_ff + 215);

    auto grr_z_xxy_yyy = pbuffer.data(idx_gr_ff + 216);

    auto grr_z_xxy_yyz = pbuffer.data(idx_gr_ff + 217);

    auto grr_z_xxy_yzz = pbuffer.data(idx_gr_ff + 218);

    auto grr_z_xxy_zzz = pbuffer.data(idx_gr_ff + 219);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xx, gr_xxy_xxx, gr_xxy_xxy, gr_xxy_xxz, gr_xxy_xy, gr_xxy_xyy, gr_xxy_xyz, gr_xxy_xz, gr_xxy_xzz, gr_xxy_yy, gr_xxy_yyy, gr_xxy_yyz, gr_xxy_yz, gr_xxy_yzz, gr_xxy_zz, gr_xxy_zzz, grr_z_xxy_xxx, grr_z_xxy_xxy, grr_z_xxy_xxz, grr_z_xxy_xyy, grr_z_xxy_xyz, grr_z_xxy_xzz, grr_z_xxy_yyy, grr_z_xxy_yyz, grr_z_xxy_yzz, grr_z_xxy_zzz, ts_xxy_xx, ts_xxy_xxx, ts_xxy_xxy, ts_xxy_xxz, ts_xxy_xy, ts_xxy_xyy, ts_xxy_xyz, ts_xxy_xz, ts_xxy_xzz, ts_xxy_yy, ts_xxy_yyy, ts_xxy_yyz, ts_xxy_yz, ts_xxy_yzz, ts_xxy_zz, ts_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxy_xxx[i] = 2.0 * ts_xxy_xxx[i] * gfe_0 * gc_z[i] + gr_xxy_xxx[i] * gc_z[i];

        grr_z_xxy_xxy[i] = 2.0 * ts_xxy_xxy[i] * gfe_0 * gc_z[i] + gr_xxy_xxy[i] * gc_z[i];

        grr_z_xxy_xxz[i] = 2.0 * ts_xxy_xx[i] * gfe2_0 + gr_xxy_xx[i] * gfe_0 + 2.0 * ts_xxy_xxz[i] * gfe_0 * gc_z[i] + gr_xxy_xxz[i] * gc_z[i];

        grr_z_xxy_xyy[i] = 2.0 * ts_xxy_xyy[i] * gfe_0 * gc_z[i] + gr_xxy_xyy[i] * gc_z[i];

        grr_z_xxy_xyz[i] = 2.0 * ts_xxy_xy[i] * gfe2_0 + gr_xxy_xy[i] * gfe_0 + 2.0 * ts_xxy_xyz[i] * gfe_0 * gc_z[i] + gr_xxy_xyz[i] * gc_z[i];

        grr_z_xxy_xzz[i] = 4.0 * ts_xxy_xz[i] * gfe2_0 + 2.0 * gr_xxy_xz[i] * gfe_0 + 2.0 * ts_xxy_xzz[i] * gfe_0 * gc_z[i] + gr_xxy_xzz[i] * gc_z[i];

        grr_z_xxy_yyy[i] = 2.0 * ts_xxy_yyy[i] * gfe_0 * gc_z[i] + gr_xxy_yyy[i] * gc_z[i];

        grr_z_xxy_yyz[i] = 2.0 * ts_xxy_yy[i] * gfe2_0 + gr_xxy_yy[i] * gfe_0 + 2.0 * ts_xxy_yyz[i] * gfe_0 * gc_z[i] + gr_xxy_yyz[i] * gc_z[i];

        grr_z_xxy_yzz[i] = 4.0 * ts_xxy_yz[i] * gfe2_0 + 2.0 * gr_xxy_yz[i] * gfe_0 + 2.0 * ts_xxy_yzz[i] * gfe_0 * gc_z[i] + gr_xxy_yzz[i] * gc_z[i];

        grr_z_xxy_zzz[i] = 6.0 * ts_xxy_zz[i] * gfe2_0 + 3.0 * gr_xxy_zz[i] * gfe_0 + 2.0 * ts_xxy_zzz[i] * gfe_0 * gc_z[i] + gr_xxy_zzz[i] * gc_z[i];
    }

    // Set up 220-230 components of targeted buffer : FF

    auto grr_z_xxz_xxx = pbuffer.data(idx_gr_ff + 220);

    auto grr_z_xxz_xxy = pbuffer.data(idx_gr_ff + 221);

    auto grr_z_xxz_xxz = pbuffer.data(idx_gr_ff + 222);

    auto grr_z_xxz_xyy = pbuffer.data(idx_gr_ff + 223);

    auto grr_z_xxz_xyz = pbuffer.data(idx_gr_ff + 224);

    auto grr_z_xxz_xzz = pbuffer.data(idx_gr_ff + 225);

    auto grr_z_xxz_yyy = pbuffer.data(idx_gr_ff + 226);

    auto grr_z_xxz_yyz = pbuffer.data(idx_gr_ff + 227);

    auto grr_z_xxz_yzz = pbuffer.data(idx_gr_ff + 228);

    auto grr_z_xxz_zzz = pbuffer.data(idx_gr_ff + 229);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xxx, gr_xx_xxy, gr_xx_xxz, gr_xx_xyy, gr_xx_xyz, gr_xx_xzz, gr_xx_yyy, gr_xx_yyz, gr_xx_yzz, gr_xx_zzz, gr_xxz_xx, gr_xxz_xxx, gr_xxz_xxy, gr_xxz_xxz, gr_xxz_xy, gr_xxz_xyy, gr_xxz_xyz, gr_xxz_xz, gr_xxz_xzz, gr_xxz_yy, gr_xxz_yyy, gr_xxz_yyz, gr_xxz_yz, gr_xxz_yzz, gr_xxz_zz, gr_xxz_zzz, grr_z_xxz_xxx, grr_z_xxz_xxy, grr_z_xxz_xxz, grr_z_xxz_xyy, grr_z_xxz_xyz, grr_z_xxz_xzz, grr_z_xxz_yyy, grr_z_xxz_yyz, grr_z_xxz_yzz, grr_z_xxz_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, ts_xxz_xx, ts_xxz_xxx, ts_xxz_xxy, ts_xxz_xxz, ts_xxz_xy, ts_xxz_xyy, ts_xxz_xyz, ts_xxz_xz, ts_xxz_xzz, ts_xxz_yy, ts_xxz_yyy, ts_xxz_yyz, ts_xxz_yz, ts_xxz_yzz, ts_xxz_zz, ts_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxz_xxx[i] = 2.0 * ts_xx_xxx[i] * gfe2_0 + gr_xx_xxx[i] * gfe_0 + 2.0 * ts_xxz_xxx[i] * gfe_0 * gc_z[i] + gr_xxz_xxx[i] * gc_z[i];

        grr_z_xxz_xxy[i] = 2.0 * ts_xx_xxy[i] * gfe2_0 + gr_xx_xxy[i] * gfe_0 + 2.0 * ts_xxz_xxy[i] * gfe_0 * gc_z[i] + gr_xxz_xxy[i] * gc_z[i];

        grr_z_xxz_xxz[i] = 2.0 * ts_xx_xxz[i] * gfe2_0 + gr_xx_xxz[i] * gfe_0 + 2.0 * ts_xxz_xx[i] * gfe2_0 + gr_xxz_xx[i] * gfe_0 + 2.0 * ts_xxz_xxz[i] * gfe_0 * gc_z[i] + gr_xxz_xxz[i] * gc_z[i];

        grr_z_xxz_xyy[i] = 2.0 * ts_xx_xyy[i] * gfe2_0 + gr_xx_xyy[i] * gfe_0 + 2.0 * ts_xxz_xyy[i] * gfe_0 * gc_z[i] + gr_xxz_xyy[i] * gc_z[i];

        grr_z_xxz_xyz[i] = 2.0 * ts_xx_xyz[i] * gfe2_0 + gr_xx_xyz[i] * gfe_0 + 2.0 * ts_xxz_xy[i] * gfe2_0 + gr_xxz_xy[i] * gfe_0 + 2.0 * ts_xxz_xyz[i] * gfe_0 * gc_z[i] + gr_xxz_xyz[i] * gc_z[i];

        grr_z_xxz_xzz[i] = 2.0 * ts_xx_xzz[i] * gfe2_0 + gr_xx_xzz[i] * gfe_0 + 4.0 * ts_xxz_xz[i] * gfe2_0 + 2.0 * gr_xxz_xz[i] * gfe_0 + 2.0 * ts_xxz_xzz[i] * gfe_0 * gc_z[i] + gr_xxz_xzz[i] * gc_z[i];

        grr_z_xxz_yyy[i] = 2.0 * ts_xx_yyy[i] * gfe2_0 + gr_xx_yyy[i] * gfe_0 + 2.0 * ts_xxz_yyy[i] * gfe_0 * gc_z[i] + gr_xxz_yyy[i] * gc_z[i];

        grr_z_xxz_yyz[i] = 2.0 * ts_xx_yyz[i] * gfe2_0 + gr_xx_yyz[i] * gfe_0 + 2.0 * ts_xxz_yy[i] * gfe2_0 + gr_xxz_yy[i] * gfe_0 + 2.0 * ts_xxz_yyz[i] * gfe_0 * gc_z[i] + gr_xxz_yyz[i] * gc_z[i];

        grr_z_xxz_yzz[i] = 2.0 * ts_xx_yzz[i] * gfe2_0 + gr_xx_yzz[i] * gfe_0 + 4.0 * ts_xxz_yz[i] * gfe2_0 + 2.0 * gr_xxz_yz[i] * gfe_0 + 2.0 * ts_xxz_yzz[i] * gfe_0 * gc_z[i] + gr_xxz_yzz[i] * gc_z[i];

        grr_z_xxz_zzz[i] = 2.0 * ts_xx_zzz[i] * gfe2_0 + gr_xx_zzz[i] * gfe_0 + 6.0 * ts_xxz_zz[i] * gfe2_0 + 3.0 * gr_xxz_zz[i] * gfe_0 + 2.0 * ts_xxz_zzz[i] * gfe_0 * gc_z[i] + gr_xxz_zzz[i] * gc_z[i];
    }

    // Set up 230-240 components of targeted buffer : FF

    auto grr_z_xyy_xxx = pbuffer.data(idx_gr_ff + 230);

    auto grr_z_xyy_xxy = pbuffer.data(idx_gr_ff + 231);

    auto grr_z_xyy_xxz = pbuffer.data(idx_gr_ff + 232);

    auto grr_z_xyy_xyy = pbuffer.data(idx_gr_ff + 233);

    auto grr_z_xyy_xyz = pbuffer.data(idx_gr_ff + 234);

    auto grr_z_xyy_xzz = pbuffer.data(idx_gr_ff + 235);

    auto grr_z_xyy_yyy = pbuffer.data(idx_gr_ff + 236);

    auto grr_z_xyy_yyz = pbuffer.data(idx_gr_ff + 237);

    auto grr_z_xyy_yzz = pbuffer.data(idx_gr_ff + 238);

    auto grr_z_xyy_zzz = pbuffer.data(idx_gr_ff + 239);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xx, gr_xyy_xxx, gr_xyy_xxy, gr_xyy_xxz, gr_xyy_xy, gr_xyy_xyy, gr_xyy_xyz, gr_xyy_xz, gr_xyy_xzz, gr_xyy_yy, gr_xyy_yyy, gr_xyy_yyz, gr_xyy_yz, gr_xyy_yzz, gr_xyy_zz, gr_xyy_zzz, grr_z_xyy_xxx, grr_z_xyy_xxy, grr_z_xyy_xxz, grr_z_xyy_xyy, grr_z_xyy_xyz, grr_z_xyy_xzz, grr_z_xyy_yyy, grr_z_xyy_yyz, grr_z_xyy_yzz, grr_z_xyy_zzz, ts_xyy_xx, ts_xyy_xxx, ts_xyy_xxy, ts_xyy_xxz, ts_xyy_xy, ts_xyy_xyy, ts_xyy_xyz, ts_xyy_xz, ts_xyy_xzz, ts_xyy_yy, ts_xyy_yyy, ts_xyy_yyz, ts_xyy_yz, ts_xyy_yzz, ts_xyy_zz, ts_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyy_xxx[i] = 2.0 * ts_xyy_xxx[i] * gfe_0 * gc_z[i] + gr_xyy_xxx[i] * gc_z[i];

        grr_z_xyy_xxy[i] = 2.0 * ts_xyy_xxy[i] * gfe_0 * gc_z[i] + gr_xyy_xxy[i] * gc_z[i];

        grr_z_xyy_xxz[i] = 2.0 * ts_xyy_xx[i] * gfe2_0 + gr_xyy_xx[i] * gfe_0 + 2.0 * ts_xyy_xxz[i] * gfe_0 * gc_z[i] + gr_xyy_xxz[i] * gc_z[i];

        grr_z_xyy_xyy[i] = 2.0 * ts_xyy_xyy[i] * gfe_0 * gc_z[i] + gr_xyy_xyy[i] * gc_z[i];

        grr_z_xyy_xyz[i] = 2.0 * ts_xyy_xy[i] * gfe2_0 + gr_xyy_xy[i] * gfe_0 + 2.0 * ts_xyy_xyz[i] * gfe_0 * gc_z[i] + gr_xyy_xyz[i] * gc_z[i];

        grr_z_xyy_xzz[i] = 4.0 * ts_xyy_xz[i] * gfe2_0 + 2.0 * gr_xyy_xz[i] * gfe_0 + 2.0 * ts_xyy_xzz[i] * gfe_0 * gc_z[i] + gr_xyy_xzz[i] * gc_z[i];

        grr_z_xyy_yyy[i] = 2.0 * ts_xyy_yyy[i] * gfe_0 * gc_z[i] + gr_xyy_yyy[i] * gc_z[i];

        grr_z_xyy_yyz[i] = 2.0 * ts_xyy_yy[i] * gfe2_0 + gr_xyy_yy[i] * gfe_0 + 2.0 * ts_xyy_yyz[i] * gfe_0 * gc_z[i] + gr_xyy_yyz[i] * gc_z[i];

        grr_z_xyy_yzz[i] = 4.0 * ts_xyy_yz[i] * gfe2_0 + 2.0 * gr_xyy_yz[i] * gfe_0 + 2.0 * ts_xyy_yzz[i] * gfe_0 * gc_z[i] + gr_xyy_yzz[i] * gc_z[i];

        grr_z_xyy_zzz[i] = 6.0 * ts_xyy_zz[i] * gfe2_0 + 3.0 * gr_xyy_zz[i] * gfe_0 + 2.0 * ts_xyy_zzz[i] * gfe_0 * gc_z[i] + gr_xyy_zzz[i] * gc_z[i];
    }

    // Set up 240-250 components of targeted buffer : FF

    auto grr_z_xyz_xxx = pbuffer.data(idx_gr_ff + 240);

    auto grr_z_xyz_xxy = pbuffer.data(idx_gr_ff + 241);

    auto grr_z_xyz_xxz = pbuffer.data(idx_gr_ff + 242);

    auto grr_z_xyz_xyy = pbuffer.data(idx_gr_ff + 243);

    auto grr_z_xyz_xyz = pbuffer.data(idx_gr_ff + 244);

    auto grr_z_xyz_xzz = pbuffer.data(idx_gr_ff + 245);

    auto grr_z_xyz_yyy = pbuffer.data(idx_gr_ff + 246);

    auto grr_z_xyz_yyz = pbuffer.data(idx_gr_ff + 247);

    auto grr_z_xyz_yzz = pbuffer.data(idx_gr_ff + 248);

    auto grr_z_xyz_zzz = pbuffer.data(idx_gr_ff + 249);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xxx, gr_xy_xxy, gr_xy_xxz, gr_xy_xyy, gr_xy_xyz, gr_xy_xzz, gr_xy_yyy, gr_xy_yyz, gr_xy_yzz, gr_xy_zzz, gr_xyz_xx, gr_xyz_xxx, gr_xyz_xxy, gr_xyz_xxz, gr_xyz_xy, gr_xyz_xyy, gr_xyz_xyz, gr_xyz_xz, gr_xyz_xzz, gr_xyz_yy, gr_xyz_yyy, gr_xyz_yyz, gr_xyz_yz, gr_xyz_yzz, gr_xyz_zz, gr_xyz_zzz, grr_z_xyz_xxx, grr_z_xyz_xxy, grr_z_xyz_xxz, grr_z_xyz_xyy, grr_z_xyz_xyz, grr_z_xyz_xzz, grr_z_xyz_yyy, grr_z_xyz_yyz, grr_z_xyz_yzz, grr_z_xyz_zzz, ts_xy_xxx, ts_xy_xxy, ts_xy_xxz, ts_xy_xyy, ts_xy_xyz, ts_xy_xzz, ts_xy_yyy, ts_xy_yyz, ts_xy_yzz, ts_xy_zzz, ts_xyz_xx, ts_xyz_xxx, ts_xyz_xxy, ts_xyz_xxz, ts_xyz_xy, ts_xyz_xyy, ts_xyz_xyz, ts_xyz_xz, ts_xyz_xzz, ts_xyz_yy, ts_xyz_yyy, ts_xyz_yyz, ts_xyz_yz, ts_xyz_yzz, ts_xyz_zz, ts_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyz_xxx[i] = 2.0 * ts_xy_xxx[i] * gfe2_0 + gr_xy_xxx[i] * gfe_0 + 2.0 * ts_xyz_xxx[i] * gfe_0 * gc_z[i] + gr_xyz_xxx[i] * gc_z[i];

        grr_z_xyz_xxy[i] = 2.0 * ts_xy_xxy[i] * gfe2_0 + gr_xy_xxy[i] * gfe_0 + 2.0 * ts_xyz_xxy[i] * gfe_0 * gc_z[i] + gr_xyz_xxy[i] * gc_z[i];

        grr_z_xyz_xxz[i] = 2.0 * ts_xy_xxz[i] * gfe2_0 + gr_xy_xxz[i] * gfe_0 + 2.0 * ts_xyz_xx[i] * gfe2_0 + gr_xyz_xx[i] * gfe_0 + 2.0 * ts_xyz_xxz[i] * gfe_0 * gc_z[i] + gr_xyz_xxz[i] * gc_z[i];

        grr_z_xyz_xyy[i] = 2.0 * ts_xy_xyy[i] * gfe2_0 + gr_xy_xyy[i] * gfe_0 + 2.0 * ts_xyz_xyy[i] * gfe_0 * gc_z[i] + gr_xyz_xyy[i] * gc_z[i];

        grr_z_xyz_xyz[i] = 2.0 * ts_xy_xyz[i] * gfe2_0 + gr_xy_xyz[i] * gfe_0 + 2.0 * ts_xyz_xy[i] * gfe2_0 + gr_xyz_xy[i] * gfe_0 + 2.0 * ts_xyz_xyz[i] * gfe_0 * gc_z[i] + gr_xyz_xyz[i] * gc_z[i];

        grr_z_xyz_xzz[i] = 2.0 * ts_xy_xzz[i] * gfe2_0 + gr_xy_xzz[i] * gfe_0 + 4.0 * ts_xyz_xz[i] * gfe2_0 + 2.0 * gr_xyz_xz[i] * gfe_0 + 2.0 * ts_xyz_xzz[i] * gfe_0 * gc_z[i] + gr_xyz_xzz[i] * gc_z[i];

        grr_z_xyz_yyy[i] = 2.0 * ts_xy_yyy[i] * gfe2_0 + gr_xy_yyy[i] * gfe_0 + 2.0 * ts_xyz_yyy[i] * gfe_0 * gc_z[i] + gr_xyz_yyy[i] * gc_z[i];

        grr_z_xyz_yyz[i] = 2.0 * ts_xy_yyz[i] * gfe2_0 + gr_xy_yyz[i] * gfe_0 + 2.0 * ts_xyz_yy[i] * gfe2_0 + gr_xyz_yy[i] * gfe_0 + 2.0 * ts_xyz_yyz[i] * gfe_0 * gc_z[i] + gr_xyz_yyz[i] * gc_z[i];

        grr_z_xyz_yzz[i] = 2.0 * ts_xy_yzz[i] * gfe2_0 + gr_xy_yzz[i] * gfe_0 + 4.0 * ts_xyz_yz[i] * gfe2_0 + 2.0 * gr_xyz_yz[i] * gfe_0 + 2.0 * ts_xyz_yzz[i] * gfe_0 * gc_z[i] + gr_xyz_yzz[i] * gc_z[i];

        grr_z_xyz_zzz[i] = 2.0 * ts_xy_zzz[i] * gfe2_0 + gr_xy_zzz[i] * gfe_0 + 6.0 * ts_xyz_zz[i] * gfe2_0 + 3.0 * gr_xyz_zz[i] * gfe_0 + 2.0 * ts_xyz_zzz[i] * gfe_0 * gc_z[i] + gr_xyz_zzz[i] * gc_z[i];
    }

    // Set up 250-260 components of targeted buffer : FF

    auto grr_z_xzz_xxx = pbuffer.data(idx_gr_ff + 250);

    auto grr_z_xzz_xxy = pbuffer.data(idx_gr_ff + 251);

    auto grr_z_xzz_xxz = pbuffer.data(idx_gr_ff + 252);

    auto grr_z_xzz_xyy = pbuffer.data(idx_gr_ff + 253);

    auto grr_z_xzz_xyz = pbuffer.data(idx_gr_ff + 254);

    auto grr_z_xzz_xzz = pbuffer.data(idx_gr_ff + 255);

    auto grr_z_xzz_yyy = pbuffer.data(idx_gr_ff + 256);

    auto grr_z_xzz_yyz = pbuffer.data(idx_gr_ff + 257);

    auto grr_z_xzz_yzz = pbuffer.data(idx_gr_ff + 258);

    auto grr_z_xzz_zzz = pbuffer.data(idx_gr_ff + 259);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xxx, gr_xz_xxy, gr_xz_xxz, gr_xz_xyy, gr_xz_xyz, gr_xz_xzz, gr_xz_yyy, gr_xz_yyz, gr_xz_yzz, gr_xz_zzz, gr_xzz_xx, gr_xzz_xxx, gr_xzz_xxy, gr_xzz_xxz, gr_xzz_xy, gr_xzz_xyy, gr_xzz_xyz, gr_xzz_xz, gr_xzz_xzz, gr_xzz_yy, gr_xzz_yyy, gr_xzz_yyz, gr_xzz_yz, gr_xzz_yzz, gr_xzz_zz, gr_xzz_zzz, grr_z_xzz_xxx, grr_z_xzz_xxy, grr_z_xzz_xxz, grr_z_xzz_xyy, grr_z_xzz_xyz, grr_z_xzz_xzz, grr_z_xzz_yyy, grr_z_xzz_yyz, grr_z_xzz_yzz, grr_z_xzz_zzz, ts_xz_xxx, ts_xz_xxy, ts_xz_xxz, ts_xz_xyy, ts_xz_xyz, ts_xz_xzz, ts_xz_yyy, ts_xz_yyz, ts_xz_yzz, ts_xz_zzz, ts_xzz_xx, ts_xzz_xxx, ts_xzz_xxy, ts_xzz_xxz, ts_xzz_xy, ts_xzz_xyy, ts_xzz_xyz, ts_xzz_xz, ts_xzz_xzz, ts_xzz_yy, ts_xzz_yyy, ts_xzz_yyz, ts_xzz_yz, ts_xzz_yzz, ts_xzz_zz, ts_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xzz_xxx[i] = 4.0 * ts_xz_xxx[i] * gfe2_0 + 2.0 * gr_xz_xxx[i] * gfe_0 + 2.0 * ts_xzz_xxx[i] * gfe_0 * gc_z[i] + gr_xzz_xxx[i] * gc_z[i];

        grr_z_xzz_xxy[i] = 4.0 * ts_xz_xxy[i] * gfe2_0 + 2.0 * gr_xz_xxy[i] * gfe_0 + 2.0 * ts_xzz_xxy[i] * gfe_0 * gc_z[i] + gr_xzz_xxy[i] * gc_z[i];

        grr_z_xzz_xxz[i] = 4.0 * ts_xz_xxz[i] * gfe2_0 + 2.0 * gr_xz_xxz[i] * gfe_0 + 2.0 * ts_xzz_xx[i] * gfe2_0 + gr_xzz_xx[i] * gfe_0 + 2.0 * ts_xzz_xxz[i] * gfe_0 * gc_z[i] + gr_xzz_xxz[i] * gc_z[i];

        grr_z_xzz_xyy[i] = 4.0 * ts_xz_xyy[i] * gfe2_0 + 2.0 * gr_xz_xyy[i] * gfe_0 + 2.0 * ts_xzz_xyy[i] * gfe_0 * gc_z[i] + gr_xzz_xyy[i] * gc_z[i];

        grr_z_xzz_xyz[i] = 4.0 * ts_xz_xyz[i] * gfe2_0 + 2.0 * gr_xz_xyz[i] * gfe_0 + 2.0 * ts_xzz_xy[i] * gfe2_0 + gr_xzz_xy[i] * gfe_0 + 2.0 * ts_xzz_xyz[i] * gfe_0 * gc_z[i] + gr_xzz_xyz[i] * gc_z[i];

        grr_z_xzz_xzz[i] = 4.0 * ts_xz_xzz[i] * gfe2_0 + 2.0 * gr_xz_xzz[i] * gfe_0 + 4.0 * ts_xzz_xz[i] * gfe2_0 + 2.0 * gr_xzz_xz[i] * gfe_0 + 2.0 * ts_xzz_xzz[i] * gfe_0 * gc_z[i] + gr_xzz_xzz[i] * gc_z[i];

        grr_z_xzz_yyy[i] = 4.0 * ts_xz_yyy[i] * gfe2_0 + 2.0 * gr_xz_yyy[i] * gfe_0 + 2.0 * ts_xzz_yyy[i] * gfe_0 * gc_z[i] + gr_xzz_yyy[i] * gc_z[i];

        grr_z_xzz_yyz[i] = 4.0 * ts_xz_yyz[i] * gfe2_0 + 2.0 * gr_xz_yyz[i] * gfe_0 + 2.0 * ts_xzz_yy[i] * gfe2_0 + gr_xzz_yy[i] * gfe_0 + 2.0 * ts_xzz_yyz[i] * gfe_0 * gc_z[i] + gr_xzz_yyz[i] * gc_z[i];

        grr_z_xzz_yzz[i] = 4.0 * ts_xz_yzz[i] * gfe2_0 + 2.0 * gr_xz_yzz[i] * gfe_0 + 4.0 * ts_xzz_yz[i] * gfe2_0 + 2.0 * gr_xzz_yz[i] * gfe_0 + 2.0 * ts_xzz_yzz[i] * gfe_0 * gc_z[i] + gr_xzz_yzz[i] * gc_z[i];

        grr_z_xzz_zzz[i] = 4.0 * ts_xz_zzz[i] * gfe2_0 + 2.0 * gr_xz_zzz[i] * gfe_0 + 6.0 * ts_xzz_zz[i] * gfe2_0 + 3.0 * gr_xzz_zz[i] * gfe_0 + 2.0 * ts_xzz_zzz[i] * gfe_0 * gc_z[i] + gr_xzz_zzz[i] * gc_z[i];
    }

    // Set up 260-270 components of targeted buffer : FF

    auto grr_z_yyy_xxx = pbuffer.data(idx_gr_ff + 260);

    auto grr_z_yyy_xxy = pbuffer.data(idx_gr_ff + 261);

    auto grr_z_yyy_xxz = pbuffer.data(idx_gr_ff + 262);

    auto grr_z_yyy_xyy = pbuffer.data(idx_gr_ff + 263);

    auto grr_z_yyy_xyz = pbuffer.data(idx_gr_ff + 264);

    auto grr_z_yyy_xzz = pbuffer.data(idx_gr_ff + 265);

    auto grr_z_yyy_yyy = pbuffer.data(idx_gr_ff + 266);

    auto grr_z_yyy_yyz = pbuffer.data(idx_gr_ff + 267);

    auto grr_z_yyy_yzz = pbuffer.data(idx_gr_ff + 268);

    auto grr_z_yyy_zzz = pbuffer.data(idx_gr_ff + 269);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xx, gr_yyy_xxx, gr_yyy_xxy, gr_yyy_xxz, gr_yyy_xy, gr_yyy_xyy, gr_yyy_xyz, gr_yyy_xz, gr_yyy_xzz, gr_yyy_yy, gr_yyy_yyy, gr_yyy_yyz, gr_yyy_yz, gr_yyy_yzz, gr_yyy_zz, gr_yyy_zzz, grr_z_yyy_xxx, grr_z_yyy_xxy, grr_z_yyy_xxz, grr_z_yyy_xyy, grr_z_yyy_xyz, grr_z_yyy_xzz, grr_z_yyy_yyy, grr_z_yyy_yyz, grr_z_yyy_yzz, grr_z_yyy_zzz, ts_yyy_xx, ts_yyy_xxx, ts_yyy_xxy, ts_yyy_xxz, ts_yyy_xy, ts_yyy_xyy, ts_yyy_xyz, ts_yyy_xz, ts_yyy_xzz, ts_yyy_yy, ts_yyy_yyy, ts_yyy_yyz, ts_yyy_yz, ts_yyy_yzz, ts_yyy_zz, ts_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyy_xxx[i] = 2.0 * ts_yyy_xxx[i] * gfe_0 * gc_z[i] + gr_yyy_xxx[i] * gc_z[i];

        grr_z_yyy_xxy[i] = 2.0 * ts_yyy_xxy[i] * gfe_0 * gc_z[i] + gr_yyy_xxy[i] * gc_z[i];

        grr_z_yyy_xxz[i] = 2.0 * ts_yyy_xx[i] * gfe2_0 + gr_yyy_xx[i] * gfe_0 + 2.0 * ts_yyy_xxz[i] * gfe_0 * gc_z[i] + gr_yyy_xxz[i] * gc_z[i];

        grr_z_yyy_xyy[i] = 2.0 * ts_yyy_xyy[i] * gfe_0 * gc_z[i] + gr_yyy_xyy[i] * gc_z[i];

        grr_z_yyy_xyz[i] = 2.0 * ts_yyy_xy[i] * gfe2_0 + gr_yyy_xy[i] * gfe_0 + 2.0 * ts_yyy_xyz[i] * gfe_0 * gc_z[i] + gr_yyy_xyz[i] * gc_z[i];

        grr_z_yyy_xzz[i] = 4.0 * ts_yyy_xz[i] * gfe2_0 + 2.0 * gr_yyy_xz[i] * gfe_0 + 2.0 * ts_yyy_xzz[i] * gfe_0 * gc_z[i] + gr_yyy_xzz[i] * gc_z[i];

        grr_z_yyy_yyy[i] = 2.0 * ts_yyy_yyy[i] * gfe_0 * gc_z[i] + gr_yyy_yyy[i] * gc_z[i];

        grr_z_yyy_yyz[i] = 2.0 * ts_yyy_yy[i] * gfe2_0 + gr_yyy_yy[i] * gfe_0 + 2.0 * ts_yyy_yyz[i] * gfe_0 * gc_z[i] + gr_yyy_yyz[i] * gc_z[i];

        grr_z_yyy_yzz[i] = 4.0 * ts_yyy_yz[i] * gfe2_0 + 2.0 * gr_yyy_yz[i] * gfe_0 + 2.0 * ts_yyy_yzz[i] * gfe_0 * gc_z[i] + gr_yyy_yzz[i] * gc_z[i];

        grr_z_yyy_zzz[i] = 6.0 * ts_yyy_zz[i] * gfe2_0 + 3.0 * gr_yyy_zz[i] * gfe_0 + 2.0 * ts_yyy_zzz[i] * gfe_0 * gc_z[i] + gr_yyy_zzz[i] * gc_z[i];
    }

    // Set up 270-280 components of targeted buffer : FF

    auto grr_z_yyz_xxx = pbuffer.data(idx_gr_ff + 270);

    auto grr_z_yyz_xxy = pbuffer.data(idx_gr_ff + 271);

    auto grr_z_yyz_xxz = pbuffer.data(idx_gr_ff + 272);

    auto grr_z_yyz_xyy = pbuffer.data(idx_gr_ff + 273);

    auto grr_z_yyz_xyz = pbuffer.data(idx_gr_ff + 274);

    auto grr_z_yyz_xzz = pbuffer.data(idx_gr_ff + 275);

    auto grr_z_yyz_yyy = pbuffer.data(idx_gr_ff + 276);

    auto grr_z_yyz_yyz = pbuffer.data(idx_gr_ff + 277);

    auto grr_z_yyz_yzz = pbuffer.data(idx_gr_ff + 278);

    auto grr_z_yyz_zzz = pbuffer.data(idx_gr_ff + 279);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xxx, gr_yy_xxy, gr_yy_xxz, gr_yy_xyy, gr_yy_xyz, gr_yy_xzz, gr_yy_yyy, gr_yy_yyz, gr_yy_yzz, gr_yy_zzz, gr_yyz_xx, gr_yyz_xxx, gr_yyz_xxy, gr_yyz_xxz, gr_yyz_xy, gr_yyz_xyy, gr_yyz_xyz, gr_yyz_xz, gr_yyz_xzz, gr_yyz_yy, gr_yyz_yyy, gr_yyz_yyz, gr_yyz_yz, gr_yyz_yzz, gr_yyz_zz, gr_yyz_zzz, grr_z_yyz_xxx, grr_z_yyz_xxy, grr_z_yyz_xxz, grr_z_yyz_xyy, grr_z_yyz_xyz, grr_z_yyz_xzz, grr_z_yyz_yyy, grr_z_yyz_yyz, grr_z_yyz_yzz, grr_z_yyz_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, ts_yyz_xx, ts_yyz_xxx, ts_yyz_xxy, ts_yyz_xxz, ts_yyz_xy, ts_yyz_xyy, ts_yyz_xyz, ts_yyz_xz, ts_yyz_xzz, ts_yyz_yy, ts_yyz_yyy, ts_yyz_yyz, ts_yyz_yz, ts_yyz_yzz, ts_yyz_zz, ts_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyz_xxx[i] = 2.0 * ts_yy_xxx[i] * gfe2_0 + gr_yy_xxx[i] * gfe_0 + 2.0 * ts_yyz_xxx[i] * gfe_0 * gc_z[i] + gr_yyz_xxx[i] * gc_z[i];

        grr_z_yyz_xxy[i] = 2.0 * ts_yy_xxy[i] * gfe2_0 + gr_yy_xxy[i] * gfe_0 + 2.0 * ts_yyz_xxy[i] * gfe_0 * gc_z[i] + gr_yyz_xxy[i] * gc_z[i];

        grr_z_yyz_xxz[i] = 2.0 * ts_yy_xxz[i] * gfe2_0 + gr_yy_xxz[i] * gfe_0 + 2.0 * ts_yyz_xx[i] * gfe2_0 + gr_yyz_xx[i] * gfe_0 + 2.0 * ts_yyz_xxz[i] * gfe_0 * gc_z[i] + gr_yyz_xxz[i] * gc_z[i];

        grr_z_yyz_xyy[i] = 2.0 * ts_yy_xyy[i] * gfe2_0 + gr_yy_xyy[i] * gfe_0 + 2.0 * ts_yyz_xyy[i] * gfe_0 * gc_z[i] + gr_yyz_xyy[i] * gc_z[i];

        grr_z_yyz_xyz[i] = 2.0 * ts_yy_xyz[i] * gfe2_0 + gr_yy_xyz[i] * gfe_0 + 2.0 * ts_yyz_xy[i] * gfe2_0 + gr_yyz_xy[i] * gfe_0 + 2.0 * ts_yyz_xyz[i] * gfe_0 * gc_z[i] + gr_yyz_xyz[i] * gc_z[i];

        grr_z_yyz_xzz[i] = 2.0 * ts_yy_xzz[i] * gfe2_0 + gr_yy_xzz[i] * gfe_0 + 4.0 * ts_yyz_xz[i] * gfe2_0 + 2.0 * gr_yyz_xz[i] * gfe_0 + 2.0 * ts_yyz_xzz[i] * gfe_0 * gc_z[i] + gr_yyz_xzz[i] * gc_z[i];

        grr_z_yyz_yyy[i] = 2.0 * ts_yy_yyy[i] * gfe2_0 + gr_yy_yyy[i] * gfe_0 + 2.0 * ts_yyz_yyy[i] * gfe_0 * gc_z[i] + gr_yyz_yyy[i] * gc_z[i];

        grr_z_yyz_yyz[i] = 2.0 * ts_yy_yyz[i] * gfe2_0 + gr_yy_yyz[i] * gfe_0 + 2.0 * ts_yyz_yy[i] * gfe2_0 + gr_yyz_yy[i] * gfe_0 + 2.0 * ts_yyz_yyz[i] * gfe_0 * gc_z[i] + gr_yyz_yyz[i] * gc_z[i];

        grr_z_yyz_yzz[i] = 2.0 * ts_yy_yzz[i] * gfe2_0 + gr_yy_yzz[i] * gfe_0 + 4.0 * ts_yyz_yz[i] * gfe2_0 + 2.0 * gr_yyz_yz[i] * gfe_0 + 2.0 * ts_yyz_yzz[i] * gfe_0 * gc_z[i] + gr_yyz_yzz[i] * gc_z[i];

        grr_z_yyz_zzz[i] = 2.0 * ts_yy_zzz[i] * gfe2_0 + gr_yy_zzz[i] * gfe_0 + 6.0 * ts_yyz_zz[i] * gfe2_0 + 3.0 * gr_yyz_zz[i] * gfe_0 + 2.0 * ts_yyz_zzz[i] * gfe_0 * gc_z[i] + gr_yyz_zzz[i] * gc_z[i];
    }

    // Set up 280-290 components of targeted buffer : FF

    auto grr_z_yzz_xxx = pbuffer.data(idx_gr_ff + 280);

    auto grr_z_yzz_xxy = pbuffer.data(idx_gr_ff + 281);

    auto grr_z_yzz_xxz = pbuffer.data(idx_gr_ff + 282);

    auto grr_z_yzz_xyy = pbuffer.data(idx_gr_ff + 283);

    auto grr_z_yzz_xyz = pbuffer.data(idx_gr_ff + 284);

    auto grr_z_yzz_xzz = pbuffer.data(idx_gr_ff + 285);

    auto grr_z_yzz_yyy = pbuffer.data(idx_gr_ff + 286);

    auto grr_z_yzz_yyz = pbuffer.data(idx_gr_ff + 287);

    auto grr_z_yzz_yzz = pbuffer.data(idx_gr_ff + 288);

    auto grr_z_yzz_zzz = pbuffer.data(idx_gr_ff + 289);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xxx, gr_yz_xxy, gr_yz_xxz, gr_yz_xyy, gr_yz_xyz, gr_yz_xzz, gr_yz_yyy, gr_yz_yyz, gr_yz_yzz, gr_yz_zzz, gr_yzz_xx, gr_yzz_xxx, gr_yzz_xxy, gr_yzz_xxz, gr_yzz_xy, gr_yzz_xyy, gr_yzz_xyz, gr_yzz_xz, gr_yzz_xzz, gr_yzz_yy, gr_yzz_yyy, gr_yzz_yyz, gr_yzz_yz, gr_yzz_yzz, gr_yzz_zz, gr_yzz_zzz, grr_z_yzz_xxx, grr_z_yzz_xxy, grr_z_yzz_xxz, grr_z_yzz_xyy, grr_z_yzz_xyz, grr_z_yzz_xzz, grr_z_yzz_yyy, grr_z_yzz_yyz, grr_z_yzz_yzz, grr_z_yzz_zzz, ts_yz_xxx, ts_yz_xxy, ts_yz_xxz, ts_yz_xyy, ts_yz_xyz, ts_yz_xzz, ts_yz_yyy, ts_yz_yyz, ts_yz_yzz, ts_yz_zzz, ts_yzz_xx, ts_yzz_xxx, ts_yzz_xxy, ts_yzz_xxz, ts_yzz_xy, ts_yzz_xyy, ts_yzz_xyz, ts_yzz_xz, ts_yzz_xzz, ts_yzz_yy, ts_yzz_yyy, ts_yzz_yyz, ts_yzz_yz, ts_yzz_yzz, ts_yzz_zz, ts_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yzz_xxx[i] = 4.0 * ts_yz_xxx[i] * gfe2_0 + 2.0 * gr_yz_xxx[i] * gfe_0 + 2.0 * ts_yzz_xxx[i] * gfe_0 * gc_z[i] + gr_yzz_xxx[i] * gc_z[i];

        grr_z_yzz_xxy[i] = 4.0 * ts_yz_xxy[i] * gfe2_0 + 2.0 * gr_yz_xxy[i] * gfe_0 + 2.0 * ts_yzz_xxy[i] * gfe_0 * gc_z[i] + gr_yzz_xxy[i] * gc_z[i];

        grr_z_yzz_xxz[i] = 4.0 * ts_yz_xxz[i] * gfe2_0 + 2.0 * gr_yz_xxz[i] * gfe_0 + 2.0 * ts_yzz_xx[i] * gfe2_0 + gr_yzz_xx[i] * gfe_0 + 2.0 * ts_yzz_xxz[i] * gfe_0 * gc_z[i] + gr_yzz_xxz[i] * gc_z[i];

        grr_z_yzz_xyy[i] = 4.0 * ts_yz_xyy[i] * gfe2_0 + 2.0 * gr_yz_xyy[i] * gfe_0 + 2.0 * ts_yzz_xyy[i] * gfe_0 * gc_z[i] + gr_yzz_xyy[i] * gc_z[i];

        grr_z_yzz_xyz[i] = 4.0 * ts_yz_xyz[i] * gfe2_0 + 2.0 * gr_yz_xyz[i] * gfe_0 + 2.0 * ts_yzz_xy[i] * gfe2_0 + gr_yzz_xy[i] * gfe_0 + 2.0 * ts_yzz_xyz[i] * gfe_0 * gc_z[i] + gr_yzz_xyz[i] * gc_z[i];

        grr_z_yzz_xzz[i] = 4.0 * ts_yz_xzz[i] * gfe2_0 + 2.0 * gr_yz_xzz[i] * gfe_0 + 4.0 * ts_yzz_xz[i] * gfe2_0 + 2.0 * gr_yzz_xz[i] * gfe_0 + 2.0 * ts_yzz_xzz[i] * gfe_0 * gc_z[i] + gr_yzz_xzz[i] * gc_z[i];

        grr_z_yzz_yyy[i] = 4.0 * ts_yz_yyy[i] * gfe2_0 + 2.0 * gr_yz_yyy[i] * gfe_0 + 2.0 * ts_yzz_yyy[i] * gfe_0 * gc_z[i] + gr_yzz_yyy[i] * gc_z[i];

        grr_z_yzz_yyz[i] = 4.0 * ts_yz_yyz[i] * gfe2_0 + 2.0 * gr_yz_yyz[i] * gfe_0 + 2.0 * ts_yzz_yy[i] * gfe2_0 + gr_yzz_yy[i] * gfe_0 + 2.0 * ts_yzz_yyz[i] * gfe_0 * gc_z[i] + gr_yzz_yyz[i] * gc_z[i];

        grr_z_yzz_yzz[i] = 4.0 * ts_yz_yzz[i] * gfe2_0 + 2.0 * gr_yz_yzz[i] * gfe_0 + 4.0 * ts_yzz_yz[i] * gfe2_0 + 2.0 * gr_yzz_yz[i] * gfe_0 + 2.0 * ts_yzz_yzz[i] * gfe_0 * gc_z[i] + gr_yzz_yzz[i] * gc_z[i];

        grr_z_yzz_zzz[i] = 4.0 * ts_yz_zzz[i] * gfe2_0 + 2.0 * gr_yz_zzz[i] * gfe_0 + 6.0 * ts_yzz_zz[i] * gfe2_0 + 3.0 * gr_yzz_zz[i] * gfe_0 + 2.0 * ts_yzz_zzz[i] * gfe_0 * gc_z[i] + gr_yzz_zzz[i] * gc_z[i];
    }

    // Set up 290-300 components of targeted buffer : FF

    auto grr_z_zzz_xxx = pbuffer.data(idx_gr_ff + 290);

    auto grr_z_zzz_xxy = pbuffer.data(idx_gr_ff + 291);

    auto grr_z_zzz_xxz = pbuffer.data(idx_gr_ff + 292);

    auto grr_z_zzz_xyy = pbuffer.data(idx_gr_ff + 293);

    auto grr_z_zzz_xyz = pbuffer.data(idx_gr_ff + 294);

    auto grr_z_zzz_xzz = pbuffer.data(idx_gr_ff + 295);

    auto grr_z_zzz_yyy = pbuffer.data(idx_gr_ff + 296);

    auto grr_z_zzz_yyz = pbuffer.data(idx_gr_ff + 297);

    auto grr_z_zzz_yzz = pbuffer.data(idx_gr_ff + 298);

    auto grr_z_zzz_zzz = pbuffer.data(idx_gr_ff + 299);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xxx, gr_zz_xxy, gr_zz_xxz, gr_zz_xyy, gr_zz_xyz, gr_zz_xzz, gr_zz_yyy, gr_zz_yyz, gr_zz_yzz, gr_zz_zzz, gr_zzz_xx, gr_zzz_xxx, gr_zzz_xxy, gr_zzz_xxz, gr_zzz_xy, gr_zzz_xyy, gr_zzz_xyz, gr_zzz_xz, gr_zzz_xzz, gr_zzz_yy, gr_zzz_yyy, gr_zzz_yyz, gr_zzz_yz, gr_zzz_yzz, gr_zzz_zz, gr_zzz_zzz, grr_z_zzz_xxx, grr_z_zzz_xxy, grr_z_zzz_xxz, grr_z_zzz_xyy, grr_z_zzz_xyz, grr_z_zzz_xzz, grr_z_zzz_yyy, grr_z_zzz_yyz, grr_z_zzz_yzz, grr_z_zzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, ts_zzz_xx, ts_zzz_xxx, ts_zzz_xxy, ts_zzz_xxz, ts_zzz_xy, ts_zzz_xyy, ts_zzz_xyz, ts_zzz_xz, ts_zzz_xzz, ts_zzz_yy, ts_zzz_yyy, ts_zzz_yyz, ts_zzz_yz, ts_zzz_yzz, ts_zzz_zz, ts_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zzz_xxx[i] = 6.0 * ts_zz_xxx[i] * gfe2_0 + 3.0 * gr_zz_xxx[i] * gfe_0 + 2.0 * ts_zzz_xxx[i] * gfe_0 * gc_z[i] + gr_zzz_xxx[i] * gc_z[i];

        grr_z_zzz_xxy[i] = 6.0 * ts_zz_xxy[i] * gfe2_0 + 3.0 * gr_zz_xxy[i] * gfe_0 + 2.0 * ts_zzz_xxy[i] * gfe_0 * gc_z[i] + gr_zzz_xxy[i] * gc_z[i];

        grr_z_zzz_xxz[i] = 6.0 * ts_zz_xxz[i] * gfe2_0 + 3.0 * gr_zz_xxz[i] * gfe_0 + 2.0 * ts_zzz_xx[i] * gfe2_0 + gr_zzz_xx[i] * gfe_0 + 2.0 * ts_zzz_xxz[i] * gfe_0 * gc_z[i] + gr_zzz_xxz[i] * gc_z[i];

        grr_z_zzz_xyy[i] = 6.0 * ts_zz_xyy[i] * gfe2_0 + 3.0 * gr_zz_xyy[i] * gfe_0 + 2.0 * ts_zzz_xyy[i] * gfe_0 * gc_z[i] + gr_zzz_xyy[i] * gc_z[i];

        grr_z_zzz_xyz[i] = 6.0 * ts_zz_xyz[i] * gfe2_0 + 3.0 * gr_zz_xyz[i] * gfe_0 + 2.0 * ts_zzz_xy[i] * gfe2_0 + gr_zzz_xy[i] * gfe_0 + 2.0 * ts_zzz_xyz[i] * gfe_0 * gc_z[i] + gr_zzz_xyz[i] * gc_z[i];

        grr_z_zzz_xzz[i] = 6.0 * ts_zz_xzz[i] * gfe2_0 + 3.0 * gr_zz_xzz[i] * gfe_0 + 4.0 * ts_zzz_xz[i] * gfe2_0 + 2.0 * gr_zzz_xz[i] * gfe_0 + 2.0 * ts_zzz_xzz[i] * gfe_0 * gc_z[i] + gr_zzz_xzz[i] * gc_z[i];

        grr_z_zzz_yyy[i] = 6.0 * ts_zz_yyy[i] * gfe2_0 + 3.0 * gr_zz_yyy[i] * gfe_0 + 2.0 * ts_zzz_yyy[i] * gfe_0 * gc_z[i] + gr_zzz_yyy[i] * gc_z[i];

        grr_z_zzz_yyz[i] = 6.0 * ts_zz_yyz[i] * gfe2_0 + 3.0 * gr_zz_yyz[i] * gfe_0 + 2.0 * ts_zzz_yy[i] * gfe2_0 + gr_zzz_yy[i] * gfe_0 + 2.0 * ts_zzz_yyz[i] * gfe_0 * gc_z[i] + gr_zzz_yyz[i] * gc_z[i];

        grr_z_zzz_yzz[i] = 6.0 * ts_zz_yzz[i] * gfe2_0 + 3.0 * gr_zz_yzz[i] * gfe_0 + 4.0 * ts_zzz_yz[i] * gfe2_0 + 2.0 * gr_zzz_yz[i] * gfe_0 + 2.0 * ts_zzz_yzz[i] * gfe_0 * gc_z[i] + gr_zzz_yzz[i] * gc_z[i];

        grr_z_zzz_zzz[i] = 6.0 * ts_zz_zzz[i] * gfe2_0 + 3.0 * gr_zz_zzz[i] * gfe_0 + 6.0 * ts_zzz_zz[i] * gfe2_0 + 3.0 * gr_zzz_zz[i] * gfe_0 + 2.0 * ts_zzz_zzz[i] * gfe_0 * gc_z[i] + gr_zzz_zzz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

