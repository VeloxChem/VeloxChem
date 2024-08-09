#include "ElectricDipoleMomentumPrimRecFF.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_ff(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_ff,
                                      const size_t idx_dip_pf,
                                      const size_t idx_dip_dd,
                                      const size_t idx_ovl_df,
                                      const size_t idx_dip_df,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : PF

    auto tr_x_x_xxx = pbuffer.data(idx_dip_pf);

    auto tr_x_x_xxy = pbuffer.data(idx_dip_pf + 1);

    auto tr_x_x_xxz = pbuffer.data(idx_dip_pf + 2);

    auto tr_x_x_xyy = pbuffer.data(idx_dip_pf + 3);

    auto tr_x_x_xyz = pbuffer.data(idx_dip_pf + 4);

    auto tr_x_x_xzz = pbuffer.data(idx_dip_pf + 5);

    auto tr_x_x_yyy = pbuffer.data(idx_dip_pf + 6);

    auto tr_x_x_yyz = pbuffer.data(idx_dip_pf + 7);

    auto tr_x_x_yzz = pbuffer.data(idx_dip_pf + 8);

    auto tr_x_x_zzz = pbuffer.data(idx_dip_pf + 9);

    auto tr_x_y_xxx = pbuffer.data(idx_dip_pf + 10);

    auto tr_x_y_xxy = pbuffer.data(idx_dip_pf + 11);

    auto tr_x_y_xxz = pbuffer.data(idx_dip_pf + 12);

    auto tr_x_y_xyy = pbuffer.data(idx_dip_pf + 13);

    auto tr_x_y_xyz = pbuffer.data(idx_dip_pf + 14);

    auto tr_x_y_xzz = pbuffer.data(idx_dip_pf + 15);

    auto tr_x_y_yyy = pbuffer.data(idx_dip_pf + 16);

    auto tr_x_y_yyz = pbuffer.data(idx_dip_pf + 17);

    auto tr_x_y_yzz = pbuffer.data(idx_dip_pf + 18);

    auto tr_x_y_zzz = pbuffer.data(idx_dip_pf + 19);

    auto tr_x_z_xxx = pbuffer.data(idx_dip_pf + 20);

    auto tr_x_z_xxy = pbuffer.data(idx_dip_pf + 21);

    auto tr_x_z_xxz = pbuffer.data(idx_dip_pf + 22);

    auto tr_x_z_xyy = pbuffer.data(idx_dip_pf + 23);

    auto tr_x_z_xyz = pbuffer.data(idx_dip_pf + 24);

    auto tr_x_z_xzz = pbuffer.data(idx_dip_pf + 25);

    auto tr_x_z_yyy = pbuffer.data(idx_dip_pf + 26);

    auto tr_x_z_yyz = pbuffer.data(idx_dip_pf + 27);

    auto tr_x_z_yzz = pbuffer.data(idx_dip_pf + 28);

    auto tr_x_z_zzz = pbuffer.data(idx_dip_pf + 29);

    auto tr_y_x_xxx = pbuffer.data(idx_dip_pf + 30);

    auto tr_y_x_xxy = pbuffer.data(idx_dip_pf + 31);

    auto tr_y_x_xxz = pbuffer.data(idx_dip_pf + 32);

    auto tr_y_x_xyy = pbuffer.data(idx_dip_pf + 33);

    auto tr_y_x_xyz = pbuffer.data(idx_dip_pf + 34);

    auto tr_y_x_xzz = pbuffer.data(idx_dip_pf + 35);

    auto tr_y_x_yyy = pbuffer.data(idx_dip_pf + 36);

    auto tr_y_x_yyz = pbuffer.data(idx_dip_pf + 37);

    auto tr_y_x_yzz = pbuffer.data(idx_dip_pf + 38);

    auto tr_y_x_zzz = pbuffer.data(idx_dip_pf + 39);

    auto tr_y_y_xxx = pbuffer.data(idx_dip_pf + 40);

    auto tr_y_y_xxy = pbuffer.data(idx_dip_pf + 41);

    auto tr_y_y_xxz = pbuffer.data(idx_dip_pf + 42);

    auto tr_y_y_xyy = pbuffer.data(idx_dip_pf + 43);

    auto tr_y_y_xyz = pbuffer.data(idx_dip_pf + 44);

    auto tr_y_y_xzz = pbuffer.data(idx_dip_pf + 45);

    auto tr_y_y_yyy = pbuffer.data(idx_dip_pf + 46);

    auto tr_y_y_yyz = pbuffer.data(idx_dip_pf + 47);

    auto tr_y_y_yzz = pbuffer.data(idx_dip_pf + 48);

    auto tr_y_y_zzz = pbuffer.data(idx_dip_pf + 49);

    auto tr_y_z_xxx = pbuffer.data(idx_dip_pf + 50);

    auto tr_y_z_xxy = pbuffer.data(idx_dip_pf + 51);

    auto tr_y_z_xxz = pbuffer.data(idx_dip_pf + 52);

    auto tr_y_z_xyy = pbuffer.data(idx_dip_pf + 53);

    auto tr_y_z_xyz = pbuffer.data(idx_dip_pf + 54);

    auto tr_y_z_xzz = pbuffer.data(idx_dip_pf + 55);

    auto tr_y_z_yyy = pbuffer.data(idx_dip_pf + 56);

    auto tr_y_z_yyz = pbuffer.data(idx_dip_pf + 57);

    auto tr_y_z_yzz = pbuffer.data(idx_dip_pf + 58);

    auto tr_y_z_zzz = pbuffer.data(idx_dip_pf + 59);

    auto tr_z_x_xxx = pbuffer.data(idx_dip_pf + 60);

    auto tr_z_x_xxy = pbuffer.data(idx_dip_pf + 61);

    auto tr_z_x_xxz = pbuffer.data(idx_dip_pf + 62);

    auto tr_z_x_xyy = pbuffer.data(idx_dip_pf + 63);

    auto tr_z_x_xyz = pbuffer.data(idx_dip_pf + 64);

    auto tr_z_x_xzz = pbuffer.data(idx_dip_pf + 65);

    auto tr_z_x_yyy = pbuffer.data(idx_dip_pf + 66);

    auto tr_z_x_yyz = pbuffer.data(idx_dip_pf + 67);

    auto tr_z_x_yzz = pbuffer.data(idx_dip_pf + 68);

    auto tr_z_x_zzz = pbuffer.data(idx_dip_pf + 69);

    auto tr_z_y_xxx = pbuffer.data(idx_dip_pf + 70);

    auto tr_z_y_xxy = pbuffer.data(idx_dip_pf + 71);

    auto tr_z_y_xxz = pbuffer.data(idx_dip_pf + 72);

    auto tr_z_y_xyy = pbuffer.data(idx_dip_pf + 73);

    auto tr_z_y_xyz = pbuffer.data(idx_dip_pf + 74);

    auto tr_z_y_xzz = pbuffer.data(idx_dip_pf + 75);

    auto tr_z_y_yyy = pbuffer.data(idx_dip_pf + 76);

    auto tr_z_y_yyz = pbuffer.data(idx_dip_pf + 77);

    auto tr_z_y_yzz = pbuffer.data(idx_dip_pf + 78);

    auto tr_z_y_zzz = pbuffer.data(idx_dip_pf + 79);

    auto tr_z_z_xxx = pbuffer.data(idx_dip_pf + 80);

    auto tr_z_z_xxy = pbuffer.data(idx_dip_pf + 81);

    auto tr_z_z_xxz = pbuffer.data(idx_dip_pf + 82);

    auto tr_z_z_xyy = pbuffer.data(idx_dip_pf + 83);

    auto tr_z_z_xyz = pbuffer.data(idx_dip_pf + 84);

    auto tr_z_z_xzz = pbuffer.data(idx_dip_pf + 85);

    auto tr_z_z_yyy = pbuffer.data(idx_dip_pf + 86);

    auto tr_z_z_yyz = pbuffer.data(idx_dip_pf + 87);

    auto tr_z_z_yzz = pbuffer.data(idx_dip_pf + 88);

    auto tr_z_z_zzz = pbuffer.data(idx_dip_pf + 89);

    // Set up components of auxiliary buffer : DD

    auto tr_x_xx_xx = pbuffer.data(idx_dip_dd);

    auto tr_x_xx_xy = pbuffer.data(idx_dip_dd + 1);

    auto tr_x_xx_xz = pbuffer.data(idx_dip_dd + 2);

    auto tr_x_xx_yy = pbuffer.data(idx_dip_dd + 3);

    auto tr_x_xx_yz = pbuffer.data(idx_dip_dd + 4);

    auto tr_x_xx_zz = pbuffer.data(idx_dip_dd + 5);

    auto tr_x_xz_xz = pbuffer.data(idx_dip_dd + 14);

    auto tr_x_yy_xx = pbuffer.data(idx_dip_dd + 18);

    auto tr_x_yy_xy = pbuffer.data(idx_dip_dd + 19);

    auto tr_x_yy_xz = pbuffer.data(idx_dip_dd + 20);

    auto tr_x_yy_yy = pbuffer.data(idx_dip_dd + 21);

    auto tr_x_yy_yz = pbuffer.data(idx_dip_dd + 22);

    auto tr_x_yy_zz = pbuffer.data(idx_dip_dd + 23);

    auto tr_x_zz_xx = pbuffer.data(idx_dip_dd + 30);

    auto tr_x_zz_xy = pbuffer.data(idx_dip_dd + 31);

    auto tr_x_zz_xz = pbuffer.data(idx_dip_dd + 32);

    auto tr_x_zz_yy = pbuffer.data(idx_dip_dd + 33);

    auto tr_x_zz_yz = pbuffer.data(idx_dip_dd + 34);

    auto tr_x_zz_zz = pbuffer.data(idx_dip_dd + 35);

    auto tr_y_xx_xx = pbuffer.data(idx_dip_dd + 36);

    auto tr_y_xx_xy = pbuffer.data(idx_dip_dd + 37);

    auto tr_y_xx_xz = pbuffer.data(idx_dip_dd + 38);

    auto tr_y_xx_yy = pbuffer.data(idx_dip_dd + 39);

    auto tr_y_xx_yz = pbuffer.data(idx_dip_dd + 40);

    auto tr_y_xx_zz = pbuffer.data(idx_dip_dd + 41);

    auto tr_y_xy_xy = pbuffer.data(idx_dip_dd + 43);

    auto tr_y_xy_yy = pbuffer.data(idx_dip_dd + 45);

    auto tr_y_xy_yz = pbuffer.data(idx_dip_dd + 46);

    auto tr_y_yy_xx = pbuffer.data(idx_dip_dd + 54);

    auto tr_y_yy_xy = pbuffer.data(idx_dip_dd + 55);

    auto tr_y_yy_xz = pbuffer.data(idx_dip_dd + 56);

    auto tr_y_yy_yy = pbuffer.data(idx_dip_dd + 57);

    auto tr_y_yy_yz = pbuffer.data(idx_dip_dd + 58);

    auto tr_y_yy_zz = pbuffer.data(idx_dip_dd + 59);

    auto tr_y_yz_xz = pbuffer.data(idx_dip_dd + 62);

    auto tr_y_yz_yz = pbuffer.data(idx_dip_dd + 64);

    auto tr_y_yz_zz = pbuffer.data(idx_dip_dd + 65);

    auto tr_y_zz_xx = pbuffer.data(idx_dip_dd + 66);

    auto tr_y_zz_xy = pbuffer.data(idx_dip_dd + 67);

    auto tr_y_zz_xz = pbuffer.data(idx_dip_dd + 68);

    auto tr_y_zz_yy = pbuffer.data(idx_dip_dd + 69);

    auto tr_y_zz_yz = pbuffer.data(idx_dip_dd + 70);

    auto tr_y_zz_zz = pbuffer.data(idx_dip_dd + 71);

    auto tr_z_xx_xx = pbuffer.data(idx_dip_dd + 72);

    auto tr_z_xx_xy = pbuffer.data(idx_dip_dd + 73);

    auto tr_z_xx_xz = pbuffer.data(idx_dip_dd + 74);

    auto tr_z_xx_yy = pbuffer.data(idx_dip_dd + 75);

    auto tr_z_xx_yz = pbuffer.data(idx_dip_dd + 76);

    auto tr_z_xx_zz = pbuffer.data(idx_dip_dd + 77);

    auto tr_z_xz_xz = pbuffer.data(idx_dip_dd + 86);

    auto tr_z_xz_yz = pbuffer.data(idx_dip_dd + 88);

    auto tr_z_xz_zz = pbuffer.data(idx_dip_dd + 89);

    auto tr_z_yy_xx = pbuffer.data(idx_dip_dd + 90);

    auto tr_z_yy_xy = pbuffer.data(idx_dip_dd + 91);

    auto tr_z_yy_xz = pbuffer.data(idx_dip_dd + 92);

    auto tr_z_yy_yy = pbuffer.data(idx_dip_dd + 93);

    auto tr_z_yy_yz = pbuffer.data(idx_dip_dd + 94);

    auto tr_z_yy_zz = pbuffer.data(idx_dip_dd + 95);

    auto tr_z_yz_xy = pbuffer.data(idx_dip_dd + 97);

    auto tr_z_yz_xz = pbuffer.data(idx_dip_dd + 98);

    auto tr_z_yz_yy = pbuffer.data(idx_dip_dd + 99);

    auto tr_z_yz_yz = pbuffer.data(idx_dip_dd + 100);

    auto tr_z_yz_zz = pbuffer.data(idx_dip_dd + 101);

    auto tr_z_zz_xx = pbuffer.data(idx_dip_dd + 102);

    auto tr_z_zz_xy = pbuffer.data(idx_dip_dd + 103);

    auto tr_z_zz_xz = pbuffer.data(idx_dip_dd + 104);

    auto tr_z_zz_yy = pbuffer.data(idx_dip_dd + 105);

    auto tr_z_zz_yz = pbuffer.data(idx_dip_dd + 106);

    auto tr_z_zz_zz = pbuffer.data(idx_dip_dd + 107);

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_ovl_df);

    auto ts_xx_xxy = pbuffer.data(idx_ovl_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_ovl_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_ovl_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_ovl_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_ovl_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_ovl_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_ovl_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_ovl_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_ovl_df + 9);

    auto ts_yy_xxx = pbuffer.data(idx_ovl_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_ovl_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_ovl_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_ovl_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_ovl_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_ovl_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_ovl_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_ovl_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_ovl_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_ovl_df + 39);

    auto ts_yz_yyz = pbuffer.data(idx_ovl_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_ovl_df + 48);

    auto ts_zz_xxx = pbuffer.data(idx_ovl_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_ovl_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_ovl_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_ovl_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_ovl_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_ovl_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_ovl_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_ovl_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_ovl_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_ovl_df + 59);

    // Set up components of auxiliary buffer : DF

    auto tr_x_xx_xxx = pbuffer.data(idx_dip_df);

    auto tr_x_xx_xxy = pbuffer.data(idx_dip_df + 1);

    auto tr_x_xx_xxz = pbuffer.data(idx_dip_df + 2);

    auto tr_x_xx_xyy = pbuffer.data(idx_dip_df + 3);

    auto tr_x_xx_xyz = pbuffer.data(idx_dip_df + 4);

    auto tr_x_xx_xzz = pbuffer.data(idx_dip_df + 5);

    auto tr_x_xx_yyy = pbuffer.data(idx_dip_df + 6);

    auto tr_x_xx_yyz = pbuffer.data(idx_dip_df + 7);

    auto tr_x_xx_yzz = pbuffer.data(idx_dip_df + 8);

    auto tr_x_xx_zzz = pbuffer.data(idx_dip_df + 9);

    auto tr_x_xy_xxx = pbuffer.data(idx_dip_df + 10);

    auto tr_x_xy_xxy = pbuffer.data(idx_dip_df + 11);

    auto tr_x_xy_xxz = pbuffer.data(idx_dip_df + 12);

    auto tr_x_xy_xyy = pbuffer.data(idx_dip_df + 13);

    auto tr_x_xy_xzz = pbuffer.data(idx_dip_df + 15);

    auto tr_x_xy_yyy = pbuffer.data(idx_dip_df + 16);

    auto tr_x_xz_xxx = pbuffer.data(idx_dip_df + 20);

    auto tr_x_xz_xxy = pbuffer.data(idx_dip_df + 21);

    auto tr_x_xz_xxz = pbuffer.data(idx_dip_df + 22);

    auto tr_x_xz_xyy = pbuffer.data(idx_dip_df + 23);

    auto tr_x_xz_xyz = pbuffer.data(idx_dip_df + 24);

    auto tr_x_xz_xzz = pbuffer.data(idx_dip_df + 25);

    auto tr_x_xz_zzz = pbuffer.data(idx_dip_df + 29);

    auto tr_x_yy_xxx = pbuffer.data(idx_dip_df + 30);

    auto tr_x_yy_xxy = pbuffer.data(idx_dip_df + 31);

    auto tr_x_yy_xxz = pbuffer.data(idx_dip_df + 32);

    auto tr_x_yy_xyy = pbuffer.data(idx_dip_df + 33);

    auto tr_x_yy_xyz = pbuffer.data(idx_dip_df + 34);

    auto tr_x_yy_xzz = pbuffer.data(idx_dip_df + 35);

    auto tr_x_yy_yyy = pbuffer.data(idx_dip_df + 36);

    auto tr_x_yy_yyz = pbuffer.data(idx_dip_df + 37);

    auto tr_x_yy_yzz = pbuffer.data(idx_dip_df + 38);

    auto tr_x_yy_zzz = pbuffer.data(idx_dip_df + 39);

    auto tr_x_yz_xxz = pbuffer.data(idx_dip_df + 42);

    auto tr_x_yz_xzz = pbuffer.data(idx_dip_df + 45);

    auto tr_x_yz_yyz = pbuffer.data(idx_dip_df + 47);

    auto tr_x_yz_yzz = pbuffer.data(idx_dip_df + 48);

    auto tr_x_yz_zzz = pbuffer.data(idx_dip_df + 49);

    auto tr_x_zz_xxx = pbuffer.data(idx_dip_df + 50);

    auto tr_x_zz_xxy = pbuffer.data(idx_dip_df + 51);

    auto tr_x_zz_xxz = pbuffer.data(idx_dip_df + 52);

    auto tr_x_zz_xyy = pbuffer.data(idx_dip_df + 53);

    auto tr_x_zz_xyz = pbuffer.data(idx_dip_df + 54);

    auto tr_x_zz_xzz = pbuffer.data(idx_dip_df + 55);

    auto tr_x_zz_yyy = pbuffer.data(idx_dip_df + 56);

    auto tr_x_zz_yyz = pbuffer.data(idx_dip_df + 57);

    auto tr_x_zz_yzz = pbuffer.data(idx_dip_df + 58);

    auto tr_x_zz_zzz = pbuffer.data(idx_dip_df + 59);

    auto tr_y_xx_xxx = pbuffer.data(idx_dip_df + 60);

    auto tr_y_xx_xxy = pbuffer.data(idx_dip_df + 61);

    auto tr_y_xx_xxz = pbuffer.data(idx_dip_df + 62);

    auto tr_y_xx_xyy = pbuffer.data(idx_dip_df + 63);

    auto tr_y_xx_xyz = pbuffer.data(idx_dip_df + 64);

    auto tr_y_xx_xzz = pbuffer.data(idx_dip_df + 65);

    auto tr_y_xx_yyy = pbuffer.data(idx_dip_df + 66);

    auto tr_y_xx_yyz = pbuffer.data(idx_dip_df + 67);

    auto tr_y_xx_yzz = pbuffer.data(idx_dip_df + 68);

    auto tr_y_xx_zzz = pbuffer.data(idx_dip_df + 69);

    auto tr_y_xy_xxx = pbuffer.data(idx_dip_df + 70);

    auto tr_y_xy_xxy = pbuffer.data(idx_dip_df + 71);

    auto tr_y_xy_xyy = pbuffer.data(idx_dip_df + 73);

    auto tr_y_xy_xyz = pbuffer.data(idx_dip_df + 74);

    auto tr_y_xy_yyy = pbuffer.data(idx_dip_df + 76);

    auto tr_y_xy_yyz = pbuffer.data(idx_dip_df + 77);

    auto tr_y_xy_yzz = pbuffer.data(idx_dip_df + 78);

    auto tr_y_xy_zzz = pbuffer.data(idx_dip_df + 79);

    auto tr_y_xz_yyz = pbuffer.data(idx_dip_df + 87);

    auto tr_y_xz_yzz = pbuffer.data(idx_dip_df + 88);

    auto tr_y_xz_zzz = pbuffer.data(idx_dip_df + 89);

    auto tr_y_yy_xxx = pbuffer.data(idx_dip_df + 90);

    auto tr_y_yy_xxy = pbuffer.data(idx_dip_df + 91);

    auto tr_y_yy_xxz = pbuffer.data(idx_dip_df + 92);

    auto tr_y_yy_xyy = pbuffer.data(idx_dip_df + 93);

    auto tr_y_yy_xyz = pbuffer.data(idx_dip_df + 94);

    auto tr_y_yy_xzz = pbuffer.data(idx_dip_df + 95);

    auto tr_y_yy_yyy = pbuffer.data(idx_dip_df + 96);

    auto tr_y_yy_yyz = pbuffer.data(idx_dip_df + 97);

    auto tr_y_yy_yzz = pbuffer.data(idx_dip_df + 98);

    auto tr_y_yy_zzz = pbuffer.data(idx_dip_df + 99);

    auto tr_y_yz_xxy = pbuffer.data(idx_dip_df + 101);

    auto tr_y_yz_xxz = pbuffer.data(idx_dip_df + 102);

    auto tr_y_yz_xyy = pbuffer.data(idx_dip_df + 103);

    auto tr_y_yz_xyz = pbuffer.data(idx_dip_df + 104);

    auto tr_y_yz_xzz = pbuffer.data(idx_dip_df + 105);

    auto tr_y_yz_yyy = pbuffer.data(idx_dip_df + 106);

    auto tr_y_yz_yyz = pbuffer.data(idx_dip_df + 107);

    auto tr_y_yz_yzz = pbuffer.data(idx_dip_df + 108);

    auto tr_y_yz_zzz = pbuffer.data(idx_dip_df + 109);

    auto tr_y_zz_xxx = pbuffer.data(idx_dip_df + 110);

    auto tr_y_zz_xxy = pbuffer.data(idx_dip_df + 111);

    auto tr_y_zz_xxz = pbuffer.data(idx_dip_df + 112);

    auto tr_y_zz_xyy = pbuffer.data(idx_dip_df + 113);

    auto tr_y_zz_xyz = pbuffer.data(idx_dip_df + 114);

    auto tr_y_zz_xzz = pbuffer.data(idx_dip_df + 115);

    auto tr_y_zz_yyy = pbuffer.data(idx_dip_df + 116);

    auto tr_y_zz_yyz = pbuffer.data(idx_dip_df + 117);

    auto tr_y_zz_yzz = pbuffer.data(idx_dip_df + 118);

    auto tr_y_zz_zzz = pbuffer.data(idx_dip_df + 119);

    auto tr_z_xx_xxx = pbuffer.data(idx_dip_df + 120);

    auto tr_z_xx_xxy = pbuffer.data(idx_dip_df + 121);

    auto tr_z_xx_xxz = pbuffer.data(idx_dip_df + 122);

    auto tr_z_xx_xyy = pbuffer.data(idx_dip_df + 123);

    auto tr_z_xx_xyz = pbuffer.data(idx_dip_df + 124);

    auto tr_z_xx_xzz = pbuffer.data(idx_dip_df + 125);

    auto tr_z_xx_yyy = pbuffer.data(idx_dip_df + 126);

    auto tr_z_xx_yyz = pbuffer.data(idx_dip_df + 127);

    auto tr_z_xx_yzz = pbuffer.data(idx_dip_df + 128);

    auto tr_z_xx_zzz = pbuffer.data(idx_dip_df + 129);

    auto tr_z_xy_yyy = pbuffer.data(idx_dip_df + 136);

    auto tr_z_xy_yyz = pbuffer.data(idx_dip_df + 137);

    auto tr_z_xy_yzz = pbuffer.data(idx_dip_df + 138);

    auto tr_z_xz_xxx = pbuffer.data(idx_dip_df + 140);

    auto tr_z_xz_xxz = pbuffer.data(idx_dip_df + 142);

    auto tr_z_xz_xyz = pbuffer.data(idx_dip_df + 144);

    auto tr_z_xz_xzz = pbuffer.data(idx_dip_df + 145);

    auto tr_z_xz_yyy = pbuffer.data(idx_dip_df + 146);

    auto tr_z_xz_yyz = pbuffer.data(idx_dip_df + 147);

    auto tr_z_xz_yzz = pbuffer.data(idx_dip_df + 148);

    auto tr_z_xz_zzz = pbuffer.data(idx_dip_df + 149);

    auto tr_z_yy_xxx = pbuffer.data(idx_dip_df + 150);

    auto tr_z_yy_xxy = pbuffer.data(idx_dip_df + 151);

    auto tr_z_yy_xxz = pbuffer.data(idx_dip_df + 152);

    auto tr_z_yy_xyy = pbuffer.data(idx_dip_df + 153);

    auto tr_z_yy_xyz = pbuffer.data(idx_dip_df + 154);

    auto tr_z_yy_xzz = pbuffer.data(idx_dip_df + 155);

    auto tr_z_yy_yyy = pbuffer.data(idx_dip_df + 156);

    auto tr_z_yy_yyz = pbuffer.data(idx_dip_df + 157);

    auto tr_z_yy_yzz = pbuffer.data(idx_dip_df + 158);

    auto tr_z_yy_zzz = pbuffer.data(idx_dip_df + 159);

    auto tr_z_yz_xxx = pbuffer.data(idx_dip_df + 160);

    auto tr_z_yz_xxy = pbuffer.data(idx_dip_df + 161);

    auto tr_z_yz_xxz = pbuffer.data(idx_dip_df + 162);

    auto tr_z_yz_xyy = pbuffer.data(idx_dip_df + 163);

    auto tr_z_yz_xyz = pbuffer.data(idx_dip_df + 164);

    auto tr_z_yz_xzz = pbuffer.data(idx_dip_df + 165);

    auto tr_z_yz_yyy = pbuffer.data(idx_dip_df + 166);

    auto tr_z_yz_yyz = pbuffer.data(idx_dip_df + 167);

    auto tr_z_yz_yzz = pbuffer.data(idx_dip_df + 168);

    auto tr_z_yz_zzz = pbuffer.data(idx_dip_df + 169);

    auto tr_z_zz_xxx = pbuffer.data(idx_dip_df + 170);

    auto tr_z_zz_xxy = pbuffer.data(idx_dip_df + 171);

    auto tr_z_zz_xxz = pbuffer.data(idx_dip_df + 172);

    auto tr_z_zz_xyy = pbuffer.data(idx_dip_df + 173);

    auto tr_z_zz_xyz = pbuffer.data(idx_dip_df + 174);

    auto tr_z_zz_xzz = pbuffer.data(idx_dip_df + 175);

    auto tr_z_zz_yyy = pbuffer.data(idx_dip_df + 176);

    auto tr_z_zz_yyz = pbuffer.data(idx_dip_df + 177);

    auto tr_z_zz_yzz = pbuffer.data(idx_dip_df + 178);

    auto tr_z_zz_zzz = pbuffer.data(idx_dip_df + 179);

    // Set up 0-10 components of targeted buffer : FF

    auto tr_x_xxx_xxx = pbuffer.data(idx_dip_ff);

    auto tr_x_xxx_xxy = pbuffer.data(idx_dip_ff + 1);

    auto tr_x_xxx_xxz = pbuffer.data(idx_dip_ff + 2);

    auto tr_x_xxx_xyy = pbuffer.data(idx_dip_ff + 3);

    auto tr_x_xxx_xyz = pbuffer.data(idx_dip_ff + 4);

    auto tr_x_xxx_xzz = pbuffer.data(idx_dip_ff + 5);

    auto tr_x_xxx_yyy = pbuffer.data(idx_dip_ff + 6);

    auto tr_x_xxx_yyz = pbuffer.data(idx_dip_ff + 7);

    auto tr_x_xxx_yzz = pbuffer.data(idx_dip_ff + 8);

    auto tr_x_xxx_zzz = pbuffer.data(idx_dip_ff + 9);

    #pragma omp simd aligned(pa_x, tr_x_x_xxx, tr_x_x_xxy, tr_x_x_xxz, tr_x_x_xyy, tr_x_x_xyz, tr_x_x_xzz, tr_x_x_yyy, tr_x_x_yyz, tr_x_x_yzz, tr_x_x_zzz, tr_x_xx_xx, tr_x_xx_xxx, tr_x_xx_xxy, tr_x_xx_xxz, tr_x_xx_xy, tr_x_xx_xyy, tr_x_xx_xyz, tr_x_xx_xz, tr_x_xx_xzz, tr_x_xx_yy, tr_x_xx_yyy, tr_x_xx_yyz, tr_x_xx_yz, tr_x_xx_yzz, tr_x_xx_zz, tr_x_xx_zzz, tr_x_xxx_xxx, tr_x_xxx_xxy, tr_x_xxx_xxz, tr_x_xxx_xyy, tr_x_xxx_xyz, tr_x_xxx_xzz, tr_x_xxx_yyy, tr_x_xxx_yyz, tr_x_xxx_yzz, tr_x_xxx_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xxz, ts_xx_xyy, ts_xx_xyz, ts_xx_xzz, ts_xx_yyy, ts_xx_yyz, ts_xx_yzz, ts_xx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxx_xxx[i] = 2.0 * tr_x_x_xxx[i] * fe_0 + 3.0 * tr_x_xx_xx[i] * fe_0 + ts_xx_xxx[i] * fe_0 + tr_x_xx_xxx[i] * pa_x[i];

        tr_x_xxx_xxy[i] = 2.0 * tr_x_x_xxy[i] * fe_0 + 2.0 * tr_x_xx_xy[i] * fe_0 + ts_xx_xxy[i] * fe_0 + tr_x_xx_xxy[i] * pa_x[i];

        tr_x_xxx_xxz[i] = 2.0 * tr_x_x_xxz[i] * fe_0 + 2.0 * tr_x_xx_xz[i] * fe_0 + ts_xx_xxz[i] * fe_0 + tr_x_xx_xxz[i] * pa_x[i];

        tr_x_xxx_xyy[i] = 2.0 * tr_x_x_xyy[i] * fe_0 + tr_x_xx_yy[i] * fe_0 + ts_xx_xyy[i] * fe_0 + tr_x_xx_xyy[i] * pa_x[i];

        tr_x_xxx_xyz[i] = 2.0 * tr_x_x_xyz[i] * fe_0 + tr_x_xx_yz[i] * fe_0 + ts_xx_xyz[i] * fe_0 + tr_x_xx_xyz[i] * pa_x[i];

        tr_x_xxx_xzz[i] = 2.0 * tr_x_x_xzz[i] * fe_0 + tr_x_xx_zz[i] * fe_0 + ts_xx_xzz[i] * fe_0 + tr_x_xx_xzz[i] * pa_x[i];

        tr_x_xxx_yyy[i] = 2.0 * tr_x_x_yyy[i] * fe_0 + ts_xx_yyy[i] * fe_0 + tr_x_xx_yyy[i] * pa_x[i];

        tr_x_xxx_yyz[i] = 2.0 * tr_x_x_yyz[i] * fe_0 + ts_xx_yyz[i] * fe_0 + tr_x_xx_yyz[i] * pa_x[i];

        tr_x_xxx_yzz[i] = 2.0 * tr_x_x_yzz[i] * fe_0 + ts_xx_yzz[i] * fe_0 + tr_x_xx_yzz[i] * pa_x[i];

        tr_x_xxx_zzz[i] = 2.0 * tr_x_x_zzz[i] * fe_0 + ts_xx_zzz[i] * fe_0 + tr_x_xx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto tr_x_xxy_xxx = pbuffer.data(idx_dip_ff + 10);

    auto tr_x_xxy_xxy = pbuffer.data(idx_dip_ff + 11);

    auto tr_x_xxy_xxz = pbuffer.data(idx_dip_ff + 12);

    auto tr_x_xxy_xyy = pbuffer.data(idx_dip_ff + 13);

    auto tr_x_xxy_xyz = pbuffer.data(idx_dip_ff + 14);

    auto tr_x_xxy_xzz = pbuffer.data(idx_dip_ff + 15);

    auto tr_x_xxy_yyy = pbuffer.data(idx_dip_ff + 16);

    auto tr_x_xxy_yyz = pbuffer.data(idx_dip_ff + 17);

    auto tr_x_xxy_yzz = pbuffer.data(idx_dip_ff + 18);

    auto tr_x_xxy_zzz = pbuffer.data(idx_dip_ff + 19);

    #pragma omp simd aligned(pa_y, tr_x_xx_xx, tr_x_xx_xxx, tr_x_xx_xxy, tr_x_xx_xxz, tr_x_xx_xy, tr_x_xx_xyy, tr_x_xx_xyz, tr_x_xx_xz, tr_x_xx_xzz, tr_x_xx_yy, tr_x_xx_yyy, tr_x_xx_yyz, tr_x_xx_yz, tr_x_xx_yzz, tr_x_xx_zz, tr_x_xx_zzz, tr_x_xxy_xxx, tr_x_xxy_xxy, tr_x_xxy_xxz, tr_x_xxy_xyy, tr_x_xxy_xyz, tr_x_xxy_xzz, tr_x_xxy_yyy, tr_x_xxy_yyz, tr_x_xxy_yzz, tr_x_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxy_xxx[i] = tr_x_xx_xxx[i] * pa_y[i];

        tr_x_xxy_xxy[i] = tr_x_xx_xx[i] * fe_0 + tr_x_xx_xxy[i] * pa_y[i];

        tr_x_xxy_xxz[i] = tr_x_xx_xxz[i] * pa_y[i];

        tr_x_xxy_xyy[i] = 2.0 * tr_x_xx_xy[i] * fe_0 + tr_x_xx_xyy[i] * pa_y[i];

        tr_x_xxy_xyz[i] = tr_x_xx_xz[i] * fe_0 + tr_x_xx_xyz[i] * pa_y[i];

        tr_x_xxy_xzz[i] = tr_x_xx_xzz[i] * pa_y[i];

        tr_x_xxy_yyy[i] = 3.0 * tr_x_xx_yy[i] * fe_0 + tr_x_xx_yyy[i] * pa_y[i];

        tr_x_xxy_yyz[i] = 2.0 * tr_x_xx_yz[i] * fe_0 + tr_x_xx_yyz[i] * pa_y[i];

        tr_x_xxy_yzz[i] = tr_x_xx_zz[i] * fe_0 + tr_x_xx_yzz[i] * pa_y[i];

        tr_x_xxy_zzz[i] = tr_x_xx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto tr_x_xxz_xxx = pbuffer.data(idx_dip_ff + 20);

    auto tr_x_xxz_xxy = pbuffer.data(idx_dip_ff + 21);

    auto tr_x_xxz_xxz = pbuffer.data(idx_dip_ff + 22);

    auto tr_x_xxz_xyy = pbuffer.data(idx_dip_ff + 23);

    auto tr_x_xxz_xyz = pbuffer.data(idx_dip_ff + 24);

    auto tr_x_xxz_xzz = pbuffer.data(idx_dip_ff + 25);

    auto tr_x_xxz_yyy = pbuffer.data(idx_dip_ff + 26);

    auto tr_x_xxz_yyz = pbuffer.data(idx_dip_ff + 27);

    auto tr_x_xxz_yzz = pbuffer.data(idx_dip_ff + 28);

    auto tr_x_xxz_zzz = pbuffer.data(idx_dip_ff + 29);

    #pragma omp simd aligned(pa_z, tr_x_xx_xx, tr_x_xx_xxx, tr_x_xx_xxy, tr_x_xx_xxz, tr_x_xx_xy, tr_x_xx_xyy, tr_x_xx_xyz, tr_x_xx_xz, tr_x_xx_xzz, tr_x_xx_yy, tr_x_xx_yyy, tr_x_xx_yyz, tr_x_xx_yz, tr_x_xx_yzz, tr_x_xx_zz, tr_x_xx_zzz, tr_x_xxz_xxx, tr_x_xxz_xxy, tr_x_xxz_xxz, tr_x_xxz_xyy, tr_x_xxz_xyz, tr_x_xxz_xzz, tr_x_xxz_yyy, tr_x_xxz_yyz, tr_x_xxz_yzz, tr_x_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxz_xxx[i] = tr_x_xx_xxx[i] * pa_z[i];

        tr_x_xxz_xxy[i] = tr_x_xx_xxy[i] * pa_z[i];

        tr_x_xxz_xxz[i] = tr_x_xx_xx[i] * fe_0 + tr_x_xx_xxz[i] * pa_z[i];

        tr_x_xxz_xyy[i] = tr_x_xx_xyy[i] * pa_z[i];

        tr_x_xxz_xyz[i] = tr_x_xx_xy[i] * fe_0 + tr_x_xx_xyz[i] * pa_z[i];

        tr_x_xxz_xzz[i] = 2.0 * tr_x_xx_xz[i] * fe_0 + tr_x_xx_xzz[i] * pa_z[i];

        tr_x_xxz_yyy[i] = tr_x_xx_yyy[i] * pa_z[i];

        tr_x_xxz_yyz[i] = tr_x_xx_yy[i] * fe_0 + tr_x_xx_yyz[i] * pa_z[i];

        tr_x_xxz_yzz[i] = 2.0 * tr_x_xx_yz[i] * fe_0 + tr_x_xx_yzz[i] * pa_z[i];

        tr_x_xxz_zzz[i] = 3.0 * tr_x_xx_zz[i] * fe_0 + tr_x_xx_zzz[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto tr_x_xyy_xxx = pbuffer.data(idx_dip_ff + 30);

    auto tr_x_xyy_xxy = pbuffer.data(idx_dip_ff + 31);

    auto tr_x_xyy_xxz = pbuffer.data(idx_dip_ff + 32);

    auto tr_x_xyy_xyy = pbuffer.data(idx_dip_ff + 33);

    auto tr_x_xyy_xyz = pbuffer.data(idx_dip_ff + 34);

    auto tr_x_xyy_xzz = pbuffer.data(idx_dip_ff + 35);

    auto tr_x_xyy_yyy = pbuffer.data(idx_dip_ff + 36);

    auto tr_x_xyy_yyz = pbuffer.data(idx_dip_ff + 37);

    auto tr_x_xyy_yzz = pbuffer.data(idx_dip_ff + 38);

    auto tr_x_xyy_zzz = pbuffer.data(idx_dip_ff + 39);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_x_xxx, tr_x_x_xxz, tr_x_x_xzz, tr_x_xy_xxx, tr_x_xy_xxz, tr_x_xy_xzz, tr_x_xyy_xxx, tr_x_xyy_xxy, tr_x_xyy_xxz, tr_x_xyy_xyy, tr_x_xyy_xyz, tr_x_xyy_xzz, tr_x_xyy_yyy, tr_x_xyy_yyz, tr_x_xyy_yzz, tr_x_xyy_zzz, tr_x_yy_xxy, tr_x_yy_xy, tr_x_yy_xyy, tr_x_yy_xyz, tr_x_yy_yy, tr_x_yy_yyy, tr_x_yy_yyz, tr_x_yy_yz, tr_x_yy_yzz, tr_x_yy_zzz, ts_yy_xxy, ts_yy_xyy, ts_yy_xyz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyy_xxx[i] = tr_x_x_xxx[i] * fe_0 + tr_x_xy_xxx[i] * pa_y[i];

        tr_x_xyy_xxy[i] = 2.0 * tr_x_yy_xy[i] * fe_0 + ts_yy_xxy[i] * fe_0 + tr_x_yy_xxy[i] * pa_x[i];

        tr_x_xyy_xxz[i] = tr_x_x_xxz[i] * fe_0 + tr_x_xy_xxz[i] * pa_y[i];

        tr_x_xyy_xyy[i] = tr_x_yy_yy[i] * fe_0 + ts_yy_xyy[i] * fe_0 + tr_x_yy_xyy[i] * pa_x[i];

        tr_x_xyy_xyz[i] = tr_x_yy_yz[i] * fe_0 + ts_yy_xyz[i] * fe_0 + tr_x_yy_xyz[i] * pa_x[i];

        tr_x_xyy_xzz[i] = tr_x_x_xzz[i] * fe_0 + tr_x_xy_xzz[i] * pa_y[i];

        tr_x_xyy_yyy[i] = ts_yy_yyy[i] * fe_0 + tr_x_yy_yyy[i] * pa_x[i];

        tr_x_xyy_yyz[i] = ts_yy_yyz[i] * fe_0 + tr_x_yy_yyz[i] * pa_x[i];

        tr_x_xyy_yzz[i] = ts_yy_yzz[i] * fe_0 + tr_x_yy_yzz[i] * pa_x[i];

        tr_x_xyy_zzz[i] = ts_yy_zzz[i] * fe_0 + tr_x_yy_zzz[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto tr_x_xyz_xxx = pbuffer.data(idx_dip_ff + 40);

    auto tr_x_xyz_xxy = pbuffer.data(idx_dip_ff + 41);

    auto tr_x_xyz_xxz = pbuffer.data(idx_dip_ff + 42);

    auto tr_x_xyz_xyy = pbuffer.data(idx_dip_ff + 43);

    auto tr_x_xyz_xyz = pbuffer.data(idx_dip_ff + 44);

    auto tr_x_xyz_xzz = pbuffer.data(idx_dip_ff + 45);

    auto tr_x_xyz_yyy = pbuffer.data(idx_dip_ff + 46);

    auto tr_x_xyz_yyz = pbuffer.data(idx_dip_ff + 47);

    auto tr_x_xyz_yzz = pbuffer.data(idx_dip_ff + 48);

    auto tr_x_xyz_zzz = pbuffer.data(idx_dip_ff + 49);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xy_xxy, tr_x_xy_xyy, tr_x_xy_yyy, tr_x_xyz_xxx, tr_x_xyz_xxy, tr_x_xyz_xxz, tr_x_xyz_xyy, tr_x_xyz_xyz, tr_x_xyz_xzz, tr_x_xyz_yyy, tr_x_xyz_yyz, tr_x_xyz_yzz, tr_x_xyz_zzz, tr_x_xz_xxx, tr_x_xz_xxz, tr_x_xz_xyz, tr_x_xz_xz, tr_x_xz_xzz, tr_x_xz_zzz, tr_x_yz_yyz, tr_x_yz_yzz, ts_yz_yyz, ts_yz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyz_xxx[i] = tr_x_xz_xxx[i] * pa_y[i];

        tr_x_xyz_xxy[i] = tr_x_xy_xxy[i] * pa_z[i];

        tr_x_xyz_xxz[i] = tr_x_xz_xxz[i] * pa_y[i];

        tr_x_xyz_xyy[i] = tr_x_xy_xyy[i] * pa_z[i];

        tr_x_xyz_xyz[i] = tr_x_xz_xz[i] * fe_0 + tr_x_xz_xyz[i] * pa_y[i];

        tr_x_xyz_xzz[i] = tr_x_xz_xzz[i] * pa_y[i];

        tr_x_xyz_yyy[i] = tr_x_xy_yyy[i] * pa_z[i];

        tr_x_xyz_yyz[i] = ts_yz_yyz[i] * fe_0 + tr_x_yz_yyz[i] * pa_x[i];

        tr_x_xyz_yzz[i] = ts_yz_yzz[i] * fe_0 + tr_x_yz_yzz[i] * pa_x[i];

        tr_x_xyz_zzz[i] = tr_x_xz_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto tr_x_xzz_xxx = pbuffer.data(idx_dip_ff + 50);

    auto tr_x_xzz_xxy = pbuffer.data(idx_dip_ff + 51);

    auto tr_x_xzz_xxz = pbuffer.data(idx_dip_ff + 52);

    auto tr_x_xzz_xyy = pbuffer.data(idx_dip_ff + 53);

    auto tr_x_xzz_xyz = pbuffer.data(idx_dip_ff + 54);

    auto tr_x_xzz_xzz = pbuffer.data(idx_dip_ff + 55);

    auto tr_x_xzz_yyy = pbuffer.data(idx_dip_ff + 56);

    auto tr_x_xzz_yyz = pbuffer.data(idx_dip_ff + 57);

    auto tr_x_xzz_yzz = pbuffer.data(idx_dip_ff + 58);

    auto tr_x_xzz_zzz = pbuffer.data(idx_dip_ff + 59);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_x_xxx, tr_x_x_xxy, tr_x_x_xyy, tr_x_xz_xxx, tr_x_xz_xxy, tr_x_xz_xyy, tr_x_xzz_xxx, tr_x_xzz_xxy, tr_x_xzz_xxz, tr_x_xzz_xyy, tr_x_xzz_xyz, tr_x_xzz_xzz, tr_x_xzz_yyy, tr_x_xzz_yyz, tr_x_xzz_yzz, tr_x_xzz_zzz, tr_x_zz_xxz, tr_x_zz_xyz, tr_x_zz_xz, tr_x_zz_xzz, tr_x_zz_yyy, tr_x_zz_yyz, tr_x_zz_yz, tr_x_zz_yzz, tr_x_zz_zz, tr_x_zz_zzz, ts_zz_xxz, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzz_xxx[i] = tr_x_x_xxx[i] * fe_0 + tr_x_xz_xxx[i] * pa_z[i];

        tr_x_xzz_xxy[i] = tr_x_x_xxy[i] * fe_0 + tr_x_xz_xxy[i] * pa_z[i];

        tr_x_xzz_xxz[i] = 2.0 * tr_x_zz_xz[i] * fe_0 + ts_zz_xxz[i] * fe_0 + tr_x_zz_xxz[i] * pa_x[i];

        tr_x_xzz_xyy[i] = tr_x_x_xyy[i] * fe_0 + tr_x_xz_xyy[i] * pa_z[i];

        tr_x_xzz_xyz[i] = tr_x_zz_yz[i] * fe_0 + ts_zz_xyz[i] * fe_0 + tr_x_zz_xyz[i] * pa_x[i];

        tr_x_xzz_xzz[i] = tr_x_zz_zz[i] * fe_0 + ts_zz_xzz[i] * fe_0 + tr_x_zz_xzz[i] * pa_x[i];

        tr_x_xzz_yyy[i] = ts_zz_yyy[i] * fe_0 + tr_x_zz_yyy[i] * pa_x[i];

        tr_x_xzz_yyz[i] = ts_zz_yyz[i] * fe_0 + tr_x_zz_yyz[i] * pa_x[i];

        tr_x_xzz_yzz[i] = ts_zz_yzz[i] * fe_0 + tr_x_zz_yzz[i] * pa_x[i];

        tr_x_xzz_zzz[i] = ts_zz_zzz[i] * fe_0 + tr_x_zz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto tr_x_yyy_xxx = pbuffer.data(idx_dip_ff + 60);

    auto tr_x_yyy_xxy = pbuffer.data(idx_dip_ff + 61);

    auto tr_x_yyy_xxz = pbuffer.data(idx_dip_ff + 62);

    auto tr_x_yyy_xyy = pbuffer.data(idx_dip_ff + 63);

    auto tr_x_yyy_xyz = pbuffer.data(idx_dip_ff + 64);

    auto tr_x_yyy_xzz = pbuffer.data(idx_dip_ff + 65);

    auto tr_x_yyy_yyy = pbuffer.data(idx_dip_ff + 66);

    auto tr_x_yyy_yyz = pbuffer.data(idx_dip_ff + 67);

    auto tr_x_yyy_yzz = pbuffer.data(idx_dip_ff + 68);

    auto tr_x_yyy_zzz = pbuffer.data(idx_dip_ff + 69);

    #pragma omp simd aligned(pa_y, tr_x_y_xxx, tr_x_y_xxy, tr_x_y_xxz, tr_x_y_xyy, tr_x_y_xyz, tr_x_y_xzz, tr_x_y_yyy, tr_x_y_yyz, tr_x_y_yzz, tr_x_y_zzz, tr_x_yy_xx, tr_x_yy_xxx, tr_x_yy_xxy, tr_x_yy_xxz, tr_x_yy_xy, tr_x_yy_xyy, tr_x_yy_xyz, tr_x_yy_xz, tr_x_yy_xzz, tr_x_yy_yy, tr_x_yy_yyy, tr_x_yy_yyz, tr_x_yy_yz, tr_x_yy_yzz, tr_x_yy_zz, tr_x_yy_zzz, tr_x_yyy_xxx, tr_x_yyy_xxy, tr_x_yyy_xxz, tr_x_yyy_xyy, tr_x_yyy_xyz, tr_x_yyy_xzz, tr_x_yyy_yyy, tr_x_yyy_yyz, tr_x_yyy_yzz, tr_x_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyy_xxx[i] = 2.0 * tr_x_y_xxx[i] * fe_0 + tr_x_yy_xxx[i] * pa_y[i];

        tr_x_yyy_xxy[i] = 2.0 * tr_x_y_xxy[i] * fe_0 + tr_x_yy_xx[i] * fe_0 + tr_x_yy_xxy[i] * pa_y[i];

        tr_x_yyy_xxz[i] = 2.0 * tr_x_y_xxz[i] * fe_0 + tr_x_yy_xxz[i] * pa_y[i];

        tr_x_yyy_xyy[i] = 2.0 * tr_x_y_xyy[i] * fe_0 + 2.0 * tr_x_yy_xy[i] * fe_0 + tr_x_yy_xyy[i] * pa_y[i];

        tr_x_yyy_xyz[i] = 2.0 * tr_x_y_xyz[i] * fe_0 + tr_x_yy_xz[i] * fe_0 + tr_x_yy_xyz[i] * pa_y[i];

        tr_x_yyy_xzz[i] = 2.0 * tr_x_y_xzz[i] * fe_0 + tr_x_yy_xzz[i] * pa_y[i];

        tr_x_yyy_yyy[i] = 2.0 * tr_x_y_yyy[i] * fe_0 + 3.0 * tr_x_yy_yy[i] * fe_0 + tr_x_yy_yyy[i] * pa_y[i];

        tr_x_yyy_yyz[i] = 2.0 * tr_x_y_yyz[i] * fe_0 + 2.0 * tr_x_yy_yz[i] * fe_0 + tr_x_yy_yyz[i] * pa_y[i];

        tr_x_yyy_yzz[i] = 2.0 * tr_x_y_yzz[i] * fe_0 + tr_x_yy_zz[i] * fe_0 + tr_x_yy_yzz[i] * pa_y[i];

        tr_x_yyy_zzz[i] = 2.0 * tr_x_y_zzz[i] * fe_0 + tr_x_yy_zzz[i] * pa_y[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto tr_x_yyz_xxx = pbuffer.data(idx_dip_ff + 70);

    auto tr_x_yyz_xxy = pbuffer.data(idx_dip_ff + 71);

    auto tr_x_yyz_xxz = pbuffer.data(idx_dip_ff + 72);

    auto tr_x_yyz_xyy = pbuffer.data(idx_dip_ff + 73);

    auto tr_x_yyz_xyz = pbuffer.data(idx_dip_ff + 74);

    auto tr_x_yyz_xzz = pbuffer.data(idx_dip_ff + 75);

    auto tr_x_yyz_yyy = pbuffer.data(idx_dip_ff + 76);

    auto tr_x_yyz_yyz = pbuffer.data(idx_dip_ff + 77);

    auto tr_x_yyz_yzz = pbuffer.data(idx_dip_ff + 78);

    auto tr_x_yyz_zzz = pbuffer.data(idx_dip_ff + 79);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yy_xxx, tr_x_yy_xxy, tr_x_yy_xy, tr_x_yy_xyy, tr_x_yy_xyz, tr_x_yy_yy, tr_x_yy_yyy, tr_x_yy_yyz, tr_x_yy_yz, tr_x_yy_yzz, tr_x_yyz_xxx, tr_x_yyz_xxy, tr_x_yyz_xxz, tr_x_yyz_xyy, tr_x_yyz_xyz, tr_x_yyz_xzz, tr_x_yyz_yyy, tr_x_yyz_yyz, tr_x_yyz_yzz, tr_x_yyz_zzz, tr_x_yz_xxz, tr_x_yz_xzz, tr_x_yz_zzz, tr_x_z_xxz, tr_x_z_xzz, tr_x_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyz_xxx[i] = tr_x_yy_xxx[i] * pa_z[i];

        tr_x_yyz_xxy[i] = tr_x_yy_xxy[i] * pa_z[i];

        tr_x_yyz_xxz[i] = tr_x_z_xxz[i] * fe_0 + tr_x_yz_xxz[i] * pa_y[i];

        tr_x_yyz_xyy[i] = tr_x_yy_xyy[i] * pa_z[i];

        tr_x_yyz_xyz[i] = tr_x_yy_xy[i] * fe_0 + tr_x_yy_xyz[i] * pa_z[i];

        tr_x_yyz_xzz[i] = tr_x_z_xzz[i] * fe_0 + tr_x_yz_xzz[i] * pa_y[i];

        tr_x_yyz_yyy[i] = tr_x_yy_yyy[i] * pa_z[i];

        tr_x_yyz_yyz[i] = tr_x_yy_yy[i] * fe_0 + tr_x_yy_yyz[i] * pa_z[i];

        tr_x_yyz_yzz[i] = 2.0 * tr_x_yy_yz[i] * fe_0 + tr_x_yy_yzz[i] * pa_z[i];

        tr_x_yyz_zzz[i] = tr_x_z_zzz[i] * fe_0 + tr_x_yz_zzz[i] * pa_y[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto tr_x_yzz_xxx = pbuffer.data(idx_dip_ff + 80);

    auto tr_x_yzz_xxy = pbuffer.data(idx_dip_ff + 81);

    auto tr_x_yzz_xxz = pbuffer.data(idx_dip_ff + 82);

    auto tr_x_yzz_xyy = pbuffer.data(idx_dip_ff + 83);

    auto tr_x_yzz_xyz = pbuffer.data(idx_dip_ff + 84);

    auto tr_x_yzz_xzz = pbuffer.data(idx_dip_ff + 85);

    auto tr_x_yzz_yyy = pbuffer.data(idx_dip_ff + 86);

    auto tr_x_yzz_yyz = pbuffer.data(idx_dip_ff + 87);

    auto tr_x_yzz_yzz = pbuffer.data(idx_dip_ff + 88);

    auto tr_x_yzz_zzz = pbuffer.data(idx_dip_ff + 89);

    #pragma omp simd aligned(pa_y, tr_x_yzz_xxx, tr_x_yzz_xxy, tr_x_yzz_xxz, tr_x_yzz_xyy, tr_x_yzz_xyz, tr_x_yzz_xzz, tr_x_yzz_yyy, tr_x_yzz_yyz, tr_x_yzz_yzz, tr_x_yzz_zzz, tr_x_zz_xx, tr_x_zz_xxx, tr_x_zz_xxy, tr_x_zz_xxz, tr_x_zz_xy, tr_x_zz_xyy, tr_x_zz_xyz, tr_x_zz_xz, tr_x_zz_xzz, tr_x_zz_yy, tr_x_zz_yyy, tr_x_zz_yyz, tr_x_zz_yz, tr_x_zz_yzz, tr_x_zz_zz, tr_x_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzz_xxx[i] = tr_x_zz_xxx[i] * pa_y[i];

        tr_x_yzz_xxy[i] = tr_x_zz_xx[i] * fe_0 + tr_x_zz_xxy[i] * pa_y[i];

        tr_x_yzz_xxz[i] = tr_x_zz_xxz[i] * pa_y[i];

        tr_x_yzz_xyy[i] = 2.0 * tr_x_zz_xy[i] * fe_0 + tr_x_zz_xyy[i] * pa_y[i];

        tr_x_yzz_xyz[i] = tr_x_zz_xz[i] * fe_0 + tr_x_zz_xyz[i] * pa_y[i];

        tr_x_yzz_xzz[i] = tr_x_zz_xzz[i] * pa_y[i];

        tr_x_yzz_yyy[i] = 3.0 * tr_x_zz_yy[i] * fe_0 + tr_x_zz_yyy[i] * pa_y[i];

        tr_x_yzz_yyz[i] = 2.0 * tr_x_zz_yz[i] * fe_0 + tr_x_zz_yyz[i] * pa_y[i];

        tr_x_yzz_yzz[i] = tr_x_zz_zz[i] * fe_0 + tr_x_zz_yzz[i] * pa_y[i];

        tr_x_yzz_zzz[i] = tr_x_zz_zzz[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto tr_x_zzz_xxx = pbuffer.data(idx_dip_ff + 90);

    auto tr_x_zzz_xxy = pbuffer.data(idx_dip_ff + 91);

    auto tr_x_zzz_xxz = pbuffer.data(idx_dip_ff + 92);

    auto tr_x_zzz_xyy = pbuffer.data(idx_dip_ff + 93);

    auto tr_x_zzz_xyz = pbuffer.data(idx_dip_ff + 94);

    auto tr_x_zzz_xzz = pbuffer.data(idx_dip_ff + 95);

    auto tr_x_zzz_yyy = pbuffer.data(idx_dip_ff + 96);

    auto tr_x_zzz_yyz = pbuffer.data(idx_dip_ff + 97);

    auto tr_x_zzz_yzz = pbuffer.data(idx_dip_ff + 98);

    auto tr_x_zzz_zzz = pbuffer.data(idx_dip_ff + 99);

    #pragma omp simd aligned(pa_z, tr_x_z_xxx, tr_x_z_xxy, tr_x_z_xxz, tr_x_z_xyy, tr_x_z_xyz, tr_x_z_xzz, tr_x_z_yyy, tr_x_z_yyz, tr_x_z_yzz, tr_x_z_zzz, tr_x_zz_xx, tr_x_zz_xxx, tr_x_zz_xxy, tr_x_zz_xxz, tr_x_zz_xy, tr_x_zz_xyy, tr_x_zz_xyz, tr_x_zz_xz, tr_x_zz_xzz, tr_x_zz_yy, tr_x_zz_yyy, tr_x_zz_yyz, tr_x_zz_yz, tr_x_zz_yzz, tr_x_zz_zz, tr_x_zz_zzz, tr_x_zzz_xxx, tr_x_zzz_xxy, tr_x_zzz_xxz, tr_x_zzz_xyy, tr_x_zzz_xyz, tr_x_zzz_xzz, tr_x_zzz_yyy, tr_x_zzz_yyz, tr_x_zzz_yzz, tr_x_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzz_xxx[i] = 2.0 * tr_x_z_xxx[i] * fe_0 + tr_x_zz_xxx[i] * pa_z[i];

        tr_x_zzz_xxy[i] = 2.0 * tr_x_z_xxy[i] * fe_0 + tr_x_zz_xxy[i] * pa_z[i];

        tr_x_zzz_xxz[i] = 2.0 * tr_x_z_xxz[i] * fe_0 + tr_x_zz_xx[i] * fe_0 + tr_x_zz_xxz[i] * pa_z[i];

        tr_x_zzz_xyy[i] = 2.0 * tr_x_z_xyy[i] * fe_0 + tr_x_zz_xyy[i] * pa_z[i];

        tr_x_zzz_xyz[i] = 2.0 * tr_x_z_xyz[i] * fe_0 + tr_x_zz_xy[i] * fe_0 + tr_x_zz_xyz[i] * pa_z[i];

        tr_x_zzz_xzz[i] = 2.0 * tr_x_z_xzz[i] * fe_0 + 2.0 * tr_x_zz_xz[i] * fe_0 + tr_x_zz_xzz[i] * pa_z[i];

        tr_x_zzz_yyy[i] = 2.0 * tr_x_z_yyy[i] * fe_0 + tr_x_zz_yyy[i] * pa_z[i];

        tr_x_zzz_yyz[i] = 2.0 * tr_x_z_yyz[i] * fe_0 + tr_x_zz_yy[i] * fe_0 + tr_x_zz_yyz[i] * pa_z[i];

        tr_x_zzz_yzz[i] = 2.0 * tr_x_z_yzz[i] * fe_0 + 2.0 * tr_x_zz_yz[i] * fe_0 + tr_x_zz_yzz[i] * pa_z[i];

        tr_x_zzz_zzz[i] = 2.0 * tr_x_z_zzz[i] * fe_0 + 3.0 * tr_x_zz_zz[i] * fe_0 + tr_x_zz_zzz[i] * pa_z[i];
    }

    // Set up 100-110 components of targeted buffer : FF

    auto tr_y_xxx_xxx = pbuffer.data(idx_dip_ff + 100);

    auto tr_y_xxx_xxy = pbuffer.data(idx_dip_ff + 101);

    auto tr_y_xxx_xxz = pbuffer.data(idx_dip_ff + 102);

    auto tr_y_xxx_xyy = pbuffer.data(idx_dip_ff + 103);

    auto tr_y_xxx_xyz = pbuffer.data(idx_dip_ff + 104);

    auto tr_y_xxx_xzz = pbuffer.data(idx_dip_ff + 105);

    auto tr_y_xxx_yyy = pbuffer.data(idx_dip_ff + 106);

    auto tr_y_xxx_yyz = pbuffer.data(idx_dip_ff + 107);

    auto tr_y_xxx_yzz = pbuffer.data(idx_dip_ff + 108);

    auto tr_y_xxx_zzz = pbuffer.data(idx_dip_ff + 109);

    #pragma omp simd aligned(pa_x, tr_y_x_xxx, tr_y_x_xxy, tr_y_x_xxz, tr_y_x_xyy, tr_y_x_xyz, tr_y_x_xzz, tr_y_x_yyy, tr_y_x_yyz, tr_y_x_yzz, tr_y_x_zzz, tr_y_xx_xx, tr_y_xx_xxx, tr_y_xx_xxy, tr_y_xx_xxz, tr_y_xx_xy, tr_y_xx_xyy, tr_y_xx_xyz, tr_y_xx_xz, tr_y_xx_xzz, tr_y_xx_yy, tr_y_xx_yyy, tr_y_xx_yyz, tr_y_xx_yz, tr_y_xx_yzz, tr_y_xx_zz, tr_y_xx_zzz, tr_y_xxx_xxx, tr_y_xxx_xxy, tr_y_xxx_xxz, tr_y_xxx_xyy, tr_y_xxx_xyz, tr_y_xxx_xzz, tr_y_xxx_yyy, tr_y_xxx_yyz, tr_y_xxx_yzz, tr_y_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxx_xxx[i] = 2.0 * tr_y_x_xxx[i] * fe_0 + 3.0 * tr_y_xx_xx[i] * fe_0 + tr_y_xx_xxx[i] * pa_x[i];

        tr_y_xxx_xxy[i] = 2.0 * tr_y_x_xxy[i] * fe_0 + 2.0 * tr_y_xx_xy[i] * fe_0 + tr_y_xx_xxy[i] * pa_x[i];

        tr_y_xxx_xxz[i] = 2.0 * tr_y_x_xxz[i] * fe_0 + 2.0 * tr_y_xx_xz[i] * fe_0 + tr_y_xx_xxz[i] * pa_x[i];

        tr_y_xxx_xyy[i] = 2.0 * tr_y_x_xyy[i] * fe_0 + tr_y_xx_yy[i] * fe_0 + tr_y_xx_xyy[i] * pa_x[i];

        tr_y_xxx_xyz[i] = 2.0 * tr_y_x_xyz[i] * fe_0 + tr_y_xx_yz[i] * fe_0 + tr_y_xx_xyz[i] * pa_x[i];

        tr_y_xxx_xzz[i] = 2.0 * tr_y_x_xzz[i] * fe_0 + tr_y_xx_zz[i] * fe_0 + tr_y_xx_xzz[i] * pa_x[i];

        tr_y_xxx_yyy[i] = 2.0 * tr_y_x_yyy[i] * fe_0 + tr_y_xx_yyy[i] * pa_x[i];

        tr_y_xxx_yyz[i] = 2.0 * tr_y_x_yyz[i] * fe_0 + tr_y_xx_yyz[i] * pa_x[i];

        tr_y_xxx_yzz[i] = 2.0 * tr_y_x_yzz[i] * fe_0 + tr_y_xx_yzz[i] * pa_x[i];

        tr_y_xxx_zzz[i] = 2.0 * tr_y_x_zzz[i] * fe_0 + tr_y_xx_zzz[i] * pa_x[i];
    }

    // Set up 110-120 components of targeted buffer : FF

    auto tr_y_xxy_xxx = pbuffer.data(idx_dip_ff + 110);

    auto tr_y_xxy_xxy = pbuffer.data(idx_dip_ff + 111);

    auto tr_y_xxy_xxz = pbuffer.data(idx_dip_ff + 112);

    auto tr_y_xxy_xyy = pbuffer.data(idx_dip_ff + 113);

    auto tr_y_xxy_xyz = pbuffer.data(idx_dip_ff + 114);

    auto tr_y_xxy_xzz = pbuffer.data(idx_dip_ff + 115);

    auto tr_y_xxy_yyy = pbuffer.data(idx_dip_ff + 116);

    auto tr_y_xxy_yyz = pbuffer.data(idx_dip_ff + 117);

    auto tr_y_xxy_yzz = pbuffer.data(idx_dip_ff + 118);

    auto tr_y_xxy_zzz = pbuffer.data(idx_dip_ff + 119);

    #pragma omp simd aligned(pa_x, pa_y, tr_y_xx_xxx, tr_y_xx_xxz, tr_y_xx_xzz, tr_y_xxy_xxx, tr_y_xxy_xxy, tr_y_xxy_xxz, tr_y_xxy_xyy, tr_y_xxy_xyz, tr_y_xxy_xzz, tr_y_xxy_yyy, tr_y_xxy_yyz, tr_y_xxy_yzz, tr_y_xxy_zzz, tr_y_xy_xxy, tr_y_xy_xy, tr_y_xy_xyy, tr_y_xy_xyz, tr_y_xy_yy, tr_y_xy_yyy, tr_y_xy_yyz, tr_y_xy_yz, tr_y_xy_yzz, tr_y_xy_zzz, tr_y_y_xxy, tr_y_y_xyy, tr_y_y_xyz, tr_y_y_yyy, tr_y_y_yyz, tr_y_y_yzz, tr_y_y_zzz, ts_xx_xxx, ts_xx_xxz, ts_xx_xzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxy_xxx[i] = ts_xx_xxx[i] * fe_0 + tr_y_xx_xxx[i] * pa_y[i];

        tr_y_xxy_xxy[i] = tr_y_y_xxy[i] * fe_0 + 2.0 * tr_y_xy_xy[i] * fe_0 + tr_y_xy_xxy[i] * pa_x[i];

        tr_y_xxy_xxz[i] = ts_xx_xxz[i] * fe_0 + tr_y_xx_xxz[i] * pa_y[i];

        tr_y_xxy_xyy[i] = tr_y_y_xyy[i] * fe_0 + tr_y_xy_yy[i] * fe_0 + tr_y_xy_xyy[i] * pa_x[i];

        tr_y_xxy_xyz[i] = tr_y_y_xyz[i] * fe_0 + tr_y_xy_yz[i] * fe_0 + tr_y_xy_xyz[i] * pa_x[i];

        tr_y_xxy_xzz[i] = ts_xx_xzz[i] * fe_0 + tr_y_xx_xzz[i] * pa_y[i];

        tr_y_xxy_yyy[i] = tr_y_y_yyy[i] * fe_0 + tr_y_xy_yyy[i] * pa_x[i];

        tr_y_xxy_yyz[i] = tr_y_y_yyz[i] * fe_0 + tr_y_xy_yyz[i] * pa_x[i];

        tr_y_xxy_yzz[i] = tr_y_y_yzz[i] * fe_0 + tr_y_xy_yzz[i] * pa_x[i];

        tr_y_xxy_zzz[i] = tr_y_y_zzz[i] * fe_0 + tr_y_xy_zzz[i] * pa_x[i];
    }

    // Set up 120-130 components of targeted buffer : FF

    auto tr_y_xxz_xxx = pbuffer.data(idx_dip_ff + 120);

    auto tr_y_xxz_xxy = pbuffer.data(idx_dip_ff + 121);

    auto tr_y_xxz_xxz = pbuffer.data(idx_dip_ff + 122);

    auto tr_y_xxz_xyy = pbuffer.data(idx_dip_ff + 123);

    auto tr_y_xxz_xyz = pbuffer.data(idx_dip_ff + 124);

    auto tr_y_xxz_xzz = pbuffer.data(idx_dip_ff + 125);

    auto tr_y_xxz_yyy = pbuffer.data(idx_dip_ff + 126);

    auto tr_y_xxz_yyz = pbuffer.data(idx_dip_ff + 127);

    auto tr_y_xxz_yzz = pbuffer.data(idx_dip_ff + 128);

    auto tr_y_xxz_zzz = pbuffer.data(idx_dip_ff + 129);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xx_xx, tr_y_xx_xxx, tr_y_xx_xxy, tr_y_xx_xxz, tr_y_xx_xy, tr_y_xx_xyy, tr_y_xx_xyz, tr_y_xx_xz, tr_y_xx_xzz, tr_y_xx_yyy, tr_y_xxz_xxx, tr_y_xxz_xxy, tr_y_xxz_xxz, tr_y_xxz_xyy, tr_y_xxz_xyz, tr_y_xxz_xzz, tr_y_xxz_yyy, tr_y_xxz_yyz, tr_y_xxz_yzz, tr_y_xxz_zzz, tr_y_xz_yyz, tr_y_xz_yzz, tr_y_xz_zzz, tr_y_z_yyz, tr_y_z_yzz, tr_y_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxz_xxx[i] = tr_y_xx_xxx[i] * pa_z[i];

        tr_y_xxz_xxy[i] = tr_y_xx_xxy[i] * pa_z[i];

        tr_y_xxz_xxz[i] = tr_y_xx_xx[i] * fe_0 + tr_y_xx_xxz[i] * pa_z[i];

        tr_y_xxz_xyy[i] = tr_y_xx_xyy[i] * pa_z[i];

        tr_y_xxz_xyz[i] = tr_y_xx_xy[i] * fe_0 + tr_y_xx_xyz[i] * pa_z[i];

        tr_y_xxz_xzz[i] = 2.0 * tr_y_xx_xz[i] * fe_0 + tr_y_xx_xzz[i] * pa_z[i];

        tr_y_xxz_yyy[i] = tr_y_xx_yyy[i] * pa_z[i];

        tr_y_xxz_yyz[i] = tr_y_z_yyz[i] * fe_0 + tr_y_xz_yyz[i] * pa_x[i];

        tr_y_xxz_yzz[i] = tr_y_z_yzz[i] * fe_0 + tr_y_xz_yzz[i] * pa_x[i];

        tr_y_xxz_zzz[i] = tr_y_z_zzz[i] * fe_0 + tr_y_xz_zzz[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : FF

    auto tr_y_xyy_xxx = pbuffer.data(idx_dip_ff + 130);

    auto tr_y_xyy_xxy = pbuffer.data(idx_dip_ff + 131);

    auto tr_y_xyy_xxz = pbuffer.data(idx_dip_ff + 132);

    auto tr_y_xyy_xyy = pbuffer.data(idx_dip_ff + 133);

    auto tr_y_xyy_xyz = pbuffer.data(idx_dip_ff + 134);

    auto tr_y_xyy_xzz = pbuffer.data(idx_dip_ff + 135);

    auto tr_y_xyy_yyy = pbuffer.data(idx_dip_ff + 136);

    auto tr_y_xyy_yyz = pbuffer.data(idx_dip_ff + 137);

    auto tr_y_xyy_yzz = pbuffer.data(idx_dip_ff + 138);

    auto tr_y_xyy_zzz = pbuffer.data(idx_dip_ff + 139);

    #pragma omp simd aligned(pa_x, tr_y_xyy_xxx, tr_y_xyy_xxy, tr_y_xyy_xxz, tr_y_xyy_xyy, tr_y_xyy_xyz, tr_y_xyy_xzz, tr_y_xyy_yyy, tr_y_xyy_yyz, tr_y_xyy_yzz, tr_y_xyy_zzz, tr_y_yy_xx, tr_y_yy_xxx, tr_y_yy_xxy, tr_y_yy_xxz, tr_y_yy_xy, tr_y_yy_xyy, tr_y_yy_xyz, tr_y_yy_xz, tr_y_yy_xzz, tr_y_yy_yy, tr_y_yy_yyy, tr_y_yy_yyz, tr_y_yy_yz, tr_y_yy_yzz, tr_y_yy_zz, tr_y_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyy_xxx[i] = 3.0 * tr_y_yy_xx[i] * fe_0 + tr_y_yy_xxx[i] * pa_x[i];

        tr_y_xyy_xxy[i] = 2.0 * tr_y_yy_xy[i] * fe_0 + tr_y_yy_xxy[i] * pa_x[i];

        tr_y_xyy_xxz[i] = 2.0 * tr_y_yy_xz[i] * fe_0 + tr_y_yy_xxz[i] * pa_x[i];

        tr_y_xyy_xyy[i] = tr_y_yy_yy[i] * fe_0 + tr_y_yy_xyy[i] * pa_x[i];

        tr_y_xyy_xyz[i] = tr_y_yy_yz[i] * fe_0 + tr_y_yy_xyz[i] * pa_x[i];

        tr_y_xyy_xzz[i] = tr_y_yy_zz[i] * fe_0 + tr_y_yy_xzz[i] * pa_x[i];

        tr_y_xyy_yyy[i] = tr_y_yy_yyy[i] * pa_x[i];

        tr_y_xyy_yyz[i] = tr_y_yy_yyz[i] * pa_x[i];

        tr_y_xyy_yzz[i] = tr_y_yy_yzz[i] * pa_x[i];

        tr_y_xyy_zzz[i] = tr_y_yy_zzz[i] * pa_x[i];
    }

    // Set up 140-150 components of targeted buffer : FF

    auto tr_y_xyz_xxx = pbuffer.data(idx_dip_ff + 140);

    auto tr_y_xyz_xxy = pbuffer.data(idx_dip_ff + 141);

    auto tr_y_xyz_xxz = pbuffer.data(idx_dip_ff + 142);

    auto tr_y_xyz_xyy = pbuffer.data(idx_dip_ff + 143);

    auto tr_y_xyz_xyz = pbuffer.data(idx_dip_ff + 144);

    auto tr_y_xyz_xzz = pbuffer.data(idx_dip_ff + 145);

    auto tr_y_xyz_yyy = pbuffer.data(idx_dip_ff + 146);

    auto tr_y_xyz_yyz = pbuffer.data(idx_dip_ff + 147);

    auto tr_y_xyz_yzz = pbuffer.data(idx_dip_ff + 148);

    auto tr_y_xyz_zzz = pbuffer.data(idx_dip_ff + 149);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xy_xxx, tr_y_xy_xxy, tr_y_xy_xyy, tr_y_xyz_xxx, tr_y_xyz_xxy, tr_y_xyz_xxz, tr_y_xyz_xyy, tr_y_xyz_xyz, tr_y_xyz_xzz, tr_y_xyz_yyy, tr_y_xyz_yyz, tr_y_xyz_yzz, tr_y_xyz_zzz, tr_y_yz_xxz, tr_y_yz_xyz, tr_y_yz_xz, tr_y_yz_xzz, tr_y_yz_yyy, tr_y_yz_yyz, tr_y_yz_yz, tr_y_yz_yzz, tr_y_yz_zz, tr_y_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyz_xxx[i] = tr_y_xy_xxx[i] * pa_z[i];

        tr_y_xyz_xxy[i] = tr_y_xy_xxy[i] * pa_z[i];

        tr_y_xyz_xxz[i] = 2.0 * tr_y_yz_xz[i] * fe_0 + tr_y_yz_xxz[i] * pa_x[i];

        tr_y_xyz_xyy[i] = tr_y_xy_xyy[i] * pa_z[i];

        tr_y_xyz_xyz[i] = tr_y_yz_yz[i] * fe_0 + tr_y_yz_xyz[i] * pa_x[i];

        tr_y_xyz_xzz[i] = tr_y_yz_zz[i] * fe_0 + tr_y_yz_xzz[i] * pa_x[i];

        tr_y_xyz_yyy[i] = tr_y_yz_yyy[i] * pa_x[i];

        tr_y_xyz_yyz[i] = tr_y_yz_yyz[i] * pa_x[i];

        tr_y_xyz_yzz[i] = tr_y_yz_yzz[i] * pa_x[i];

        tr_y_xyz_zzz[i] = tr_y_yz_zzz[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : FF

    auto tr_y_xzz_xxx = pbuffer.data(idx_dip_ff + 150);

    auto tr_y_xzz_xxy = pbuffer.data(idx_dip_ff + 151);

    auto tr_y_xzz_xxz = pbuffer.data(idx_dip_ff + 152);

    auto tr_y_xzz_xyy = pbuffer.data(idx_dip_ff + 153);

    auto tr_y_xzz_xyz = pbuffer.data(idx_dip_ff + 154);

    auto tr_y_xzz_xzz = pbuffer.data(idx_dip_ff + 155);

    auto tr_y_xzz_yyy = pbuffer.data(idx_dip_ff + 156);

    auto tr_y_xzz_yyz = pbuffer.data(idx_dip_ff + 157);

    auto tr_y_xzz_yzz = pbuffer.data(idx_dip_ff + 158);

    auto tr_y_xzz_zzz = pbuffer.data(idx_dip_ff + 159);

    #pragma omp simd aligned(pa_x, tr_y_xzz_xxx, tr_y_xzz_xxy, tr_y_xzz_xxz, tr_y_xzz_xyy, tr_y_xzz_xyz, tr_y_xzz_xzz, tr_y_xzz_yyy, tr_y_xzz_yyz, tr_y_xzz_yzz, tr_y_xzz_zzz, tr_y_zz_xx, tr_y_zz_xxx, tr_y_zz_xxy, tr_y_zz_xxz, tr_y_zz_xy, tr_y_zz_xyy, tr_y_zz_xyz, tr_y_zz_xz, tr_y_zz_xzz, tr_y_zz_yy, tr_y_zz_yyy, tr_y_zz_yyz, tr_y_zz_yz, tr_y_zz_yzz, tr_y_zz_zz, tr_y_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzz_xxx[i] = 3.0 * tr_y_zz_xx[i] * fe_0 + tr_y_zz_xxx[i] * pa_x[i];

        tr_y_xzz_xxy[i] = 2.0 * tr_y_zz_xy[i] * fe_0 + tr_y_zz_xxy[i] * pa_x[i];

        tr_y_xzz_xxz[i] = 2.0 * tr_y_zz_xz[i] * fe_0 + tr_y_zz_xxz[i] * pa_x[i];

        tr_y_xzz_xyy[i] = tr_y_zz_yy[i] * fe_0 + tr_y_zz_xyy[i] * pa_x[i];

        tr_y_xzz_xyz[i] = tr_y_zz_yz[i] * fe_0 + tr_y_zz_xyz[i] * pa_x[i];

        tr_y_xzz_xzz[i] = tr_y_zz_zz[i] * fe_0 + tr_y_zz_xzz[i] * pa_x[i];

        tr_y_xzz_yyy[i] = tr_y_zz_yyy[i] * pa_x[i];

        tr_y_xzz_yyz[i] = tr_y_zz_yyz[i] * pa_x[i];

        tr_y_xzz_yzz[i] = tr_y_zz_yzz[i] * pa_x[i];

        tr_y_xzz_zzz[i] = tr_y_zz_zzz[i] * pa_x[i];
    }

    // Set up 160-170 components of targeted buffer : FF

    auto tr_y_yyy_xxx = pbuffer.data(idx_dip_ff + 160);

    auto tr_y_yyy_xxy = pbuffer.data(idx_dip_ff + 161);

    auto tr_y_yyy_xxz = pbuffer.data(idx_dip_ff + 162);

    auto tr_y_yyy_xyy = pbuffer.data(idx_dip_ff + 163);

    auto tr_y_yyy_xyz = pbuffer.data(idx_dip_ff + 164);

    auto tr_y_yyy_xzz = pbuffer.data(idx_dip_ff + 165);

    auto tr_y_yyy_yyy = pbuffer.data(idx_dip_ff + 166);

    auto tr_y_yyy_yyz = pbuffer.data(idx_dip_ff + 167);

    auto tr_y_yyy_yzz = pbuffer.data(idx_dip_ff + 168);

    auto tr_y_yyy_zzz = pbuffer.data(idx_dip_ff + 169);

    #pragma omp simd aligned(pa_y, tr_y_y_xxx, tr_y_y_xxy, tr_y_y_xxz, tr_y_y_xyy, tr_y_y_xyz, tr_y_y_xzz, tr_y_y_yyy, tr_y_y_yyz, tr_y_y_yzz, tr_y_y_zzz, tr_y_yy_xx, tr_y_yy_xxx, tr_y_yy_xxy, tr_y_yy_xxz, tr_y_yy_xy, tr_y_yy_xyy, tr_y_yy_xyz, tr_y_yy_xz, tr_y_yy_xzz, tr_y_yy_yy, tr_y_yy_yyy, tr_y_yy_yyz, tr_y_yy_yz, tr_y_yy_yzz, tr_y_yy_zz, tr_y_yy_zzz, tr_y_yyy_xxx, tr_y_yyy_xxy, tr_y_yyy_xxz, tr_y_yyy_xyy, tr_y_yyy_xyz, tr_y_yyy_xzz, tr_y_yyy_yyy, tr_y_yyy_yyz, tr_y_yyy_yzz, tr_y_yyy_zzz, ts_yy_xxx, ts_yy_xxy, ts_yy_xxz, ts_yy_xyy, ts_yy_xyz, ts_yy_xzz, ts_yy_yyy, ts_yy_yyz, ts_yy_yzz, ts_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyy_xxx[i] = 2.0 * tr_y_y_xxx[i] * fe_0 + ts_yy_xxx[i] * fe_0 + tr_y_yy_xxx[i] * pa_y[i];

        tr_y_yyy_xxy[i] = 2.0 * tr_y_y_xxy[i] * fe_0 + tr_y_yy_xx[i] * fe_0 + ts_yy_xxy[i] * fe_0 + tr_y_yy_xxy[i] * pa_y[i];

        tr_y_yyy_xxz[i] = 2.0 * tr_y_y_xxz[i] * fe_0 + ts_yy_xxz[i] * fe_0 + tr_y_yy_xxz[i] * pa_y[i];

        tr_y_yyy_xyy[i] = 2.0 * tr_y_y_xyy[i] * fe_0 + 2.0 * tr_y_yy_xy[i] * fe_0 + ts_yy_xyy[i] * fe_0 + tr_y_yy_xyy[i] * pa_y[i];

        tr_y_yyy_xyz[i] = 2.0 * tr_y_y_xyz[i] * fe_0 + tr_y_yy_xz[i] * fe_0 + ts_yy_xyz[i] * fe_0 + tr_y_yy_xyz[i] * pa_y[i];

        tr_y_yyy_xzz[i] = 2.0 * tr_y_y_xzz[i] * fe_0 + ts_yy_xzz[i] * fe_0 + tr_y_yy_xzz[i] * pa_y[i];

        tr_y_yyy_yyy[i] = 2.0 * tr_y_y_yyy[i] * fe_0 + 3.0 * tr_y_yy_yy[i] * fe_0 + ts_yy_yyy[i] * fe_0 + tr_y_yy_yyy[i] * pa_y[i];

        tr_y_yyy_yyz[i] = 2.0 * tr_y_y_yyz[i] * fe_0 + 2.0 * tr_y_yy_yz[i] * fe_0 + ts_yy_yyz[i] * fe_0 + tr_y_yy_yyz[i] * pa_y[i];

        tr_y_yyy_yzz[i] = 2.0 * tr_y_y_yzz[i] * fe_0 + tr_y_yy_zz[i] * fe_0 + ts_yy_yzz[i] * fe_0 + tr_y_yy_yzz[i] * pa_y[i];

        tr_y_yyy_zzz[i] = 2.0 * tr_y_y_zzz[i] * fe_0 + ts_yy_zzz[i] * fe_0 + tr_y_yy_zzz[i] * pa_y[i];
    }

    // Set up 170-180 components of targeted buffer : FF

    auto tr_y_yyz_xxx = pbuffer.data(idx_dip_ff + 170);

    auto tr_y_yyz_xxy = pbuffer.data(idx_dip_ff + 171);

    auto tr_y_yyz_xxz = pbuffer.data(idx_dip_ff + 172);

    auto tr_y_yyz_xyy = pbuffer.data(idx_dip_ff + 173);

    auto tr_y_yyz_xyz = pbuffer.data(idx_dip_ff + 174);

    auto tr_y_yyz_xzz = pbuffer.data(idx_dip_ff + 175);

    auto tr_y_yyz_yyy = pbuffer.data(idx_dip_ff + 176);

    auto tr_y_yyz_yyz = pbuffer.data(idx_dip_ff + 177);

    auto tr_y_yyz_yzz = pbuffer.data(idx_dip_ff + 178);

    auto tr_y_yyz_zzz = pbuffer.data(idx_dip_ff + 179);

    #pragma omp simd aligned(pa_z, tr_y_yy_xx, tr_y_yy_xxx, tr_y_yy_xxy, tr_y_yy_xxz, tr_y_yy_xy, tr_y_yy_xyy, tr_y_yy_xyz, tr_y_yy_xz, tr_y_yy_xzz, tr_y_yy_yy, tr_y_yy_yyy, tr_y_yy_yyz, tr_y_yy_yz, tr_y_yy_yzz, tr_y_yy_zz, tr_y_yy_zzz, tr_y_yyz_xxx, tr_y_yyz_xxy, tr_y_yyz_xxz, tr_y_yyz_xyy, tr_y_yyz_xyz, tr_y_yyz_xzz, tr_y_yyz_yyy, tr_y_yyz_yyz, tr_y_yyz_yzz, tr_y_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyz_xxx[i] = tr_y_yy_xxx[i] * pa_z[i];

        tr_y_yyz_xxy[i] = tr_y_yy_xxy[i] * pa_z[i];

        tr_y_yyz_xxz[i] = tr_y_yy_xx[i] * fe_0 + tr_y_yy_xxz[i] * pa_z[i];

        tr_y_yyz_xyy[i] = tr_y_yy_xyy[i] * pa_z[i];

        tr_y_yyz_xyz[i] = tr_y_yy_xy[i] * fe_0 + tr_y_yy_xyz[i] * pa_z[i];

        tr_y_yyz_xzz[i] = 2.0 * tr_y_yy_xz[i] * fe_0 + tr_y_yy_xzz[i] * pa_z[i];

        tr_y_yyz_yyy[i] = tr_y_yy_yyy[i] * pa_z[i];

        tr_y_yyz_yyz[i] = tr_y_yy_yy[i] * fe_0 + tr_y_yy_yyz[i] * pa_z[i];

        tr_y_yyz_yzz[i] = 2.0 * tr_y_yy_yz[i] * fe_0 + tr_y_yy_yzz[i] * pa_z[i];

        tr_y_yyz_zzz[i] = 3.0 * tr_y_yy_zz[i] * fe_0 + tr_y_yy_zzz[i] * pa_z[i];
    }

    // Set up 180-190 components of targeted buffer : FF

    auto tr_y_yzz_xxx = pbuffer.data(idx_dip_ff + 180);

    auto tr_y_yzz_xxy = pbuffer.data(idx_dip_ff + 181);

    auto tr_y_yzz_xxz = pbuffer.data(idx_dip_ff + 182);

    auto tr_y_yzz_xyy = pbuffer.data(idx_dip_ff + 183);

    auto tr_y_yzz_xyz = pbuffer.data(idx_dip_ff + 184);

    auto tr_y_yzz_xzz = pbuffer.data(idx_dip_ff + 185);

    auto tr_y_yzz_yyy = pbuffer.data(idx_dip_ff + 186);

    auto tr_y_yzz_yyz = pbuffer.data(idx_dip_ff + 187);

    auto tr_y_yzz_yzz = pbuffer.data(idx_dip_ff + 188);

    auto tr_y_yzz_zzz = pbuffer.data(idx_dip_ff + 189);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_y_xxy, tr_y_y_xyy, tr_y_y_yyy, tr_y_yz_xxy, tr_y_yz_xyy, tr_y_yz_yyy, tr_y_yzz_xxx, tr_y_yzz_xxy, tr_y_yzz_xxz, tr_y_yzz_xyy, tr_y_yzz_xyz, tr_y_yzz_xzz, tr_y_yzz_yyy, tr_y_yzz_yyz, tr_y_yzz_yzz, tr_y_yzz_zzz, tr_y_zz_xxx, tr_y_zz_xxz, tr_y_zz_xyz, tr_y_zz_xz, tr_y_zz_xzz, tr_y_zz_yyz, tr_y_zz_yz, tr_y_zz_yzz, tr_y_zz_zz, tr_y_zz_zzz, ts_zz_xxx, ts_zz_xxz, ts_zz_xyz, ts_zz_xzz, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzz_xxx[i] = ts_zz_xxx[i] * fe_0 + tr_y_zz_xxx[i] * pa_y[i];

        tr_y_yzz_xxy[i] = tr_y_y_xxy[i] * fe_0 + tr_y_yz_xxy[i] * pa_z[i];

        tr_y_yzz_xxz[i] = ts_zz_xxz[i] * fe_0 + tr_y_zz_xxz[i] * pa_y[i];

        tr_y_yzz_xyy[i] = tr_y_y_xyy[i] * fe_0 + tr_y_yz_xyy[i] * pa_z[i];

        tr_y_yzz_xyz[i] = tr_y_zz_xz[i] * fe_0 + ts_zz_xyz[i] * fe_0 + tr_y_zz_xyz[i] * pa_y[i];

        tr_y_yzz_xzz[i] = ts_zz_xzz[i] * fe_0 + tr_y_zz_xzz[i] * pa_y[i];

        tr_y_yzz_yyy[i] = tr_y_y_yyy[i] * fe_0 + tr_y_yz_yyy[i] * pa_z[i];

        tr_y_yzz_yyz[i] = 2.0 * tr_y_zz_yz[i] * fe_0 + ts_zz_yyz[i] * fe_0 + tr_y_zz_yyz[i] * pa_y[i];

        tr_y_yzz_yzz[i] = tr_y_zz_zz[i] * fe_0 + ts_zz_yzz[i] * fe_0 + tr_y_zz_yzz[i] * pa_y[i];

        tr_y_yzz_zzz[i] = ts_zz_zzz[i] * fe_0 + tr_y_zz_zzz[i] * pa_y[i];
    }

    // Set up 190-200 components of targeted buffer : FF

    auto tr_y_zzz_xxx = pbuffer.data(idx_dip_ff + 190);

    auto tr_y_zzz_xxy = pbuffer.data(idx_dip_ff + 191);

    auto tr_y_zzz_xxz = pbuffer.data(idx_dip_ff + 192);

    auto tr_y_zzz_xyy = pbuffer.data(idx_dip_ff + 193);

    auto tr_y_zzz_xyz = pbuffer.data(idx_dip_ff + 194);

    auto tr_y_zzz_xzz = pbuffer.data(idx_dip_ff + 195);

    auto tr_y_zzz_yyy = pbuffer.data(idx_dip_ff + 196);

    auto tr_y_zzz_yyz = pbuffer.data(idx_dip_ff + 197);

    auto tr_y_zzz_yzz = pbuffer.data(idx_dip_ff + 198);

    auto tr_y_zzz_zzz = pbuffer.data(idx_dip_ff + 199);

    #pragma omp simd aligned(pa_z, tr_y_z_xxx, tr_y_z_xxy, tr_y_z_xxz, tr_y_z_xyy, tr_y_z_xyz, tr_y_z_xzz, tr_y_z_yyy, tr_y_z_yyz, tr_y_z_yzz, tr_y_z_zzz, tr_y_zz_xx, tr_y_zz_xxx, tr_y_zz_xxy, tr_y_zz_xxz, tr_y_zz_xy, tr_y_zz_xyy, tr_y_zz_xyz, tr_y_zz_xz, tr_y_zz_xzz, tr_y_zz_yy, tr_y_zz_yyy, tr_y_zz_yyz, tr_y_zz_yz, tr_y_zz_yzz, tr_y_zz_zz, tr_y_zz_zzz, tr_y_zzz_xxx, tr_y_zzz_xxy, tr_y_zzz_xxz, tr_y_zzz_xyy, tr_y_zzz_xyz, tr_y_zzz_xzz, tr_y_zzz_yyy, tr_y_zzz_yyz, tr_y_zzz_yzz, tr_y_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzz_xxx[i] = 2.0 * tr_y_z_xxx[i] * fe_0 + tr_y_zz_xxx[i] * pa_z[i];

        tr_y_zzz_xxy[i] = 2.0 * tr_y_z_xxy[i] * fe_0 + tr_y_zz_xxy[i] * pa_z[i];

        tr_y_zzz_xxz[i] = 2.0 * tr_y_z_xxz[i] * fe_0 + tr_y_zz_xx[i] * fe_0 + tr_y_zz_xxz[i] * pa_z[i];

        tr_y_zzz_xyy[i] = 2.0 * tr_y_z_xyy[i] * fe_0 + tr_y_zz_xyy[i] * pa_z[i];

        tr_y_zzz_xyz[i] = 2.0 * tr_y_z_xyz[i] * fe_0 + tr_y_zz_xy[i] * fe_0 + tr_y_zz_xyz[i] * pa_z[i];

        tr_y_zzz_xzz[i] = 2.0 * tr_y_z_xzz[i] * fe_0 + 2.0 * tr_y_zz_xz[i] * fe_0 + tr_y_zz_xzz[i] * pa_z[i];

        tr_y_zzz_yyy[i] = 2.0 * tr_y_z_yyy[i] * fe_0 + tr_y_zz_yyy[i] * pa_z[i];

        tr_y_zzz_yyz[i] = 2.0 * tr_y_z_yyz[i] * fe_0 + tr_y_zz_yy[i] * fe_0 + tr_y_zz_yyz[i] * pa_z[i];

        tr_y_zzz_yzz[i] = 2.0 * tr_y_z_yzz[i] * fe_0 + 2.0 * tr_y_zz_yz[i] * fe_0 + tr_y_zz_yzz[i] * pa_z[i];

        tr_y_zzz_zzz[i] = 2.0 * tr_y_z_zzz[i] * fe_0 + 3.0 * tr_y_zz_zz[i] * fe_0 + tr_y_zz_zzz[i] * pa_z[i];
    }

    // Set up 200-210 components of targeted buffer : FF

    auto tr_z_xxx_xxx = pbuffer.data(idx_dip_ff + 200);

    auto tr_z_xxx_xxy = pbuffer.data(idx_dip_ff + 201);

    auto tr_z_xxx_xxz = pbuffer.data(idx_dip_ff + 202);

    auto tr_z_xxx_xyy = pbuffer.data(idx_dip_ff + 203);

    auto tr_z_xxx_xyz = pbuffer.data(idx_dip_ff + 204);

    auto tr_z_xxx_xzz = pbuffer.data(idx_dip_ff + 205);

    auto tr_z_xxx_yyy = pbuffer.data(idx_dip_ff + 206);

    auto tr_z_xxx_yyz = pbuffer.data(idx_dip_ff + 207);

    auto tr_z_xxx_yzz = pbuffer.data(idx_dip_ff + 208);

    auto tr_z_xxx_zzz = pbuffer.data(idx_dip_ff + 209);

    #pragma omp simd aligned(pa_x, tr_z_x_xxx, tr_z_x_xxy, tr_z_x_xxz, tr_z_x_xyy, tr_z_x_xyz, tr_z_x_xzz, tr_z_x_yyy, tr_z_x_yyz, tr_z_x_yzz, tr_z_x_zzz, tr_z_xx_xx, tr_z_xx_xxx, tr_z_xx_xxy, tr_z_xx_xxz, tr_z_xx_xy, tr_z_xx_xyy, tr_z_xx_xyz, tr_z_xx_xz, tr_z_xx_xzz, tr_z_xx_yy, tr_z_xx_yyy, tr_z_xx_yyz, tr_z_xx_yz, tr_z_xx_yzz, tr_z_xx_zz, tr_z_xx_zzz, tr_z_xxx_xxx, tr_z_xxx_xxy, tr_z_xxx_xxz, tr_z_xxx_xyy, tr_z_xxx_xyz, tr_z_xxx_xzz, tr_z_xxx_yyy, tr_z_xxx_yyz, tr_z_xxx_yzz, tr_z_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxx_xxx[i] = 2.0 * tr_z_x_xxx[i] * fe_0 + 3.0 * tr_z_xx_xx[i] * fe_0 + tr_z_xx_xxx[i] * pa_x[i];

        tr_z_xxx_xxy[i] = 2.0 * tr_z_x_xxy[i] * fe_0 + 2.0 * tr_z_xx_xy[i] * fe_0 + tr_z_xx_xxy[i] * pa_x[i];

        tr_z_xxx_xxz[i] = 2.0 * tr_z_x_xxz[i] * fe_0 + 2.0 * tr_z_xx_xz[i] * fe_0 + tr_z_xx_xxz[i] * pa_x[i];

        tr_z_xxx_xyy[i] = 2.0 * tr_z_x_xyy[i] * fe_0 + tr_z_xx_yy[i] * fe_0 + tr_z_xx_xyy[i] * pa_x[i];

        tr_z_xxx_xyz[i] = 2.0 * tr_z_x_xyz[i] * fe_0 + tr_z_xx_yz[i] * fe_0 + tr_z_xx_xyz[i] * pa_x[i];

        tr_z_xxx_xzz[i] = 2.0 * tr_z_x_xzz[i] * fe_0 + tr_z_xx_zz[i] * fe_0 + tr_z_xx_xzz[i] * pa_x[i];

        tr_z_xxx_yyy[i] = 2.0 * tr_z_x_yyy[i] * fe_0 + tr_z_xx_yyy[i] * pa_x[i];

        tr_z_xxx_yyz[i] = 2.0 * tr_z_x_yyz[i] * fe_0 + tr_z_xx_yyz[i] * pa_x[i];

        tr_z_xxx_yzz[i] = 2.0 * tr_z_x_yzz[i] * fe_0 + tr_z_xx_yzz[i] * pa_x[i];

        tr_z_xxx_zzz[i] = 2.0 * tr_z_x_zzz[i] * fe_0 + tr_z_xx_zzz[i] * pa_x[i];
    }

    // Set up 210-220 components of targeted buffer : FF

    auto tr_z_xxy_xxx = pbuffer.data(idx_dip_ff + 210);

    auto tr_z_xxy_xxy = pbuffer.data(idx_dip_ff + 211);

    auto tr_z_xxy_xxz = pbuffer.data(idx_dip_ff + 212);

    auto tr_z_xxy_xyy = pbuffer.data(idx_dip_ff + 213);

    auto tr_z_xxy_xyz = pbuffer.data(idx_dip_ff + 214);

    auto tr_z_xxy_xzz = pbuffer.data(idx_dip_ff + 215);

    auto tr_z_xxy_yyy = pbuffer.data(idx_dip_ff + 216);

    auto tr_z_xxy_yyz = pbuffer.data(idx_dip_ff + 217);

    auto tr_z_xxy_yzz = pbuffer.data(idx_dip_ff + 218);

    auto tr_z_xxy_zzz = pbuffer.data(idx_dip_ff + 219);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xx_xx, tr_z_xx_xxx, tr_z_xx_xxy, tr_z_xx_xxz, tr_z_xx_xy, tr_z_xx_xyy, tr_z_xx_xyz, tr_z_xx_xz, tr_z_xx_xzz, tr_z_xx_zzz, tr_z_xxy_xxx, tr_z_xxy_xxy, tr_z_xxy_xxz, tr_z_xxy_xyy, tr_z_xxy_xyz, tr_z_xxy_xzz, tr_z_xxy_yyy, tr_z_xxy_yyz, tr_z_xxy_yzz, tr_z_xxy_zzz, tr_z_xy_yyy, tr_z_xy_yyz, tr_z_xy_yzz, tr_z_y_yyy, tr_z_y_yyz, tr_z_y_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxy_xxx[i] = tr_z_xx_xxx[i] * pa_y[i];

        tr_z_xxy_xxy[i] = tr_z_xx_xx[i] * fe_0 + tr_z_xx_xxy[i] * pa_y[i];

        tr_z_xxy_xxz[i] = tr_z_xx_xxz[i] * pa_y[i];

        tr_z_xxy_xyy[i] = 2.0 * tr_z_xx_xy[i] * fe_0 + tr_z_xx_xyy[i] * pa_y[i];

        tr_z_xxy_xyz[i] = tr_z_xx_xz[i] * fe_0 + tr_z_xx_xyz[i] * pa_y[i];

        tr_z_xxy_xzz[i] = tr_z_xx_xzz[i] * pa_y[i];

        tr_z_xxy_yyy[i] = tr_z_y_yyy[i] * fe_0 + tr_z_xy_yyy[i] * pa_x[i];

        tr_z_xxy_yyz[i] = tr_z_y_yyz[i] * fe_0 + tr_z_xy_yyz[i] * pa_x[i];

        tr_z_xxy_yzz[i] = tr_z_y_yzz[i] * fe_0 + tr_z_xy_yzz[i] * pa_x[i];

        tr_z_xxy_zzz[i] = tr_z_xx_zzz[i] * pa_y[i];
    }

    // Set up 220-230 components of targeted buffer : FF

    auto tr_z_xxz_xxx = pbuffer.data(idx_dip_ff + 220);

    auto tr_z_xxz_xxy = pbuffer.data(idx_dip_ff + 221);

    auto tr_z_xxz_xxz = pbuffer.data(idx_dip_ff + 222);

    auto tr_z_xxz_xyy = pbuffer.data(idx_dip_ff + 223);

    auto tr_z_xxz_xyz = pbuffer.data(idx_dip_ff + 224);

    auto tr_z_xxz_xzz = pbuffer.data(idx_dip_ff + 225);

    auto tr_z_xxz_yyy = pbuffer.data(idx_dip_ff + 226);

    auto tr_z_xxz_yyz = pbuffer.data(idx_dip_ff + 227);

    auto tr_z_xxz_yzz = pbuffer.data(idx_dip_ff + 228);

    auto tr_z_xxz_zzz = pbuffer.data(idx_dip_ff + 229);

    #pragma omp simd aligned(pa_x, pa_z, tr_z_xx_xxx, tr_z_xx_xxy, tr_z_xx_xyy, tr_z_xxz_xxx, tr_z_xxz_xxy, tr_z_xxz_xxz, tr_z_xxz_xyy, tr_z_xxz_xyz, tr_z_xxz_xzz, tr_z_xxz_yyy, tr_z_xxz_yyz, tr_z_xxz_yzz, tr_z_xxz_zzz, tr_z_xz_xxz, tr_z_xz_xyz, tr_z_xz_xz, tr_z_xz_xzz, tr_z_xz_yyy, tr_z_xz_yyz, tr_z_xz_yz, tr_z_xz_yzz, tr_z_xz_zz, tr_z_xz_zzz, tr_z_z_xxz, tr_z_z_xyz, tr_z_z_xzz, tr_z_z_yyy, tr_z_z_yyz, tr_z_z_yzz, tr_z_z_zzz, ts_xx_xxx, ts_xx_xxy, ts_xx_xyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxz_xxx[i] = ts_xx_xxx[i] * fe_0 + tr_z_xx_xxx[i] * pa_z[i];

        tr_z_xxz_xxy[i] = ts_xx_xxy[i] * fe_0 + tr_z_xx_xxy[i] * pa_z[i];

        tr_z_xxz_xxz[i] = tr_z_z_xxz[i] * fe_0 + 2.0 * tr_z_xz_xz[i] * fe_0 + tr_z_xz_xxz[i] * pa_x[i];

        tr_z_xxz_xyy[i] = ts_xx_xyy[i] * fe_0 + tr_z_xx_xyy[i] * pa_z[i];

        tr_z_xxz_xyz[i] = tr_z_z_xyz[i] * fe_0 + tr_z_xz_yz[i] * fe_0 + tr_z_xz_xyz[i] * pa_x[i];

        tr_z_xxz_xzz[i] = tr_z_z_xzz[i] * fe_0 + tr_z_xz_zz[i] * fe_0 + tr_z_xz_xzz[i] * pa_x[i];

        tr_z_xxz_yyy[i] = tr_z_z_yyy[i] * fe_0 + tr_z_xz_yyy[i] * pa_x[i];

        tr_z_xxz_yyz[i] = tr_z_z_yyz[i] * fe_0 + tr_z_xz_yyz[i] * pa_x[i];

        tr_z_xxz_yzz[i] = tr_z_z_yzz[i] * fe_0 + tr_z_xz_yzz[i] * pa_x[i];

        tr_z_xxz_zzz[i] = tr_z_z_zzz[i] * fe_0 + tr_z_xz_zzz[i] * pa_x[i];
    }

    // Set up 230-240 components of targeted buffer : FF

    auto tr_z_xyy_xxx = pbuffer.data(idx_dip_ff + 230);

    auto tr_z_xyy_xxy = pbuffer.data(idx_dip_ff + 231);

    auto tr_z_xyy_xxz = pbuffer.data(idx_dip_ff + 232);

    auto tr_z_xyy_xyy = pbuffer.data(idx_dip_ff + 233);

    auto tr_z_xyy_xyz = pbuffer.data(idx_dip_ff + 234);

    auto tr_z_xyy_xzz = pbuffer.data(idx_dip_ff + 235);

    auto tr_z_xyy_yyy = pbuffer.data(idx_dip_ff + 236);

    auto tr_z_xyy_yyz = pbuffer.data(idx_dip_ff + 237);

    auto tr_z_xyy_yzz = pbuffer.data(idx_dip_ff + 238);

    auto tr_z_xyy_zzz = pbuffer.data(idx_dip_ff + 239);

    #pragma omp simd aligned(pa_x, tr_z_xyy_xxx, tr_z_xyy_xxy, tr_z_xyy_xxz, tr_z_xyy_xyy, tr_z_xyy_xyz, tr_z_xyy_xzz, tr_z_xyy_yyy, tr_z_xyy_yyz, tr_z_xyy_yzz, tr_z_xyy_zzz, tr_z_yy_xx, tr_z_yy_xxx, tr_z_yy_xxy, tr_z_yy_xxz, tr_z_yy_xy, tr_z_yy_xyy, tr_z_yy_xyz, tr_z_yy_xz, tr_z_yy_xzz, tr_z_yy_yy, tr_z_yy_yyy, tr_z_yy_yyz, tr_z_yy_yz, tr_z_yy_yzz, tr_z_yy_zz, tr_z_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyy_xxx[i] = 3.0 * tr_z_yy_xx[i] * fe_0 + tr_z_yy_xxx[i] * pa_x[i];

        tr_z_xyy_xxy[i] = 2.0 * tr_z_yy_xy[i] * fe_0 + tr_z_yy_xxy[i] * pa_x[i];

        tr_z_xyy_xxz[i] = 2.0 * tr_z_yy_xz[i] * fe_0 + tr_z_yy_xxz[i] * pa_x[i];

        tr_z_xyy_xyy[i] = tr_z_yy_yy[i] * fe_0 + tr_z_yy_xyy[i] * pa_x[i];

        tr_z_xyy_xyz[i] = tr_z_yy_yz[i] * fe_0 + tr_z_yy_xyz[i] * pa_x[i];

        tr_z_xyy_xzz[i] = tr_z_yy_zz[i] * fe_0 + tr_z_yy_xzz[i] * pa_x[i];

        tr_z_xyy_yyy[i] = tr_z_yy_yyy[i] * pa_x[i];

        tr_z_xyy_yyz[i] = tr_z_yy_yyz[i] * pa_x[i];

        tr_z_xyy_yzz[i] = tr_z_yy_yzz[i] * pa_x[i];

        tr_z_xyy_zzz[i] = tr_z_yy_zzz[i] * pa_x[i];
    }

    // Set up 240-250 components of targeted buffer : FF

    auto tr_z_xyz_xxx = pbuffer.data(idx_dip_ff + 240);

    auto tr_z_xyz_xxy = pbuffer.data(idx_dip_ff + 241);

    auto tr_z_xyz_xxz = pbuffer.data(idx_dip_ff + 242);

    auto tr_z_xyz_xyy = pbuffer.data(idx_dip_ff + 243);

    auto tr_z_xyz_xyz = pbuffer.data(idx_dip_ff + 244);

    auto tr_z_xyz_xzz = pbuffer.data(idx_dip_ff + 245);

    auto tr_z_xyz_yyy = pbuffer.data(idx_dip_ff + 246);

    auto tr_z_xyz_yyz = pbuffer.data(idx_dip_ff + 247);

    auto tr_z_xyz_yzz = pbuffer.data(idx_dip_ff + 248);

    auto tr_z_xyz_zzz = pbuffer.data(idx_dip_ff + 249);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xyz_xxx, tr_z_xyz_xxy, tr_z_xyz_xxz, tr_z_xyz_xyy, tr_z_xyz_xyz, tr_z_xyz_xzz, tr_z_xyz_yyy, tr_z_xyz_yyz, tr_z_xyz_yzz, tr_z_xyz_zzz, tr_z_xz_xxx, tr_z_xz_xxz, tr_z_xz_xzz, tr_z_yz_xxy, tr_z_yz_xy, tr_z_yz_xyy, tr_z_yz_xyz, tr_z_yz_yy, tr_z_yz_yyy, tr_z_yz_yyz, tr_z_yz_yz, tr_z_yz_yzz, tr_z_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyz_xxx[i] = tr_z_xz_xxx[i] * pa_y[i];

        tr_z_xyz_xxy[i] = 2.0 * tr_z_yz_xy[i] * fe_0 + tr_z_yz_xxy[i] * pa_x[i];

        tr_z_xyz_xxz[i] = tr_z_xz_xxz[i] * pa_y[i];

        tr_z_xyz_xyy[i] = tr_z_yz_yy[i] * fe_0 + tr_z_yz_xyy[i] * pa_x[i];

        tr_z_xyz_xyz[i] = tr_z_yz_yz[i] * fe_0 + tr_z_yz_xyz[i] * pa_x[i];

        tr_z_xyz_xzz[i] = tr_z_xz_xzz[i] * pa_y[i];

        tr_z_xyz_yyy[i] = tr_z_yz_yyy[i] * pa_x[i];

        tr_z_xyz_yyz[i] = tr_z_yz_yyz[i] * pa_x[i];

        tr_z_xyz_yzz[i] = tr_z_yz_yzz[i] * pa_x[i];

        tr_z_xyz_zzz[i] = tr_z_yz_zzz[i] * pa_x[i];
    }

    // Set up 250-260 components of targeted buffer : FF

    auto tr_z_xzz_xxx = pbuffer.data(idx_dip_ff + 250);

    auto tr_z_xzz_xxy = pbuffer.data(idx_dip_ff + 251);

    auto tr_z_xzz_xxz = pbuffer.data(idx_dip_ff + 252);

    auto tr_z_xzz_xyy = pbuffer.data(idx_dip_ff + 253);

    auto tr_z_xzz_xyz = pbuffer.data(idx_dip_ff + 254);

    auto tr_z_xzz_xzz = pbuffer.data(idx_dip_ff + 255);

    auto tr_z_xzz_yyy = pbuffer.data(idx_dip_ff + 256);

    auto tr_z_xzz_yyz = pbuffer.data(idx_dip_ff + 257);

    auto tr_z_xzz_yzz = pbuffer.data(idx_dip_ff + 258);

    auto tr_z_xzz_zzz = pbuffer.data(idx_dip_ff + 259);

    #pragma omp simd aligned(pa_x, tr_z_xzz_xxx, tr_z_xzz_xxy, tr_z_xzz_xxz, tr_z_xzz_xyy, tr_z_xzz_xyz, tr_z_xzz_xzz, tr_z_xzz_yyy, tr_z_xzz_yyz, tr_z_xzz_yzz, tr_z_xzz_zzz, tr_z_zz_xx, tr_z_zz_xxx, tr_z_zz_xxy, tr_z_zz_xxz, tr_z_zz_xy, tr_z_zz_xyy, tr_z_zz_xyz, tr_z_zz_xz, tr_z_zz_xzz, tr_z_zz_yy, tr_z_zz_yyy, tr_z_zz_yyz, tr_z_zz_yz, tr_z_zz_yzz, tr_z_zz_zz, tr_z_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzz_xxx[i] = 3.0 * tr_z_zz_xx[i] * fe_0 + tr_z_zz_xxx[i] * pa_x[i];

        tr_z_xzz_xxy[i] = 2.0 * tr_z_zz_xy[i] * fe_0 + tr_z_zz_xxy[i] * pa_x[i];

        tr_z_xzz_xxz[i] = 2.0 * tr_z_zz_xz[i] * fe_0 + tr_z_zz_xxz[i] * pa_x[i];

        tr_z_xzz_xyy[i] = tr_z_zz_yy[i] * fe_0 + tr_z_zz_xyy[i] * pa_x[i];

        tr_z_xzz_xyz[i] = tr_z_zz_yz[i] * fe_0 + tr_z_zz_xyz[i] * pa_x[i];

        tr_z_xzz_xzz[i] = tr_z_zz_zz[i] * fe_0 + tr_z_zz_xzz[i] * pa_x[i];

        tr_z_xzz_yyy[i] = tr_z_zz_yyy[i] * pa_x[i];

        tr_z_xzz_yyz[i] = tr_z_zz_yyz[i] * pa_x[i];

        tr_z_xzz_yzz[i] = tr_z_zz_yzz[i] * pa_x[i];

        tr_z_xzz_zzz[i] = tr_z_zz_zzz[i] * pa_x[i];
    }

    // Set up 260-270 components of targeted buffer : FF

    auto tr_z_yyy_xxx = pbuffer.data(idx_dip_ff + 260);

    auto tr_z_yyy_xxy = pbuffer.data(idx_dip_ff + 261);

    auto tr_z_yyy_xxz = pbuffer.data(idx_dip_ff + 262);

    auto tr_z_yyy_xyy = pbuffer.data(idx_dip_ff + 263);

    auto tr_z_yyy_xyz = pbuffer.data(idx_dip_ff + 264);

    auto tr_z_yyy_xzz = pbuffer.data(idx_dip_ff + 265);

    auto tr_z_yyy_yyy = pbuffer.data(idx_dip_ff + 266);

    auto tr_z_yyy_yyz = pbuffer.data(idx_dip_ff + 267);

    auto tr_z_yyy_yzz = pbuffer.data(idx_dip_ff + 268);

    auto tr_z_yyy_zzz = pbuffer.data(idx_dip_ff + 269);

    #pragma omp simd aligned(pa_y, tr_z_y_xxx, tr_z_y_xxy, tr_z_y_xxz, tr_z_y_xyy, tr_z_y_xyz, tr_z_y_xzz, tr_z_y_yyy, tr_z_y_yyz, tr_z_y_yzz, tr_z_y_zzz, tr_z_yy_xx, tr_z_yy_xxx, tr_z_yy_xxy, tr_z_yy_xxz, tr_z_yy_xy, tr_z_yy_xyy, tr_z_yy_xyz, tr_z_yy_xz, tr_z_yy_xzz, tr_z_yy_yy, tr_z_yy_yyy, tr_z_yy_yyz, tr_z_yy_yz, tr_z_yy_yzz, tr_z_yy_zz, tr_z_yy_zzz, tr_z_yyy_xxx, tr_z_yyy_xxy, tr_z_yyy_xxz, tr_z_yyy_xyy, tr_z_yyy_xyz, tr_z_yyy_xzz, tr_z_yyy_yyy, tr_z_yyy_yyz, tr_z_yyy_yzz, tr_z_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyy_xxx[i] = 2.0 * tr_z_y_xxx[i] * fe_0 + tr_z_yy_xxx[i] * pa_y[i];

        tr_z_yyy_xxy[i] = 2.0 * tr_z_y_xxy[i] * fe_0 + tr_z_yy_xx[i] * fe_0 + tr_z_yy_xxy[i] * pa_y[i];

        tr_z_yyy_xxz[i] = 2.0 * tr_z_y_xxz[i] * fe_0 + tr_z_yy_xxz[i] * pa_y[i];

        tr_z_yyy_xyy[i] = 2.0 * tr_z_y_xyy[i] * fe_0 + 2.0 * tr_z_yy_xy[i] * fe_0 + tr_z_yy_xyy[i] * pa_y[i];

        tr_z_yyy_xyz[i] = 2.0 * tr_z_y_xyz[i] * fe_0 + tr_z_yy_xz[i] * fe_0 + tr_z_yy_xyz[i] * pa_y[i];

        tr_z_yyy_xzz[i] = 2.0 * tr_z_y_xzz[i] * fe_0 + tr_z_yy_xzz[i] * pa_y[i];

        tr_z_yyy_yyy[i] = 2.0 * tr_z_y_yyy[i] * fe_0 + 3.0 * tr_z_yy_yy[i] * fe_0 + tr_z_yy_yyy[i] * pa_y[i];

        tr_z_yyy_yyz[i] = 2.0 * tr_z_y_yyz[i] * fe_0 + 2.0 * tr_z_yy_yz[i] * fe_0 + tr_z_yy_yyz[i] * pa_y[i];

        tr_z_yyy_yzz[i] = 2.0 * tr_z_y_yzz[i] * fe_0 + tr_z_yy_zz[i] * fe_0 + tr_z_yy_yzz[i] * pa_y[i];

        tr_z_yyy_zzz[i] = 2.0 * tr_z_y_zzz[i] * fe_0 + tr_z_yy_zzz[i] * pa_y[i];
    }

    // Set up 270-280 components of targeted buffer : FF

    auto tr_z_yyz_xxx = pbuffer.data(idx_dip_ff + 270);

    auto tr_z_yyz_xxy = pbuffer.data(idx_dip_ff + 271);

    auto tr_z_yyz_xxz = pbuffer.data(idx_dip_ff + 272);

    auto tr_z_yyz_xyy = pbuffer.data(idx_dip_ff + 273);

    auto tr_z_yyz_xyz = pbuffer.data(idx_dip_ff + 274);

    auto tr_z_yyz_xzz = pbuffer.data(idx_dip_ff + 275);

    auto tr_z_yyz_yyy = pbuffer.data(idx_dip_ff + 276);

    auto tr_z_yyz_yyz = pbuffer.data(idx_dip_ff + 277);

    auto tr_z_yyz_yzz = pbuffer.data(idx_dip_ff + 278);

    auto tr_z_yyz_zzz = pbuffer.data(idx_dip_ff + 279);

    #pragma omp simd aligned(pa_y, pa_z, tr_z_yy_xxy, tr_z_yy_xyy, tr_z_yy_yyy, tr_z_yyz_xxx, tr_z_yyz_xxy, tr_z_yyz_xxz, tr_z_yyz_xyy, tr_z_yyz_xyz, tr_z_yyz_xzz, tr_z_yyz_yyy, tr_z_yyz_yyz, tr_z_yyz_yzz, tr_z_yyz_zzz, tr_z_yz_xxx, tr_z_yz_xxz, tr_z_yz_xyz, tr_z_yz_xz, tr_z_yz_xzz, tr_z_yz_yyz, tr_z_yz_yz, tr_z_yz_yzz, tr_z_yz_zz, tr_z_yz_zzz, tr_z_z_xxx, tr_z_z_xxz, tr_z_z_xyz, tr_z_z_xzz, tr_z_z_yyz, tr_z_z_yzz, tr_z_z_zzz, ts_yy_xxy, ts_yy_xyy, ts_yy_yyy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyz_xxx[i] = tr_z_z_xxx[i] * fe_0 + tr_z_yz_xxx[i] * pa_y[i];

        tr_z_yyz_xxy[i] = ts_yy_xxy[i] * fe_0 + tr_z_yy_xxy[i] * pa_z[i];

        tr_z_yyz_xxz[i] = tr_z_z_xxz[i] * fe_0 + tr_z_yz_xxz[i] * pa_y[i];

        tr_z_yyz_xyy[i] = ts_yy_xyy[i] * fe_0 + tr_z_yy_xyy[i] * pa_z[i];

        tr_z_yyz_xyz[i] = tr_z_z_xyz[i] * fe_0 + tr_z_yz_xz[i] * fe_0 + tr_z_yz_xyz[i] * pa_y[i];

        tr_z_yyz_xzz[i] = tr_z_z_xzz[i] * fe_0 + tr_z_yz_xzz[i] * pa_y[i];

        tr_z_yyz_yyy[i] = ts_yy_yyy[i] * fe_0 + tr_z_yy_yyy[i] * pa_z[i];

        tr_z_yyz_yyz[i] = tr_z_z_yyz[i] * fe_0 + 2.0 * tr_z_yz_yz[i] * fe_0 + tr_z_yz_yyz[i] * pa_y[i];

        tr_z_yyz_yzz[i] = tr_z_z_yzz[i] * fe_0 + tr_z_yz_zz[i] * fe_0 + tr_z_yz_yzz[i] * pa_y[i];

        tr_z_yyz_zzz[i] = tr_z_z_zzz[i] * fe_0 + tr_z_yz_zzz[i] * pa_y[i];
    }

    // Set up 280-290 components of targeted buffer : FF

    auto tr_z_yzz_xxx = pbuffer.data(idx_dip_ff + 280);

    auto tr_z_yzz_xxy = pbuffer.data(idx_dip_ff + 281);

    auto tr_z_yzz_xxz = pbuffer.data(idx_dip_ff + 282);

    auto tr_z_yzz_xyy = pbuffer.data(idx_dip_ff + 283);

    auto tr_z_yzz_xyz = pbuffer.data(idx_dip_ff + 284);

    auto tr_z_yzz_xzz = pbuffer.data(idx_dip_ff + 285);

    auto tr_z_yzz_yyy = pbuffer.data(idx_dip_ff + 286);

    auto tr_z_yzz_yyz = pbuffer.data(idx_dip_ff + 287);

    auto tr_z_yzz_yzz = pbuffer.data(idx_dip_ff + 288);

    auto tr_z_yzz_zzz = pbuffer.data(idx_dip_ff + 289);

    #pragma omp simd aligned(pa_y, tr_z_yzz_xxx, tr_z_yzz_xxy, tr_z_yzz_xxz, tr_z_yzz_xyy, tr_z_yzz_xyz, tr_z_yzz_xzz, tr_z_yzz_yyy, tr_z_yzz_yyz, tr_z_yzz_yzz, tr_z_yzz_zzz, tr_z_zz_xx, tr_z_zz_xxx, tr_z_zz_xxy, tr_z_zz_xxz, tr_z_zz_xy, tr_z_zz_xyy, tr_z_zz_xyz, tr_z_zz_xz, tr_z_zz_xzz, tr_z_zz_yy, tr_z_zz_yyy, tr_z_zz_yyz, tr_z_zz_yz, tr_z_zz_yzz, tr_z_zz_zz, tr_z_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzz_xxx[i] = tr_z_zz_xxx[i] * pa_y[i];

        tr_z_yzz_xxy[i] = tr_z_zz_xx[i] * fe_0 + tr_z_zz_xxy[i] * pa_y[i];

        tr_z_yzz_xxz[i] = tr_z_zz_xxz[i] * pa_y[i];

        tr_z_yzz_xyy[i] = 2.0 * tr_z_zz_xy[i] * fe_0 + tr_z_zz_xyy[i] * pa_y[i];

        tr_z_yzz_xyz[i] = tr_z_zz_xz[i] * fe_0 + tr_z_zz_xyz[i] * pa_y[i];

        tr_z_yzz_xzz[i] = tr_z_zz_xzz[i] * pa_y[i];

        tr_z_yzz_yyy[i] = 3.0 * tr_z_zz_yy[i] * fe_0 + tr_z_zz_yyy[i] * pa_y[i];

        tr_z_yzz_yyz[i] = 2.0 * tr_z_zz_yz[i] * fe_0 + tr_z_zz_yyz[i] * pa_y[i];

        tr_z_yzz_yzz[i] = tr_z_zz_zz[i] * fe_0 + tr_z_zz_yzz[i] * pa_y[i];

        tr_z_yzz_zzz[i] = tr_z_zz_zzz[i] * pa_y[i];
    }

    // Set up 290-300 components of targeted buffer : FF

    auto tr_z_zzz_xxx = pbuffer.data(idx_dip_ff + 290);

    auto tr_z_zzz_xxy = pbuffer.data(idx_dip_ff + 291);

    auto tr_z_zzz_xxz = pbuffer.data(idx_dip_ff + 292);

    auto tr_z_zzz_xyy = pbuffer.data(idx_dip_ff + 293);

    auto tr_z_zzz_xyz = pbuffer.data(idx_dip_ff + 294);

    auto tr_z_zzz_xzz = pbuffer.data(idx_dip_ff + 295);

    auto tr_z_zzz_yyy = pbuffer.data(idx_dip_ff + 296);

    auto tr_z_zzz_yyz = pbuffer.data(idx_dip_ff + 297);

    auto tr_z_zzz_yzz = pbuffer.data(idx_dip_ff + 298);

    auto tr_z_zzz_zzz = pbuffer.data(idx_dip_ff + 299);

    #pragma omp simd aligned(pa_z, tr_z_z_xxx, tr_z_z_xxy, tr_z_z_xxz, tr_z_z_xyy, tr_z_z_xyz, tr_z_z_xzz, tr_z_z_yyy, tr_z_z_yyz, tr_z_z_yzz, tr_z_z_zzz, tr_z_zz_xx, tr_z_zz_xxx, tr_z_zz_xxy, tr_z_zz_xxz, tr_z_zz_xy, tr_z_zz_xyy, tr_z_zz_xyz, tr_z_zz_xz, tr_z_zz_xzz, tr_z_zz_yy, tr_z_zz_yyy, tr_z_zz_yyz, tr_z_zz_yz, tr_z_zz_yzz, tr_z_zz_zz, tr_z_zz_zzz, tr_z_zzz_xxx, tr_z_zzz_xxy, tr_z_zzz_xxz, tr_z_zzz_xyy, tr_z_zzz_xyz, tr_z_zzz_xzz, tr_z_zzz_yyy, tr_z_zzz_yyz, tr_z_zzz_yzz, tr_z_zzz_zzz, ts_zz_xxx, ts_zz_xxy, ts_zz_xxz, ts_zz_xyy, ts_zz_xyz, ts_zz_xzz, ts_zz_yyy, ts_zz_yyz, ts_zz_yzz, ts_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzz_xxx[i] = 2.0 * tr_z_z_xxx[i] * fe_0 + ts_zz_xxx[i] * fe_0 + tr_z_zz_xxx[i] * pa_z[i];

        tr_z_zzz_xxy[i] = 2.0 * tr_z_z_xxy[i] * fe_0 + ts_zz_xxy[i] * fe_0 + tr_z_zz_xxy[i] * pa_z[i];

        tr_z_zzz_xxz[i] = 2.0 * tr_z_z_xxz[i] * fe_0 + tr_z_zz_xx[i] * fe_0 + ts_zz_xxz[i] * fe_0 + tr_z_zz_xxz[i] * pa_z[i];

        tr_z_zzz_xyy[i] = 2.0 * tr_z_z_xyy[i] * fe_0 + ts_zz_xyy[i] * fe_0 + tr_z_zz_xyy[i] * pa_z[i];

        tr_z_zzz_xyz[i] = 2.0 * tr_z_z_xyz[i] * fe_0 + tr_z_zz_xy[i] * fe_0 + ts_zz_xyz[i] * fe_0 + tr_z_zz_xyz[i] * pa_z[i];

        tr_z_zzz_xzz[i] = 2.0 * tr_z_z_xzz[i] * fe_0 + 2.0 * tr_z_zz_xz[i] * fe_0 + ts_zz_xzz[i] * fe_0 + tr_z_zz_xzz[i] * pa_z[i];

        tr_z_zzz_yyy[i] = 2.0 * tr_z_z_yyy[i] * fe_0 + ts_zz_yyy[i] * fe_0 + tr_z_zz_yyy[i] * pa_z[i];

        tr_z_zzz_yyz[i] = 2.0 * tr_z_z_yyz[i] * fe_0 + tr_z_zz_yy[i] * fe_0 + ts_zz_yyz[i] * fe_0 + tr_z_zz_yyz[i] * pa_z[i];

        tr_z_zzz_yzz[i] = 2.0 * tr_z_z_yzz[i] * fe_0 + 2.0 * tr_z_zz_yz[i] * fe_0 + ts_zz_yzz[i] * fe_0 + tr_z_zz_yzz[i] * pa_z[i];

        tr_z_zzz_zzz[i] = 2.0 * tr_z_z_zzz[i] * fe_0 + 3.0 * tr_z_zz_zz[i] * fe_0 + ts_zz_zzz[i] * fe_0 + tr_z_zz_zzz[i] * pa_z[i];
    }

}

} // diprec namespace

