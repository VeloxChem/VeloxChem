#include "TwoCenterElectronRepulsionPrimRecHF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_hf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_hf,
                                const size_t idx_eri_0_ff,
                                const size_t idx_eri_1_ff,
                                const size_t idx_eri_1_gd,
                                const size_t idx_eri_1_gf,
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

    // Set up components of auxiliary buffer : FF

    auto g_xxx_xxx_0 = pbuffer.data(idx_eri_0_ff);

    auto g_xxx_xxy_0 = pbuffer.data(idx_eri_0_ff + 1);

    auto g_xxx_xxz_0 = pbuffer.data(idx_eri_0_ff + 2);

    auto g_xxx_xyy_0 = pbuffer.data(idx_eri_0_ff + 3);

    auto g_xxx_xyz_0 = pbuffer.data(idx_eri_0_ff + 4);

    auto g_xxx_xzz_0 = pbuffer.data(idx_eri_0_ff + 5);

    auto g_xxx_yyy_0 = pbuffer.data(idx_eri_0_ff + 6);

    auto g_xxx_yyz_0 = pbuffer.data(idx_eri_0_ff + 7);

    auto g_xxx_yzz_0 = pbuffer.data(idx_eri_0_ff + 8);

    auto g_xxx_zzz_0 = pbuffer.data(idx_eri_0_ff + 9);

    auto g_xxy_xxx_0 = pbuffer.data(idx_eri_0_ff + 10);

    auto g_xxy_xxz_0 = pbuffer.data(idx_eri_0_ff + 12);

    auto g_xxy_xzz_0 = pbuffer.data(idx_eri_0_ff + 15);

    auto g_xxz_xxx_0 = pbuffer.data(idx_eri_0_ff + 20);

    auto g_xxz_xxy_0 = pbuffer.data(idx_eri_0_ff + 21);

    auto g_xxz_xyy_0 = pbuffer.data(idx_eri_0_ff + 23);

    auto g_xyy_xxy_0 = pbuffer.data(idx_eri_0_ff + 31);

    auto g_xyy_xyy_0 = pbuffer.data(idx_eri_0_ff + 33);

    auto g_xyy_xyz_0 = pbuffer.data(idx_eri_0_ff + 34);

    auto g_xyy_yyy_0 = pbuffer.data(idx_eri_0_ff + 36);

    auto g_xyy_yyz_0 = pbuffer.data(idx_eri_0_ff + 37);

    auto g_xyy_yzz_0 = pbuffer.data(idx_eri_0_ff + 38);

    auto g_xyy_zzz_0 = pbuffer.data(idx_eri_0_ff + 39);

    auto g_xzz_xxz_0 = pbuffer.data(idx_eri_0_ff + 52);

    auto g_xzz_xyz_0 = pbuffer.data(idx_eri_0_ff + 54);

    auto g_xzz_xzz_0 = pbuffer.data(idx_eri_0_ff + 55);

    auto g_xzz_yyy_0 = pbuffer.data(idx_eri_0_ff + 56);

    auto g_xzz_yyz_0 = pbuffer.data(idx_eri_0_ff + 57);

    auto g_xzz_yzz_0 = pbuffer.data(idx_eri_0_ff + 58);

    auto g_xzz_zzz_0 = pbuffer.data(idx_eri_0_ff + 59);

    auto g_yyy_xxx_0 = pbuffer.data(idx_eri_0_ff + 60);

    auto g_yyy_xxy_0 = pbuffer.data(idx_eri_0_ff + 61);

    auto g_yyy_xxz_0 = pbuffer.data(idx_eri_0_ff + 62);

    auto g_yyy_xyy_0 = pbuffer.data(idx_eri_0_ff + 63);

    auto g_yyy_xyz_0 = pbuffer.data(idx_eri_0_ff + 64);

    auto g_yyy_xzz_0 = pbuffer.data(idx_eri_0_ff + 65);

    auto g_yyy_yyy_0 = pbuffer.data(idx_eri_0_ff + 66);

    auto g_yyy_yyz_0 = pbuffer.data(idx_eri_0_ff + 67);

    auto g_yyy_yzz_0 = pbuffer.data(idx_eri_0_ff + 68);

    auto g_yyy_zzz_0 = pbuffer.data(idx_eri_0_ff + 69);

    auto g_yyz_xxy_0 = pbuffer.data(idx_eri_0_ff + 71);

    auto g_yyz_xyy_0 = pbuffer.data(idx_eri_0_ff + 73);

    auto g_yyz_yyy_0 = pbuffer.data(idx_eri_0_ff + 76);

    auto g_yzz_xxx_0 = pbuffer.data(idx_eri_0_ff + 80);

    auto g_yzz_xxz_0 = pbuffer.data(idx_eri_0_ff + 82);

    auto g_yzz_xyz_0 = pbuffer.data(idx_eri_0_ff + 84);

    auto g_yzz_xzz_0 = pbuffer.data(idx_eri_0_ff + 85);

    auto g_yzz_yyz_0 = pbuffer.data(idx_eri_0_ff + 87);

    auto g_yzz_yzz_0 = pbuffer.data(idx_eri_0_ff + 88);

    auto g_yzz_zzz_0 = pbuffer.data(idx_eri_0_ff + 89);

    auto g_zzz_xxx_0 = pbuffer.data(idx_eri_0_ff + 90);

    auto g_zzz_xxy_0 = pbuffer.data(idx_eri_0_ff + 91);

    auto g_zzz_xxz_0 = pbuffer.data(idx_eri_0_ff + 92);

    auto g_zzz_xyy_0 = pbuffer.data(idx_eri_0_ff + 93);

    auto g_zzz_xyz_0 = pbuffer.data(idx_eri_0_ff + 94);

    auto g_zzz_xzz_0 = pbuffer.data(idx_eri_0_ff + 95);

    auto g_zzz_yyy_0 = pbuffer.data(idx_eri_0_ff + 96);

    auto g_zzz_yyz_0 = pbuffer.data(idx_eri_0_ff + 97);

    auto g_zzz_yzz_0 = pbuffer.data(idx_eri_0_ff + 98);

    auto g_zzz_zzz_0 = pbuffer.data(idx_eri_0_ff + 99);

    // Set up components of auxiliary buffer : FF

    auto g_xxx_xxx_1 = pbuffer.data(idx_eri_1_ff);

    auto g_xxx_xxy_1 = pbuffer.data(idx_eri_1_ff + 1);

    auto g_xxx_xxz_1 = pbuffer.data(idx_eri_1_ff + 2);

    auto g_xxx_xyy_1 = pbuffer.data(idx_eri_1_ff + 3);

    auto g_xxx_xyz_1 = pbuffer.data(idx_eri_1_ff + 4);

    auto g_xxx_xzz_1 = pbuffer.data(idx_eri_1_ff + 5);

    auto g_xxx_yyy_1 = pbuffer.data(idx_eri_1_ff + 6);

    auto g_xxx_yyz_1 = pbuffer.data(idx_eri_1_ff + 7);

    auto g_xxx_yzz_1 = pbuffer.data(idx_eri_1_ff + 8);

    auto g_xxx_zzz_1 = pbuffer.data(idx_eri_1_ff + 9);

    auto g_xxy_xxx_1 = pbuffer.data(idx_eri_1_ff + 10);

    auto g_xxy_xxz_1 = pbuffer.data(idx_eri_1_ff + 12);

    auto g_xxy_xzz_1 = pbuffer.data(idx_eri_1_ff + 15);

    auto g_xxz_xxx_1 = pbuffer.data(idx_eri_1_ff + 20);

    auto g_xxz_xxy_1 = pbuffer.data(idx_eri_1_ff + 21);

    auto g_xxz_xyy_1 = pbuffer.data(idx_eri_1_ff + 23);

    auto g_xyy_xxy_1 = pbuffer.data(idx_eri_1_ff + 31);

    auto g_xyy_xyy_1 = pbuffer.data(idx_eri_1_ff + 33);

    auto g_xyy_xyz_1 = pbuffer.data(idx_eri_1_ff + 34);

    auto g_xyy_yyy_1 = pbuffer.data(idx_eri_1_ff + 36);

    auto g_xyy_yyz_1 = pbuffer.data(idx_eri_1_ff + 37);

    auto g_xyy_yzz_1 = pbuffer.data(idx_eri_1_ff + 38);

    auto g_xyy_zzz_1 = pbuffer.data(idx_eri_1_ff + 39);

    auto g_xzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 52);

    auto g_xzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 54);

    auto g_xzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 55);

    auto g_xzz_yyy_1 = pbuffer.data(idx_eri_1_ff + 56);

    auto g_xzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 57);

    auto g_xzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 58);

    auto g_xzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 59);

    auto g_yyy_xxx_1 = pbuffer.data(idx_eri_1_ff + 60);

    auto g_yyy_xxy_1 = pbuffer.data(idx_eri_1_ff + 61);

    auto g_yyy_xxz_1 = pbuffer.data(idx_eri_1_ff + 62);

    auto g_yyy_xyy_1 = pbuffer.data(idx_eri_1_ff + 63);

    auto g_yyy_xyz_1 = pbuffer.data(idx_eri_1_ff + 64);

    auto g_yyy_xzz_1 = pbuffer.data(idx_eri_1_ff + 65);

    auto g_yyy_yyy_1 = pbuffer.data(idx_eri_1_ff + 66);

    auto g_yyy_yyz_1 = pbuffer.data(idx_eri_1_ff + 67);

    auto g_yyy_yzz_1 = pbuffer.data(idx_eri_1_ff + 68);

    auto g_yyy_zzz_1 = pbuffer.data(idx_eri_1_ff + 69);

    auto g_yyz_xxy_1 = pbuffer.data(idx_eri_1_ff + 71);

    auto g_yyz_xyy_1 = pbuffer.data(idx_eri_1_ff + 73);

    auto g_yyz_yyy_1 = pbuffer.data(idx_eri_1_ff + 76);

    auto g_yzz_xxx_1 = pbuffer.data(idx_eri_1_ff + 80);

    auto g_yzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 82);

    auto g_yzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 84);

    auto g_yzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 85);

    auto g_yzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 87);

    auto g_yzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 88);

    auto g_yzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 89);

    auto g_zzz_xxx_1 = pbuffer.data(idx_eri_1_ff + 90);

    auto g_zzz_xxy_1 = pbuffer.data(idx_eri_1_ff + 91);

    auto g_zzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 92);

    auto g_zzz_xyy_1 = pbuffer.data(idx_eri_1_ff + 93);

    auto g_zzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 94);

    auto g_zzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 95);

    auto g_zzz_yyy_1 = pbuffer.data(idx_eri_1_ff + 96);

    auto g_zzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 97);

    auto g_zzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 98);

    auto g_zzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 99);

    // Set up components of auxiliary buffer : GD

    auto g_xxxx_xx_1 = pbuffer.data(idx_eri_1_gd);

    auto g_xxxx_xy_1 = pbuffer.data(idx_eri_1_gd + 1);

    auto g_xxxx_xz_1 = pbuffer.data(idx_eri_1_gd + 2);

    auto g_xxxx_yy_1 = pbuffer.data(idx_eri_1_gd + 3);

    auto g_xxxx_yz_1 = pbuffer.data(idx_eri_1_gd + 4);

    auto g_xxxx_zz_1 = pbuffer.data(idx_eri_1_gd + 5);

    auto g_xxxz_xz_1 = pbuffer.data(idx_eri_1_gd + 14);

    auto g_xxxz_yz_1 = pbuffer.data(idx_eri_1_gd + 16);

    auto g_xxxz_zz_1 = pbuffer.data(idx_eri_1_gd + 17);

    auto g_xxyy_xx_1 = pbuffer.data(idx_eri_1_gd + 18);

    auto g_xxyy_xy_1 = pbuffer.data(idx_eri_1_gd + 19);

    auto g_xxyy_xz_1 = pbuffer.data(idx_eri_1_gd + 20);

    auto g_xxyy_yy_1 = pbuffer.data(idx_eri_1_gd + 21);

    auto g_xxyy_yz_1 = pbuffer.data(idx_eri_1_gd + 22);

    auto g_xxyy_zz_1 = pbuffer.data(idx_eri_1_gd + 23);

    auto g_xxzz_xx_1 = pbuffer.data(idx_eri_1_gd + 30);

    auto g_xxzz_xy_1 = pbuffer.data(idx_eri_1_gd + 31);

    auto g_xxzz_xz_1 = pbuffer.data(idx_eri_1_gd + 32);

    auto g_xxzz_yy_1 = pbuffer.data(idx_eri_1_gd + 33);

    auto g_xxzz_yz_1 = pbuffer.data(idx_eri_1_gd + 34);

    auto g_xxzz_zz_1 = pbuffer.data(idx_eri_1_gd + 35);

    auto g_xyyy_xy_1 = pbuffer.data(idx_eri_1_gd + 37);

    auto g_xyyy_yy_1 = pbuffer.data(idx_eri_1_gd + 39);

    auto g_xyyy_yz_1 = pbuffer.data(idx_eri_1_gd + 40);

    auto g_xzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 56);

    auto g_xzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 58);

    auto g_xzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 59);

    auto g_yyyy_xx_1 = pbuffer.data(idx_eri_1_gd + 60);

    auto g_yyyy_xy_1 = pbuffer.data(idx_eri_1_gd + 61);

    auto g_yyyy_xz_1 = pbuffer.data(idx_eri_1_gd + 62);

    auto g_yyyy_yy_1 = pbuffer.data(idx_eri_1_gd + 63);

    auto g_yyyy_yz_1 = pbuffer.data(idx_eri_1_gd + 64);

    auto g_yyyy_zz_1 = pbuffer.data(idx_eri_1_gd + 65);

    auto g_yyyz_xz_1 = pbuffer.data(idx_eri_1_gd + 68);

    auto g_yyyz_yz_1 = pbuffer.data(idx_eri_1_gd + 70);

    auto g_yyyz_zz_1 = pbuffer.data(idx_eri_1_gd + 71);

    auto g_yyzz_xx_1 = pbuffer.data(idx_eri_1_gd + 72);

    auto g_yyzz_xy_1 = pbuffer.data(idx_eri_1_gd + 73);

    auto g_yyzz_xz_1 = pbuffer.data(idx_eri_1_gd + 74);

    auto g_yyzz_yy_1 = pbuffer.data(idx_eri_1_gd + 75);

    auto g_yyzz_yz_1 = pbuffer.data(idx_eri_1_gd + 76);

    auto g_yyzz_zz_1 = pbuffer.data(idx_eri_1_gd + 77);

    auto g_yzzz_xy_1 = pbuffer.data(idx_eri_1_gd + 79);

    auto g_yzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 80);

    auto g_yzzz_yy_1 = pbuffer.data(idx_eri_1_gd + 81);

    auto g_yzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 82);

    auto g_yzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 83);

    auto g_zzzz_xx_1 = pbuffer.data(idx_eri_1_gd + 84);

    auto g_zzzz_xy_1 = pbuffer.data(idx_eri_1_gd + 85);

    auto g_zzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 86);

    auto g_zzzz_yy_1 = pbuffer.data(idx_eri_1_gd + 87);

    auto g_zzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 88);

    auto g_zzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 89);

    // Set up components of auxiliary buffer : GF

    auto g_xxxx_xxx_1 = pbuffer.data(idx_eri_1_gf);

    auto g_xxxx_xxy_1 = pbuffer.data(idx_eri_1_gf + 1);

    auto g_xxxx_xxz_1 = pbuffer.data(idx_eri_1_gf + 2);

    auto g_xxxx_xyy_1 = pbuffer.data(idx_eri_1_gf + 3);

    auto g_xxxx_xyz_1 = pbuffer.data(idx_eri_1_gf + 4);

    auto g_xxxx_xzz_1 = pbuffer.data(idx_eri_1_gf + 5);

    auto g_xxxx_yyy_1 = pbuffer.data(idx_eri_1_gf + 6);

    auto g_xxxx_yyz_1 = pbuffer.data(idx_eri_1_gf + 7);

    auto g_xxxx_yzz_1 = pbuffer.data(idx_eri_1_gf + 8);

    auto g_xxxx_zzz_1 = pbuffer.data(idx_eri_1_gf + 9);

    auto g_xxxy_xxx_1 = pbuffer.data(idx_eri_1_gf + 10);

    auto g_xxxy_xxy_1 = pbuffer.data(idx_eri_1_gf + 11);

    auto g_xxxy_xxz_1 = pbuffer.data(idx_eri_1_gf + 12);

    auto g_xxxy_xyy_1 = pbuffer.data(idx_eri_1_gf + 13);

    auto g_xxxy_xzz_1 = pbuffer.data(idx_eri_1_gf + 15);

    auto g_xxxy_yyy_1 = pbuffer.data(idx_eri_1_gf + 16);

    auto g_xxxz_xxx_1 = pbuffer.data(idx_eri_1_gf + 20);

    auto g_xxxz_xxy_1 = pbuffer.data(idx_eri_1_gf + 21);

    auto g_xxxz_xxz_1 = pbuffer.data(idx_eri_1_gf + 22);

    auto g_xxxz_xyy_1 = pbuffer.data(idx_eri_1_gf + 23);

    auto g_xxxz_xyz_1 = pbuffer.data(idx_eri_1_gf + 24);

    auto g_xxxz_xzz_1 = pbuffer.data(idx_eri_1_gf + 25);

    auto g_xxxz_yyz_1 = pbuffer.data(idx_eri_1_gf + 27);

    auto g_xxxz_yzz_1 = pbuffer.data(idx_eri_1_gf + 28);

    auto g_xxxz_zzz_1 = pbuffer.data(idx_eri_1_gf + 29);

    auto g_xxyy_xxx_1 = pbuffer.data(idx_eri_1_gf + 30);

    auto g_xxyy_xxy_1 = pbuffer.data(idx_eri_1_gf + 31);

    auto g_xxyy_xxz_1 = pbuffer.data(idx_eri_1_gf + 32);

    auto g_xxyy_xyy_1 = pbuffer.data(idx_eri_1_gf + 33);

    auto g_xxyy_xyz_1 = pbuffer.data(idx_eri_1_gf + 34);

    auto g_xxyy_xzz_1 = pbuffer.data(idx_eri_1_gf + 35);

    auto g_xxyy_yyy_1 = pbuffer.data(idx_eri_1_gf + 36);

    auto g_xxyy_yyz_1 = pbuffer.data(idx_eri_1_gf + 37);

    auto g_xxyy_yzz_1 = pbuffer.data(idx_eri_1_gf + 38);

    auto g_xxyy_zzz_1 = pbuffer.data(idx_eri_1_gf + 39);

    auto g_xxzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 50);

    auto g_xxzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 51);

    auto g_xxzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 52);

    auto g_xxzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 53);

    auto g_xxzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 54);

    auto g_xxzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 55);

    auto g_xxzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 56);

    auto g_xxzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 57);

    auto g_xxzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 58);

    auto g_xxzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 59);

    auto g_xyyy_xxx_1 = pbuffer.data(idx_eri_1_gf + 60);

    auto g_xyyy_xxy_1 = pbuffer.data(idx_eri_1_gf + 61);

    auto g_xyyy_xyy_1 = pbuffer.data(idx_eri_1_gf + 63);

    auto g_xyyy_xyz_1 = pbuffer.data(idx_eri_1_gf + 64);

    auto g_xyyy_yyy_1 = pbuffer.data(idx_eri_1_gf + 66);

    auto g_xyyy_yyz_1 = pbuffer.data(idx_eri_1_gf + 67);

    auto g_xyyy_yzz_1 = pbuffer.data(idx_eri_1_gf + 68);

    auto g_xyyy_zzz_1 = pbuffer.data(idx_eri_1_gf + 69);

    auto g_xzzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 90);

    auto g_xzzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 92);

    auto g_xzzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 94);

    auto g_xzzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 95);

    auto g_xzzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 96);

    auto g_xzzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 97);

    auto g_xzzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 98);

    auto g_xzzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 99);

    auto g_yyyy_xxx_1 = pbuffer.data(idx_eri_1_gf + 100);

    auto g_yyyy_xxy_1 = pbuffer.data(idx_eri_1_gf + 101);

    auto g_yyyy_xxz_1 = pbuffer.data(idx_eri_1_gf + 102);

    auto g_yyyy_xyy_1 = pbuffer.data(idx_eri_1_gf + 103);

    auto g_yyyy_xyz_1 = pbuffer.data(idx_eri_1_gf + 104);

    auto g_yyyy_xzz_1 = pbuffer.data(idx_eri_1_gf + 105);

    auto g_yyyy_yyy_1 = pbuffer.data(idx_eri_1_gf + 106);

    auto g_yyyy_yyz_1 = pbuffer.data(idx_eri_1_gf + 107);

    auto g_yyyy_yzz_1 = pbuffer.data(idx_eri_1_gf + 108);

    auto g_yyyy_zzz_1 = pbuffer.data(idx_eri_1_gf + 109);

    auto g_yyyz_xxy_1 = pbuffer.data(idx_eri_1_gf + 111);

    auto g_yyyz_xxz_1 = pbuffer.data(idx_eri_1_gf + 112);

    auto g_yyyz_xyy_1 = pbuffer.data(idx_eri_1_gf + 113);

    auto g_yyyz_xyz_1 = pbuffer.data(idx_eri_1_gf + 114);

    auto g_yyyz_xzz_1 = pbuffer.data(idx_eri_1_gf + 115);

    auto g_yyyz_yyy_1 = pbuffer.data(idx_eri_1_gf + 116);

    auto g_yyyz_yyz_1 = pbuffer.data(idx_eri_1_gf + 117);

    auto g_yyyz_yzz_1 = pbuffer.data(idx_eri_1_gf + 118);

    auto g_yyyz_zzz_1 = pbuffer.data(idx_eri_1_gf + 119);

    auto g_yyzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 120);

    auto g_yyzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 121);

    auto g_yyzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 122);

    auto g_yyzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 123);

    auto g_yyzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 124);

    auto g_yyzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 125);

    auto g_yyzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 126);

    auto g_yyzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 127);

    auto g_yyzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 128);

    auto g_yyzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 129);

    auto g_yzzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 130);

    auto g_yzzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 131);

    auto g_yzzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 132);

    auto g_yzzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 133);

    auto g_yzzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 134);

    auto g_yzzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 135);

    auto g_yzzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 136);

    auto g_yzzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 137);

    auto g_yzzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 138);

    auto g_yzzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 139);

    auto g_zzzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 140);

    auto g_zzzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 141);

    auto g_zzzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 142);

    auto g_zzzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 143);

    auto g_zzzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 144);

    auto g_zzzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 145);

    auto g_zzzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 146);

    auto g_zzzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 147);

    auto g_zzzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 148);

    auto g_zzzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 149);

    // Set up 0-10 components of targeted buffer : HF

    auto g_xxxxx_xxx_0 = pbuffer.data(idx_eri_0_hf);

    auto g_xxxxx_xxy_0 = pbuffer.data(idx_eri_0_hf + 1);

    auto g_xxxxx_xxz_0 = pbuffer.data(idx_eri_0_hf + 2);

    auto g_xxxxx_xyy_0 = pbuffer.data(idx_eri_0_hf + 3);

    auto g_xxxxx_xyz_0 = pbuffer.data(idx_eri_0_hf + 4);

    auto g_xxxxx_xzz_0 = pbuffer.data(idx_eri_0_hf + 5);

    auto g_xxxxx_yyy_0 = pbuffer.data(idx_eri_0_hf + 6);

    auto g_xxxxx_yyz_0 = pbuffer.data(idx_eri_0_hf + 7);

    auto g_xxxxx_yzz_0 = pbuffer.data(idx_eri_0_hf + 8);

    auto g_xxxxx_zzz_0 = pbuffer.data(idx_eri_0_hf + 9);

    #pragma omp simd aligned(g_xxx_xxx_0, g_xxx_xxx_1, g_xxx_xxy_0, g_xxx_xxy_1, g_xxx_xxz_0, g_xxx_xxz_1, g_xxx_xyy_0, g_xxx_xyy_1, g_xxx_xyz_0, g_xxx_xyz_1, g_xxx_xzz_0, g_xxx_xzz_1, g_xxx_yyy_0, g_xxx_yyy_1, g_xxx_yyz_0, g_xxx_yyz_1, g_xxx_yzz_0, g_xxx_yzz_1, g_xxx_zzz_0, g_xxx_zzz_1, g_xxxx_xx_1, g_xxxx_xxx_1, g_xxxx_xxy_1, g_xxxx_xxz_1, g_xxxx_xy_1, g_xxxx_xyy_1, g_xxxx_xyz_1, g_xxxx_xz_1, g_xxxx_xzz_1, g_xxxx_yy_1, g_xxxx_yyy_1, g_xxxx_yyz_1, g_xxxx_yz_1, g_xxxx_yzz_1, g_xxxx_zz_1, g_xxxx_zzz_1, g_xxxxx_xxx_0, g_xxxxx_xxy_0, g_xxxxx_xxz_0, g_xxxxx_xyy_0, g_xxxxx_xyz_0, g_xxxxx_xzz_0, g_xxxxx_yyy_0, g_xxxxx_yyz_0, g_xxxxx_yzz_0, g_xxxxx_zzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxx_xxx_0[i] = 4.0 * g_xxx_xxx_0[i] * fbe_0 - 4.0 * g_xxx_xxx_1[i] * fz_be_0 + 3.0 * g_xxxx_xx_1[i] * fe_0 + g_xxxx_xxx_1[i] * pa_x[i];

        g_xxxxx_xxy_0[i] = 4.0 * g_xxx_xxy_0[i] * fbe_0 - 4.0 * g_xxx_xxy_1[i] * fz_be_0 + 2.0 * g_xxxx_xy_1[i] * fe_0 + g_xxxx_xxy_1[i] * pa_x[i];

        g_xxxxx_xxz_0[i] = 4.0 * g_xxx_xxz_0[i] * fbe_0 - 4.0 * g_xxx_xxz_1[i] * fz_be_0 + 2.0 * g_xxxx_xz_1[i] * fe_0 + g_xxxx_xxz_1[i] * pa_x[i];

        g_xxxxx_xyy_0[i] = 4.0 * g_xxx_xyy_0[i] * fbe_0 - 4.0 * g_xxx_xyy_1[i] * fz_be_0 + g_xxxx_yy_1[i] * fe_0 + g_xxxx_xyy_1[i] * pa_x[i];

        g_xxxxx_xyz_0[i] = 4.0 * g_xxx_xyz_0[i] * fbe_0 - 4.0 * g_xxx_xyz_1[i] * fz_be_0 + g_xxxx_yz_1[i] * fe_0 + g_xxxx_xyz_1[i] * pa_x[i];

        g_xxxxx_xzz_0[i] = 4.0 * g_xxx_xzz_0[i] * fbe_0 - 4.0 * g_xxx_xzz_1[i] * fz_be_0 + g_xxxx_zz_1[i] * fe_0 + g_xxxx_xzz_1[i] * pa_x[i];

        g_xxxxx_yyy_0[i] = 4.0 * g_xxx_yyy_0[i] * fbe_0 - 4.0 * g_xxx_yyy_1[i] * fz_be_0 + g_xxxx_yyy_1[i] * pa_x[i];

        g_xxxxx_yyz_0[i] = 4.0 * g_xxx_yyz_0[i] * fbe_0 - 4.0 * g_xxx_yyz_1[i] * fz_be_0 + g_xxxx_yyz_1[i] * pa_x[i];

        g_xxxxx_yzz_0[i] = 4.0 * g_xxx_yzz_0[i] * fbe_0 - 4.0 * g_xxx_yzz_1[i] * fz_be_0 + g_xxxx_yzz_1[i] * pa_x[i];

        g_xxxxx_zzz_0[i] = 4.0 * g_xxx_zzz_0[i] * fbe_0 - 4.0 * g_xxx_zzz_1[i] * fz_be_0 + g_xxxx_zzz_1[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : HF

    auto g_xxxxy_xxx_0 = pbuffer.data(idx_eri_0_hf + 10);

    auto g_xxxxy_xxy_0 = pbuffer.data(idx_eri_0_hf + 11);

    auto g_xxxxy_xxz_0 = pbuffer.data(idx_eri_0_hf + 12);

    auto g_xxxxy_xyy_0 = pbuffer.data(idx_eri_0_hf + 13);

    auto g_xxxxy_xyz_0 = pbuffer.data(idx_eri_0_hf + 14);

    auto g_xxxxy_xzz_0 = pbuffer.data(idx_eri_0_hf + 15);

    auto g_xxxxy_yyy_0 = pbuffer.data(idx_eri_0_hf + 16);

    auto g_xxxxy_yyz_0 = pbuffer.data(idx_eri_0_hf + 17);

    auto g_xxxxy_yzz_0 = pbuffer.data(idx_eri_0_hf + 18);

    auto g_xxxxy_zzz_0 = pbuffer.data(idx_eri_0_hf + 19);

    #pragma omp simd aligned(g_xxxx_xx_1, g_xxxx_xxx_1, g_xxxx_xxy_1, g_xxxx_xxz_1, g_xxxx_xy_1, g_xxxx_xyy_1, g_xxxx_xyz_1, g_xxxx_xz_1, g_xxxx_xzz_1, g_xxxx_yy_1, g_xxxx_yyy_1, g_xxxx_yyz_1, g_xxxx_yz_1, g_xxxx_yzz_1, g_xxxx_zz_1, g_xxxx_zzz_1, g_xxxxy_xxx_0, g_xxxxy_xxy_0, g_xxxxy_xxz_0, g_xxxxy_xyy_0, g_xxxxy_xyz_0, g_xxxxy_xzz_0, g_xxxxy_yyy_0, g_xxxxy_yyz_0, g_xxxxy_yzz_0, g_xxxxy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxy_xxx_0[i] = g_xxxx_xxx_1[i] * pa_y[i];

        g_xxxxy_xxy_0[i] = g_xxxx_xx_1[i] * fe_0 + g_xxxx_xxy_1[i] * pa_y[i];

        g_xxxxy_xxz_0[i] = g_xxxx_xxz_1[i] * pa_y[i];

        g_xxxxy_xyy_0[i] = 2.0 * g_xxxx_xy_1[i] * fe_0 + g_xxxx_xyy_1[i] * pa_y[i];

        g_xxxxy_xyz_0[i] = g_xxxx_xz_1[i] * fe_0 + g_xxxx_xyz_1[i] * pa_y[i];

        g_xxxxy_xzz_0[i] = g_xxxx_xzz_1[i] * pa_y[i];

        g_xxxxy_yyy_0[i] = 3.0 * g_xxxx_yy_1[i] * fe_0 + g_xxxx_yyy_1[i] * pa_y[i];

        g_xxxxy_yyz_0[i] = 2.0 * g_xxxx_yz_1[i] * fe_0 + g_xxxx_yyz_1[i] * pa_y[i];

        g_xxxxy_yzz_0[i] = g_xxxx_zz_1[i] * fe_0 + g_xxxx_yzz_1[i] * pa_y[i];

        g_xxxxy_zzz_0[i] = g_xxxx_zzz_1[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : HF

    auto g_xxxxz_xxx_0 = pbuffer.data(idx_eri_0_hf + 20);

    auto g_xxxxz_xxy_0 = pbuffer.data(idx_eri_0_hf + 21);

    auto g_xxxxz_xxz_0 = pbuffer.data(idx_eri_0_hf + 22);

    auto g_xxxxz_xyy_0 = pbuffer.data(idx_eri_0_hf + 23);

    auto g_xxxxz_xyz_0 = pbuffer.data(idx_eri_0_hf + 24);

    auto g_xxxxz_xzz_0 = pbuffer.data(idx_eri_0_hf + 25);

    auto g_xxxxz_yyy_0 = pbuffer.data(idx_eri_0_hf + 26);

    auto g_xxxxz_yyz_0 = pbuffer.data(idx_eri_0_hf + 27);

    auto g_xxxxz_yzz_0 = pbuffer.data(idx_eri_0_hf + 28);

    auto g_xxxxz_zzz_0 = pbuffer.data(idx_eri_0_hf + 29);

    #pragma omp simd aligned(g_xxxx_xx_1, g_xxxx_xxx_1, g_xxxx_xxy_1, g_xxxx_xxz_1, g_xxxx_xy_1, g_xxxx_xyy_1, g_xxxx_xyz_1, g_xxxx_xz_1, g_xxxx_xzz_1, g_xxxx_yy_1, g_xxxx_yyy_1, g_xxxx_yyz_1, g_xxxx_yz_1, g_xxxx_yzz_1, g_xxxx_zz_1, g_xxxx_zzz_1, g_xxxxz_xxx_0, g_xxxxz_xxy_0, g_xxxxz_xxz_0, g_xxxxz_xyy_0, g_xxxxz_xyz_0, g_xxxxz_xzz_0, g_xxxxz_yyy_0, g_xxxxz_yyz_0, g_xxxxz_yzz_0, g_xxxxz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxz_xxx_0[i] = g_xxxx_xxx_1[i] * pa_z[i];

        g_xxxxz_xxy_0[i] = g_xxxx_xxy_1[i] * pa_z[i];

        g_xxxxz_xxz_0[i] = g_xxxx_xx_1[i] * fe_0 + g_xxxx_xxz_1[i] * pa_z[i];

        g_xxxxz_xyy_0[i] = g_xxxx_xyy_1[i] * pa_z[i];

        g_xxxxz_xyz_0[i] = g_xxxx_xy_1[i] * fe_0 + g_xxxx_xyz_1[i] * pa_z[i];

        g_xxxxz_xzz_0[i] = 2.0 * g_xxxx_xz_1[i] * fe_0 + g_xxxx_xzz_1[i] * pa_z[i];

        g_xxxxz_yyy_0[i] = g_xxxx_yyy_1[i] * pa_z[i];

        g_xxxxz_yyz_0[i] = g_xxxx_yy_1[i] * fe_0 + g_xxxx_yyz_1[i] * pa_z[i];

        g_xxxxz_yzz_0[i] = 2.0 * g_xxxx_yz_1[i] * fe_0 + g_xxxx_yzz_1[i] * pa_z[i];

        g_xxxxz_zzz_0[i] = 3.0 * g_xxxx_zz_1[i] * fe_0 + g_xxxx_zzz_1[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : HF

    auto g_xxxyy_xxx_0 = pbuffer.data(idx_eri_0_hf + 30);

    auto g_xxxyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 31);

    auto g_xxxyy_xxz_0 = pbuffer.data(idx_eri_0_hf + 32);

    auto g_xxxyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 33);

    auto g_xxxyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 34);

    auto g_xxxyy_xzz_0 = pbuffer.data(idx_eri_0_hf + 35);

    auto g_xxxyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 36);

    auto g_xxxyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 37);

    auto g_xxxyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 38);

    auto g_xxxyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 39);

    #pragma omp simd aligned(g_xxx_xxx_0, g_xxx_xxx_1, g_xxx_xxz_0, g_xxx_xxz_1, g_xxx_xzz_0, g_xxx_xzz_1, g_xxxy_xxx_1, g_xxxy_xxz_1, g_xxxy_xzz_1, g_xxxyy_xxx_0, g_xxxyy_xxy_0, g_xxxyy_xxz_0, g_xxxyy_xyy_0, g_xxxyy_xyz_0, g_xxxyy_xzz_0, g_xxxyy_yyy_0, g_xxxyy_yyz_0, g_xxxyy_yzz_0, g_xxxyy_zzz_0, g_xxyy_xxy_1, g_xxyy_xy_1, g_xxyy_xyy_1, g_xxyy_xyz_1, g_xxyy_yy_1, g_xxyy_yyy_1, g_xxyy_yyz_1, g_xxyy_yz_1, g_xxyy_yzz_1, g_xxyy_zzz_1, g_xyy_xxy_0, g_xyy_xxy_1, g_xyy_xyy_0, g_xyy_xyy_1, g_xyy_xyz_0, g_xyy_xyz_1, g_xyy_yyy_0, g_xyy_yyy_1, g_xyy_yyz_0, g_xyy_yyz_1, g_xyy_yzz_0, g_xyy_yzz_1, g_xyy_zzz_0, g_xyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyy_xxx_0[i] = g_xxx_xxx_0[i] * fbe_0 - g_xxx_xxx_1[i] * fz_be_0 + g_xxxy_xxx_1[i] * pa_y[i];

        g_xxxyy_xxy_0[i] = 2.0 * g_xyy_xxy_0[i] * fbe_0 - 2.0 * g_xyy_xxy_1[i] * fz_be_0 + 2.0 * g_xxyy_xy_1[i] * fe_0 + g_xxyy_xxy_1[i] * pa_x[i];

        g_xxxyy_xxz_0[i] = g_xxx_xxz_0[i] * fbe_0 - g_xxx_xxz_1[i] * fz_be_0 + g_xxxy_xxz_1[i] * pa_y[i];

        g_xxxyy_xyy_0[i] = 2.0 * g_xyy_xyy_0[i] * fbe_0 - 2.0 * g_xyy_xyy_1[i] * fz_be_0 + g_xxyy_yy_1[i] * fe_0 + g_xxyy_xyy_1[i] * pa_x[i];

        g_xxxyy_xyz_0[i] = 2.0 * g_xyy_xyz_0[i] * fbe_0 - 2.0 * g_xyy_xyz_1[i] * fz_be_0 + g_xxyy_yz_1[i] * fe_0 + g_xxyy_xyz_1[i] * pa_x[i];

        g_xxxyy_xzz_0[i] = g_xxx_xzz_0[i] * fbe_0 - g_xxx_xzz_1[i] * fz_be_0 + g_xxxy_xzz_1[i] * pa_y[i];

        g_xxxyy_yyy_0[i] = 2.0 * g_xyy_yyy_0[i] * fbe_0 - 2.0 * g_xyy_yyy_1[i] * fz_be_0 + g_xxyy_yyy_1[i] * pa_x[i];

        g_xxxyy_yyz_0[i] = 2.0 * g_xyy_yyz_0[i] * fbe_0 - 2.0 * g_xyy_yyz_1[i] * fz_be_0 + g_xxyy_yyz_1[i] * pa_x[i];

        g_xxxyy_yzz_0[i] = 2.0 * g_xyy_yzz_0[i] * fbe_0 - 2.0 * g_xyy_yzz_1[i] * fz_be_0 + g_xxyy_yzz_1[i] * pa_x[i];

        g_xxxyy_zzz_0[i] = 2.0 * g_xyy_zzz_0[i] * fbe_0 - 2.0 * g_xyy_zzz_1[i] * fz_be_0 + g_xxyy_zzz_1[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : HF

    auto g_xxxyz_xxx_0 = pbuffer.data(idx_eri_0_hf + 40);

    auto g_xxxyz_xxy_0 = pbuffer.data(idx_eri_0_hf + 41);

    auto g_xxxyz_xxz_0 = pbuffer.data(idx_eri_0_hf + 42);

    auto g_xxxyz_xyy_0 = pbuffer.data(idx_eri_0_hf + 43);

    auto g_xxxyz_xyz_0 = pbuffer.data(idx_eri_0_hf + 44);

    auto g_xxxyz_xzz_0 = pbuffer.data(idx_eri_0_hf + 45);

    auto g_xxxyz_yyy_0 = pbuffer.data(idx_eri_0_hf + 46);

    auto g_xxxyz_yyz_0 = pbuffer.data(idx_eri_0_hf + 47);

    auto g_xxxyz_yzz_0 = pbuffer.data(idx_eri_0_hf + 48);

    auto g_xxxyz_zzz_0 = pbuffer.data(idx_eri_0_hf + 49);

    #pragma omp simd aligned(g_xxxy_xxy_1, g_xxxy_xyy_1, g_xxxy_yyy_1, g_xxxyz_xxx_0, g_xxxyz_xxy_0, g_xxxyz_xxz_0, g_xxxyz_xyy_0, g_xxxyz_xyz_0, g_xxxyz_xzz_0, g_xxxyz_yyy_0, g_xxxyz_yyz_0, g_xxxyz_yzz_0, g_xxxyz_zzz_0, g_xxxz_xxx_1, g_xxxz_xxz_1, g_xxxz_xyz_1, g_xxxz_xz_1, g_xxxz_xzz_1, g_xxxz_yyz_1, g_xxxz_yz_1, g_xxxz_yzz_1, g_xxxz_zz_1, g_xxxz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyz_xxx_0[i] = g_xxxz_xxx_1[i] * pa_y[i];

        g_xxxyz_xxy_0[i] = g_xxxy_xxy_1[i] * pa_z[i];

        g_xxxyz_xxz_0[i] = g_xxxz_xxz_1[i] * pa_y[i];

        g_xxxyz_xyy_0[i] = g_xxxy_xyy_1[i] * pa_z[i];

        g_xxxyz_xyz_0[i] = g_xxxz_xz_1[i] * fe_0 + g_xxxz_xyz_1[i] * pa_y[i];

        g_xxxyz_xzz_0[i] = g_xxxz_xzz_1[i] * pa_y[i];

        g_xxxyz_yyy_0[i] = g_xxxy_yyy_1[i] * pa_z[i];

        g_xxxyz_yyz_0[i] = 2.0 * g_xxxz_yz_1[i] * fe_0 + g_xxxz_yyz_1[i] * pa_y[i];

        g_xxxyz_yzz_0[i] = g_xxxz_zz_1[i] * fe_0 + g_xxxz_yzz_1[i] * pa_y[i];

        g_xxxyz_zzz_0[i] = g_xxxz_zzz_1[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : HF

    auto g_xxxzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 50);

    auto g_xxxzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 51);

    auto g_xxxzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 52);

    auto g_xxxzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 53);

    auto g_xxxzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 54);

    auto g_xxxzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 55);

    auto g_xxxzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 56);

    auto g_xxxzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 57);

    auto g_xxxzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 58);

    auto g_xxxzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 59);

    #pragma omp simd aligned(g_xxx_xxx_0, g_xxx_xxx_1, g_xxx_xxy_0, g_xxx_xxy_1, g_xxx_xyy_0, g_xxx_xyy_1, g_xxxz_xxx_1, g_xxxz_xxy_1, g_xxxz_xyy_1, g_xxxzz_xxx_0, g_xxxzz_xxy_0, g_xxxzz_xxz_0, g_xxxzz_xyy_0, g_xxxzz_xyz_0, g_xxxzz_xzz_0, g_xxxzz_yyy_0, g_xxxzz_yyz_0, g_xxxzz_yzz_0, g_xxxzz_zzz_0, g_xxzz_xxz_1, g_xxzz_xyz_1, g_xxzz_xz_1, g_xxzz_xzz_1, g_xxzz_yyy_1, g_xxzz_yyz_1, g_xxzz_yz_1, g_xxzz_yzz_1, g_xxzz_zz_1, g_xxzz_zzz_1, g_xzz_xxz_0, g_xzz_xxz_1, g_xzz_xyz_0, g_xzz_xyz_1, g_xzz_xzz_0, g_xzz_xzz_1, g_xzz_yyy_0, g_xzz_yyy_1, g_xzz_yyz_0, g_xzz_yyz_1, g_xzz_yzz_0, g_xzz_yzz_1, g_xzz_zzz_0, g_xzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzz_xxx_0[i] = g_xxx_xxx_0[i] * fbe_0 - g_xxx_xxx_1[i] * fz_be_0 + g_xxxz_xxx_1[i] * pa_z[i];

        g_xxxzz_xxy_0[i] = g_xxx_xxy_0[i] * fbe_0 - g_xxx_xxy_1[i] * fz_be_0 + g_xxxz_xxy_1[i] * pa_z[i];

        g_xxxzz_xxz_0[i] = 2.0 * g_xzz_xxz_0[i] * fbe_0 - 2.0 * g_xzz_xxz_1[i] * fz_be_0 + 2.0 * g_xxzz_xz_1[i] * fe_0 + g_xxzz_xxz_1[i] * pa_x[i];

        g_xxxzz_xyy_0[i] = g_xxx_xyy_0[i] * fbe_0 - g_xxx_xyy_1[i] * fz_be_0 + g_xxxz_xyy_1[i] * pa_z[i];

        g_xxxzz_xyz_0[i] = 2.0 * g_xzz_xyz_0[i] * fbe_0 - 2.0 * g_xzz_xyz_1[i] * fz_be_0 + g_xxzz_yz_1[i] * fe_0 + g_xxzz_xyz_1[i] * pa_x[i];

        g_xxxzz_xzz_0[i] = 2.0 * g_xzz_xzz_0[i] * fbe_0 - 2.0 * g_xzz_xzz_1[i] * fz_be_0 + g_xxzz_zz_1[i] * fe_0 + g_xxzz_xzz_1[i] * pa_x[i];

        g_xxxzz_yyy_0[i] = 2.0 * g_xzz_yyy_0[i] * fbe_0 - 2.0 * g_xzz_yyy_1[i] * fz_be_0 + g_xxzz_yyy_1[i] * pa_x[i];

        g_xxxzz_yyz_0[i] = 2.0 * g_xzz_yyz_0[i] * fbe_0 - 2.0 * g_xzz_yyz_1[i] * fz_be_0 + g_xxzz_yyz_1[i] * pa_x[i];

        g_xxxzz_yzz_0[i] = 2.0 * g_xzz_yzz_0[i] * fbe_0 - 2.0 * g_xzz_yzz_1[i] * fz_be_0 + g_xxzz_yzz_1[i] * pa_x[i];

        g_xxxzz_zzz_0[i] = 2.0 * g_xzz_zzz_0[i] * fbe_0 - 2.0 * g_xzz_zzz_1[i] * fz_be_0 + g_xxzz_zzz_1[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : HF

    auto g_xxyyy_xxx_0 = pbuffer.data(idx_eri_0_hf + 60);

    auto g_xxyyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 61);

    auto g_xxyyy_xxz_0 = pbuffer.data(idx_eri_0_hf + 62);

    auto g_xxyyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 63);

    auto g_xxyyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 64);

    auto g_xxyyy_xzz_0 = pbuffer.data(idx_eri_0_hf + 65);

    auto g_xxyyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 66);

    auto g_xxyyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 67);

    auto g_xxyyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 68);

    auto g_xxyyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 69);

    #pragma omp simd aligned(g_xxy_xxx_0, g_xxy_xxx_1, g_xxy_xxz_0, g_xxy_xxz_1, g_xxy_xzz_0, g_xxy_xzz_1, g_xxyy_xxx_1, g_xxyy_xxz_1, g_xxyy_xzz_1, g_xxyyy_xxx_0, g_xxyyy_xxy_0, g_xxyyy_xxz_0, g_xxyyy_xyy_0, g_xxyyy_xyz_0, g_xxyyy_xzz_0, g_xxyyy_yyy_0, g_xxyyy_yyz_0, g_xxyyy_yzz_0, g_xxyyy_zzz_0, g_xyyy_xxy_1, g_xyyy_xy_1, g_xyyy_xyy_1, g_xyyy_xyz_1, g_xyyy_yy_1, g_xyyy_yyy_1, g_xyyy_yyz_1, g_xyyy_yz_1, g_xyyy_yzz_1, g_xyyy_zzz_1, g_yyy_xxy_0, g_yyy_xxy_1, g_yyy_xyy_0, g_yyy_xyy_1, g_yyy_xyz_0, g_yyy_xyz_1, g_yyy_yyy_0, g_yyy_yyy_1, g_yyy_yyz_0, g_yyy_yyz_1, g_yyy_yzz_0, g_yyy_yzz_1, g_yyy_zzz_0, g_yyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyy_xxx_0[i] = 2.0 * g_xxy_xxx_0[i] * fbe_0 - 2.0 * g_xxy_xxx_1[i] * fz_be_0 + g_xxyy_xxx_1[i] * pa_y[i];

        g_xxyyy_xxy_0[i] = g_yyy_xxy_0[i] * fbe_0 - g_yyy_xxy_1[i] * fz_be_0 + 2.0 * g_xyyy_xy_1[i] * fe_0 + g_xyyy_xxy_1[i] * pa_x[i];

        g_xxyyy_xxz_0[i] = 2.0 * g_xxy_xxz_0[i] * fbe_0 - 2.0 * g_xxy_xxz_1[i] * fz_be_0 + g_xxyy_xxz_1[i] * pa_y[i];

        g_xxyyy_xyy_0[i] = g_yyy_xyy_0[i] * fbe_0 - g_yyy_xyy_1[i] * fz_be_0 + g_xyyy_yy_1[i] * fe_0 + g_xyyy_xyy_1[i] * pa_x[i];

        g_xxyyy_xyz_0[i] = g_yyy_xyz_0[i] * fbe_0 - g_yyy_xyz_1[i] * fz_be_0 + g_xyyy_yz_1[i] * fe_0 + g_xyyy_xyz_1[i] * pa_x[i];

        g_xxyyy_xzz_0[i] = 2.0 * g_xxy_xzz_0[i] * fbe_0 - 2.0 * g_xxy_xzz_1[i] * fz_be_0 + g_xxyy_xzz_1[i] * pa_y[i];

        g_xxyyy_yyy_0[i] = g_yyy_yyy_0[i] * fbe_0 - g_yyy_yyy_1[i] * fz_be_0 + g_xyyy_yyy_1[i] * pa_x[i];

        g_xxyyy_yyz_0[i] = g_yyy_yyz_0[i] * fbe_0 - g_yyy_yyz_1[i] * fz_be_0 + g_xyyy_yyz_1[i] * pa_x[i];

        g_xxyyy_yzz_0[i] = g_yyy_yzz_0[i] * fbe_0 - g_yyy_yzz_1[i] * fz_be_0 + g_xyyy_yzz_1[i] * pa_x[i];

        g_xxyyy_zzz_0[i] = g_yyy_zzz_0[i] * fbe_0 - g_yyy_zzz_1[i] * fz_be_0 + g_xyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : HF

    auto g_xxyyz_xxx_0 = pbuffer.data(idx_eri_0_hf + 70);

    auto g_xxyyz_xxy_0 = pbuffer.data(idx_eri_0_hf + 71);

    auto g_xxyyz_xxz_0 = pbuffer.data(idx_eri_0_hf + 72);

    auto g_xxyyz_xyy_0 = pbuffer.data(idx_eri_0_hf + 73);

    auto g_xxyyz_xyz_0 = pbuffer.data(idx_eri_0_hf + 74);

    auto g_xxyyz_xzz_0 = pbuffer.data(idx_eri_0_hf + 75);

    auto g_xxyyz_yyy_0 = pbuffer.data(idx_eri_0_hf + 76);

    auto g_xxyyz_yyz_0 = pbuffer.data(idx_eri_0_hf + 77);

    auto g_xxyyz_yzz_0 = pbuffer.data(idx_eri_0_hf + 78);

    auto g_xxyyz_zzz_0 = pbuffer.data(idx_eri_0_hf + 79);

    #pragma omp simd aligned(g_xxyy_xx_1, g_xxyy_xxx_1, g_xxyy_xxy_1, g_xxyy_xxz_1, g_xxyy_xy_1, g_xxyy_xyy_1, g_xxyy_xyz_1, g_xxyy_xz_1, g_xxyy_xzz_1, g_xxyy_yy_1, g_xxyy_yyy_1, g_xxyy_yyz_1, g_xxyy_yz_1, g_xxyy_yzz_1, g_xxyy_zz_1, g_xxyy_zzz_1, g_xxyyz_xxx_0, g_xxyyz_xxy_0, g_xxyyz_xxz_0, g_xxyyz_xyy_0, g_xxyyz_xyz_0, g_xxyyz_xzz_0, g_xxyyz_yyy_0, g_xxyyz_yyz_0, g_xxyyz_yzz_0, g_xxyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyz_xxx_0[i] = g_xxyy_xxx_1[i] * pa_z[i];

        g_xxyyz_xxy_0[i] = g_xxyy_xxy_1[i] * pa_z[i];

        g_xxyyz_xxz_0[i] = g_xxyy_xx_1[i] * fe_0 + g_xxyy_xxz_1[i] * pa_z[i];

        g_xxyyz_xyy_0[i] = g_xxyy_xyy_1[i] * pa_z[i];

        g_xxyyz_xyz_0[i] = g_xxyy_xy_1[i] * fe_0 + g_xxyy_xyz_1[i] * pa_z[i];

        g_xxyyz_xzz_0[i] = 2.0 * g_xxyy_xz_1[i] * fe_0 + g_xxyy_xzz_1[i] * pa_z[i];

        g_xxyyz_yyy_0[i] = g_xxyy_yyy_1[i] * pa_z[i];

        g_xxyyz_yyz_0[i] = g_xxyy_yy_1[i] * fe_0 + g_xxyy_yyz_1[i] * pa_z[i];

        g_xxyyz_yzz_0[i] = 2.0 * g_xxyy_yz_1[i] * fe_0 + g_xxyy_yzz_1[i] * pa_z[i];

        g_xxyyz_zzz_0[i] = 3.0 * g_xxyy_zz_1[i] * fe_0 + g_xxyy_zzz_1[i] * pa_z[i];
    }

    // Set up 80-90 components of targeted buffer : HF

    auto g_xxyzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 80);

    auto g_xxyzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 81);

    auto g_xxyzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 82);

    auto g_xxyzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 83);

    auto g_xxyzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 84);

    auto g_xxyzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 85);

    auto g_xxyzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 86);

    auto g_xxyzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 87);

    auto g_xxyzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 88);

    auto g_xxyzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 89);

    #pragma omp simd aligned(g_xxyzz_xxx_0, g_xxyzz_xxy_0, g_xxyzz_xxz_0, g_xxyzz_xyy_0, g_xxyzz_xyz_0, g_xxyzz_xzz_0, g_xxyzz_yyy_0, g_xxyzz_yyz_0, g_xxyzz_yzz_0, g_xxyzz_zzz_0, g_xxzz_xx_1, g_xxzz_xxx_1, g_xxzz_xxy_1, g_xxzz_xxz_1, g_xxzz_xy_1, g_xxzz_xyy_1, g_xxzz_xyz_1, g_xxzz_xz_1, g_xxzz_xzz_1, g_xxzz_yy_1, g_xxzz_yyy_1, g_xxzz_yyz_1, g_xxzz_yz_1, g_xxzz_yzz_1, g_xxzz_zz_1, g_xxzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzz_xxx_0[i] = g_xxzz_xxx_1[i] * pa_y[i];

        g_xxyzz_xxy_0[i] = g_xxzz_xx_1[i] * fe_0 + g_xxzz_xxy_1[i] * pa_y[i];

        g_xxyzz_xxz_0[i] = g_xxzz_xxz_1[i] * pa_y[i];

        g_xxyzz_xyy_0[i] = 2.0 * g_xxzz_xy_1[i] * fe_0 + g_xxzz_xyy_1[i] * pa_y[i];

        g_xxyzz_xyz_0[i] = g_xxzz_xz_1[i] * fe_0 + g_xxzz_xyz_1[i] * pa_y[i];

        g_xxyzz_xzz_0[i] = g_xxzz_xzz_1[i] * pa_y[i];

        g_xxyzz_yyy_0[i] = 3.0 * g_xxzz_yy_1[i] * fe_0 + g_xxzz_yyy_1[i] * pa_y[i];

        g_xxyzz_yyz_0[i] = 2.0 * g_xxzz_yz_1[i] * fe_0 + g_xxzz_yyz_1[i] * pa_y[i];

        g_xxyzz_yzz_0[i] = g_xxzz_zz_1[i] * fe_0 + g_xxzz_yzz_1[i] * pa_y[i];

        g_xxyzz_zzz_0[i] = g_xxzz_zzz_1[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : HF

    auto g_xxzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 90);

    auto g_xxzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 91);

    auto g_xxzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 92);

    auto g_xxzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 93);

    auto g_xxzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 94);

    auto g_xxzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 95);

    auto g_xxzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 96);

    auto g_xxzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 97);

    auto g_xxzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 98);

    auto g_xxzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 99);

    #pragma omp simd aligned(g_xxz_xxx_0, g_xxz_xxx_1, g_xxz_xxy_0, g_xxz_xxy_1, g_xxz_xyy_0, g_xxz_xyy_1, g_xxzz_xxx_1, g_xxzz_xxy_1, g_xxzz_xyy_1, g_xxzzz_xxx_0, g_xxzzz_xxy_0, g_xxzzz_xxz_0, g_xxzzz_xyy_0, g_xxzzz_xyz_0, g_xxzzz_xzz_0, g_xxzzz_yyy_0, g_xxzzz_yyz_0, g_xxzzz_yzz_0, g_xxzzz_zzz_0, g_xzzz_xxz_1, g_xzzz_xyz_1, g_xzzz_xz_1, g_xzzz_xzz_1, g_xzzz_yyy_1, g_xzzz_yyz_1, g_xzzz_yz_1, g_xzzz_yzz_1, g_xzzz_zz_1, g_xzzz_zzz_1, g_zzz_xxz_0, g_zzz_xxz_1, g_zzz_xyz_0, g_zzz_xyz_1, g_zzz_xzz_0, g_zzz_xzz_1, g_zzz_yyy_0, g_zzz_yyy_1, g_zzz_yyz_0, g_zzz_yyz_1, g_zzz_yzz_0, g_zzz_yzz_1, g_zzz_zzz_0, g_zzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzz_xxx_0[i] = 2.0 * g_xxz_xxx_0[i] * fbe_0 - 2.0 * g_xxz_xxx_1[i] * fz_be_0 + g_xxzz_xxx_1[i] * pa_z[i];

        g_xxzzz_xxy_0[i] = 2.0 * g_xxz_xxy_0[i] * fbe_0 - 2.0 * g_xxz_xxy_1[i] * fz_be_0 + g_xxzz_xxy_1[i] * pa_z[i];

        g_xxzzz_xxz_0[i] = g_zzz_xxz_0[i] * fbe_0 - g_zzz_xxz_1[i] * fz_be_0 + 2.0 * g_xzzz_xz_1[i] * fe_0 + g_xzzz_xxz_1[i] * pa_x[i];

        g_xxzzz_xyy_0[i] = 2.0 * g_xxz_xyy_0[i] * fbe_0 - 2.0 * g_xxz_xyy_1[i] * fz_be_0 + g_xxzz_xyy_1[i] * pa_z[i];

        g_xxzzz_xyz_0[i] = g_zzz_xyz_0[i] * fbe_0 - g_zzz_xyz_1[i] * fz_be_0 + g_xzzz_yz_1[i] * fe_0 + g_xzzz_xyz_1[i] * pa_x[i];

        g_xxzzz_xzz_0[i] = g_zzz_xzz_0[i] * fbe_0 - g_zzz_xzz_1[i] * fz_be_0 + g_xzzz_zz_1[i] * fe_0 + g_xzzz_xzz_1[i] * pa_x[i];

        g_xxzzz_yyy_0[i] = g_zzz_yyy_0[i] * fbe_0 - g_zzz_yyy_1[i] * fz_be_0 + g_xzzz_yyy_1[i] * pa_x[i];

        g_xxzzz_yyz_0[i] = g_zzz_yyz_0[i] * fbe_0 - g_zzz_yyz_1[i] * fz_be_0 + g_xzzz_yyz_1[i] * pa_x[i];

        g_xxzzz_yzz_0[i] = g_zzz_yzz_0[i] * fbe_0 - g_zzz_yzz_1[i] * fz_be_0 + g_xzzz_yzz_1[i] * pa_x[i];

        g_xxzzz_zzz_0[i] = g_zzz_zzz_0[i] * fbe_0 - g_zzz_zzz_1[i] * fz_be_0 + g_xzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : HF

    auto g_xyyyy_xxx_0 = pbuffer.data(idx_eri_0_hf + 100);

    auto g_xyyyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 101);

    auto g_xyyyy_xxz_0 = pbuffer.data(idx_eri_0_hf + 102);

    auto g_xyyyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 103);

    auto g_xyyyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 104);

    auto g_xyyyy_xzz_0 = pbuffer.data(idx_eri_0_hf + 105);

    auto g_xyyyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 106);

    auto g_xyyyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 107);

    auto g_xyyyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 108);

    auto g_xyyyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 109);

    #pragma omp simd aligned(g_xyyyy_xxx_0, g_xyyyy_xxy_0, g_xyyyy_xxz_0, g_xyyyy_xyy_0, g_xyyyy_xyz_0, g_xyyyy_xzz_0, g_xyyyy_yyy_0, g_xyyyy_yyz_0, g_xyyyy_yzz_0, g_xyyyy_zzz_0, g_yyyy_xx_1, g_yyyy_xxx_1, g_yyyy_xxy_1, g_yyyy_xxz_1, g_yyyy_xy_1, g_yyyy_xyy_1, g_yyyy_xyz_1, g_yyyy_xz_1, g_yyyy_xzz_1, g_yyyy_yy_1, g_yyyy_yyy_1, g_yyyy_yyz_1, g_yyyy_yz_1, g_yyyy_yzz_1, g_yyyy_zz_1, g_yyyy_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyy_xxx_0[i] = 3.0 * g_yyyy_xx_1[i] * fe_0 + g_yyyy_xxx_1[i] * pa_x[i];

        g_xyyyy_xxy_0[i] = 2.0 * g_yyyy_xy_1[i] * fe_0 + g_yyyy_xxy_1[i] * pa_x[i];

        g_xyyyy_xxz_0[i] = 2.0 * g_yyyy_xz_1[i] * fe_0 + g_yyyy_xxz_1[i] * pa_x[i];

        g_xyyyy_xyy_0[i] = g_yyyy_yy_1[i] * fe_0 + g_yyyy_xyy_1[i] * pa_x[i];

        g_xyyyy_xyz_0[i] = g_yyyy_yz_1[i] * fe_0 + g_yyyy_xyz_1[i] * pa_x[i];

        g_xyyyy_xzz_0[i] = g_yyyy_zz_1[i] * fe_0 + g_yyyy_xzz_1[i] * pa_x[i];

        g_xyyyy_yyy_0[i] = g_yyyy_yyy_1[i] * pa_x[i];

        g_xyyyy_yyz_0[i] = g_yyyy_yyz_1[i] * pa_x[i];

        g_xyyyy_yzz_0[i] = g_yyyy_yzz_1[i] * pa_x[i];

        g_xyyyy_zzz_0[i] = g_yyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 110-120 components of targeted buffer : HF

    auto g_xyyyz_xxx_0 = pbuffer.data(idx_eri_0_hf + 110);

    auto g_xyyyz_xxy_0 = pbuffer.data(idx_eri_0_hf + 111);

    auto g_xyyyz_xxz_0 = pbuffer.data(idx_eri_0_hf + 112);

    auto g_xyyyz_xyy_0 = pbuffer.data(idx_eri_0_hf + 113);

    auto g_xyyyz_xyz_0 = pbuffer.data(idx_eri_0_hf + 114);

    auto g_xyyyz_xzz_0 = pbuffer.data(idx_eri_0_hf + 115);

    auto g_xyyyz_yyy_0 = pbuffer.data(idx_eri_0_hf + 116);

    auto g_xyyyz_yyz_0 = pbuffer.data(idx_eri_0_hf + 117);

    auto g_xyyyz_yzz_0 = pbuffer.data(idx_eri_0_hf + 118);

    auto g_xyyyz_zzz_0 = pbuffer.data(idx_eri_0_hf + 119);

    #pragma omp simd aligned(g_xyyy_xxx_1, g_xyyy_xxy_1, g_xyyy_xyy_1, g_xyyyz_xxx_0, g_xyyyz_xxy_0, g_xyyyz_xxz_0, g_xyyyz_xyy_0, g_xyyyz_xyz_0, g_xyyyz_xzz_0, g_xyyyz_yyy_0, g_xyyyz_yyz_0, g_xyyyz_yzz_0, g_xyyyz_zzz_0, g_yyyz_xxz_1, g_yyyz_xyz_1, g_yyyz_xz_1, g_yyyz_xzz_1, g_yyyz_yyy_1, g_yyyz_yyz_1, g_yyyz_yz_1, g_yyyz_yzz_1, g_yyyz_zz_1, g_yyyz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyz_xxx_0[i] = g_xyyy_xxx_1[i] * pa_z[i];

        g_xyyyz_xxy_0[i] = g_xyyy_xxy_1[i] * pa_z[i];

        g_xyyyz_xxz_0[i] = 2.0 * g_yyyz_xz_1[i] * fe_0 + g_yyyz_xxz_1[i] * pa_x[i];

        g_xyyyz_xyy_0[i] = g_xyyy_xyy_1[i] * pa_z[i];

        g_xyyyz_xyz_0[i] = g_yyyz_yz_1[i] * fe_0 + g_yyyz_xyz_1[i] * pa_x[i];

        g_xyyyz_xzz_0[i] = g_yyyz_zz_1[i] * fe_0 + g_yyyz_xzz_1[i] * pa_x[i];

        g_xyyyz_yyy_0[i] = g_yyyz_yyy_1[i] * pa_x[i];

        g_xyyyz_yyz_0[i] = g_yyyz_yyz_1[i] * pa_x[i];

        g_xyyyz_yzz_0[i] = g_yyyz_yzz_1[i] * pa_x[i];

        g_xyyyz_zzz_0[i] = g_yyyz_zzz_1[i] * pa_x[i];
    }

    // Set up 120-130 components of targeted buffer : HF

    auto g_xyyzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 120);

    auto g_xyyzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 121);

    auto g_xyyzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 122);

    auto g_xyyzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 123);

    auto g_xyyzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 124);

    auto g_xyyzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 125);

    auto g_xyyzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 126);

    auto g_xyyzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 127);

    auto g_xyyzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 128);

    auto g_xyyzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 129);

    #pragma omp simd aligned(g_xyyzz_xxx_0, g_xyyzz_xxy_0, g_xyyzz_xxz_0, g_xyyzz_xyy_0, g_xyyzz_xyz_0, g_xyyzz_xzz_0, g_xyyzz_yyy_0, g_xyyzz_yyz_0, g_xyyzz_yzz_0, g_xyyzz_zzz_0, g_yyzz_xx_1, g_yyzz_xxx_1, g_yyzz_xxy_1, g_yyzz_xxz_1, g_yyzz_xy_1, g_yyzz_xyy_1, g_yyzz_xyz_1, g_yyzz_xz_1, g_yyzz_xzz_1, g_yyzz_yy_1, g_yyzz_yyy_1, g_yyzz_yyz_1, g_yyzz_yz_1, g_yyzz_yzz_1, g_yyzz_zz_1, g_yyzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzz_xxx_0[i] = 3.0 * g_yyzz_xx_1[i] * fe_0 + g_yyzz_xxx_1[i] * pa_x[i];

        g_xyyzz_xxy_0[i] = 2.0 * g_yyzz_xy_1[i] * fe_0 + g_yyzz_xxy_1[i] * pa_x[i];

        g_xyyzz_xxz_0[i] = 2.0 * g_yyzz_xz_1[i] * fe_0 + g_yyzz_xxz_1[i] * pa_x[i];

        g_xyyzz_xyy_0[i] = g_yyzz_yy_1[i] * fe_0 + g_yyzz_xyy_1[i] * pa_x[i];

        g_xyyzz_xyz_0[i] = g_yyzz_yz_1[i] * fe_0 + g_yyzz_xyz_1[i] * pa_x[i];

        g_xyyzz_xzz_0[i] = g_yyzz_zz_1[i] * fe_0 + g_yyzz_xzz_1[i] * pa_x[i];

        g_xyyzz_yyy_0[i] = g_yyzz_yyy_1[i] * pa_x[i];

        g_xyyzz_yyz_0[i] = g_yyzz_yyz_1[i] * pa_x[i];

        g_xyyzz_yzz_0[i] = g_yyzz_yzz_1[i] * pa_x[i];

        g_xyyzz_zzz_0[i] = g_yyzz_zzz_1[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : HF

    auto g_xyzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 130);

    auto g_xyzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 131);

    auto g_xyzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 132);

    auto g_xyzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 133);

    auto g_xyzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 134);

    auto g_xyzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 135);

    auto g_xyzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 136);

    auto g_xyzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 137);

    auto g_xyzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 138);

    auto g_xyzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 139);

    #pragma omp simd aligned(g_xyzzz_xxx_0, g_xyzzz_xxy_0, g_xyzzz_xxz_0, g_xyzzz_xyy_0, g_xyzzz_xyz_0, g_xyzzz_xzz_0, g_xyzzz_yyy_0, g_xyzzz_yyz_0, g_xyzzz_yzz_0, g_xyzzz_zzz_0, g_xzzz_xxx_1, g_xzzz_xxz_1, g_xzzz_xzz_1, g_yzzz_xxy_1, g_yzzz_xy_1, g_yzzz_xyy_1, g_yzzz_xyz_1, g_yzzz_yy_1, g_yzzz_yyy_1, g_yzzz_yyz_1, g_yzzz_yz_1, g_yzzz_yzz_1, g_yzzz_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzz_xxx_0[i] = g_xzzz_xxx_1[i] * pa_y[i];

        g_xyzzz_xxy_0[i] = 2.0 * g_yzzz_xy_1[i] * fe_0 + g_yzzz_xxy_1[i] * pa_x[i];

        g_xyzzz_xxz_0[i] = g_xzzz_xxz_1[i] * pa_y[i];

        g_xyzzz_xyy_0[i] = g_yzzz_yy_1[i] * fe_0 + g_yzzz_xyy_1[i] * pa_x[i];

        g_xyzzz_xyz_0[i] = g_yzzz_yz_1[i] * fe_0 + g_yzzz_xyz_1[i] * pa_x[i];

        g_xyzzz_xzz_0[i] = g_xzzz_xzz_1[i] * pa_y[i];

        g_xyzzz_yyy_0[i] = g_yzzz_yyy_1[i] * pa_x[i];

        g_xyzzz_yyz_0[i] = g_yzzz_yyz_1[i] * pa_x[i];

        g_xyzzz_yzz_0[i] = g_yzzz_yzz_1[i] * pa_x[i];

        g_xyzzz_zzz_0[i] = g_yzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 140-150 components of targeted buffer : HF

    auto g_xzzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 140);

    auto g_xzzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 141);

    auto g_xzzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 142);

    auto g_xzzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 143);

    auto g_xzzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 144);

    auto g_xzzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 145);

    auto g_xzzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 146);

    auto g_xzzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 147);

    auto g_xzzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 148);

    auto g_xzzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 149);

    #pragma omp simd aligned(g_xzzzz_xxx_0, g_xzzzz_xxy_0, g_xzzzz_xxz_0, g_xzzzz_xyy_0, g_xzzzz_xyz_0, g_xzzzz_xzz_0, g_xzzzz_yyy_0, g_xzzzz_yyz_0, g_xzzzz_yzz_0, g_xzzzz_zzz_0, g_zzzz_xx_1, g_zzzz_xxx_1, g_zzzz_xxy_1, g_zzzz_xxz_1, g_zzzz_xy_1, g_zzzz_xyy_1, g_zzzz_xyz_1, g_zzzz_xz_1, g_zzzz_xzz_1, g_zzzz_yy_1, g_zzzz_yyy_1, g_zzzz_yyz_1, g_zzzz_yz_1, g_zzzz_yzz_1, g_zzzz_zz_1, g_zzzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzz_xxx_0[i] = 3.0 * g_zzzz_xx_1[i] * fe_0 + g_zzzz_xxx_1[i] * pa_x[i];

        g_xzzzz_xxy_0[i] = 2.0 * g_zzzz_xy_1[i] * fe_0 + g_zzzz_xxy_1[i] * pa_x[i];

        g_xzzzz_xxz_0[i] = 2.0 * g_zzzz_xz_1[i] * fe_0 + g_zzzz_xxz_1[i] * pa_x[i];

        g_xzzzz_xyy_0[i] = g_zzzz_yy_1[i] * fe_0 + g_zzzz_xyy_1[i] * pa_x[i];

        g_xzzzz_xyz_0[i] = g_zzzz_yz_1[i] * fe_0 + g_zzzz_xyz_1[i] * pa_x[i];

        g_xzzzz_xzz_0[i] = g_zzzz_zz_1[i] * fe_0 + g_zzzz_xzz_1[i] * pa_x[i];

        g_xzzzz_yyy_0[i] = g_zzzz_yyy_1[i] * pa_x[i];

        g_xzzzz_yyz_0[i] = g_zzzz_yyz_1[i] * pa_x[i];

        g_xzzzz_yzz_0[i] = g_zzzz_yzz_1[i] * pa_x[i];

        g_xzzzz_zzz_0[i] = g_zzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : HF

    auto g_yyyyy_xxx_0 = pbuffer.data(idx_eri_0_hf + 150);

    auto g_yyyyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 151);

    auto g_yyyyy_xxz_0 = pbuffer.data(idx_eri_0_hf + 152);

    auto g_yyyyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 153);

    auto g_yyyyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 154);

    auto g_yyyyy_xzz_0 = pbuffer.data(idx_eri_0_hf + 155);

    auto g_yyyyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 156);

    auto g_yyyyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 157);

    auto g_yyyyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 158);

    auto g_yyyyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 159);

    #pragma omp simd aligned(g_yyy_xxx_0, g_yyy_xxx_1, g_yyy_xxy_0, g_yyy_xxy_1, g_yyy_xxz_0, g_yyy_xxz_1, g_yyy_xyy_0, g_yyy_xyy_1, g_yyy_xyz_0, g_yyy_xyz_1, g_yyy_xzz_0, g_yyy_xzz_1, g_yyy_yyy_0, g_yyy_yyy_1, g_yyy_yyz_0, g_yyy_yyz_1, g_yyy_yzz_0, g_yyy_yzz_1, g_yyy_zzz_0, g_yyy_zzz_1, g_yyyy_xx_1, g_yyyy_xxx_1, g_yyyy_xxy_1, g_yyyy_xxz_1, g_yyyy_xy_1, g_yyyy_xyy_1, g_yyyy_xyz_1, g_yyyy_xz_1, g_yyyy_xzz_1, g_yyyy_yy_1, g_yyyy_yyy_1, g_yyyy_yyz_1, g_yyyy_yz_1, g_yyyy_yzz_1, g_yyyy_zz_1, g_yyyy_zzz_1, g_yyyyy_xxx_0, g_yyyyy_xxy_0, g_yyyyy_xxz_0, g_yyyyy_xyy_0, g_yyyyy_xyz_0, g_yyyyy_xzz_0, g_yyyyy_yyy_0, g_yyyyy_yyz_0, g_yyyyy_yzz_0, g_yyyyy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyy_xxx_0[i] = 4.0 * g_yyy_xxx_0[i] * fbe_0 - 4.0 * g_yyy_xxx_1[i] * fz_be_0 + g_yyyy_xxx_1[i] * pa_y[i];

        g_yyyyy_xxy_0[i] = 4.0 * g_yyy_xxy_0[i] * fbe_0 - 4.0 * g_yyy_xxy_1[i] * fz_be_0 + g_yyyy_xx_1[i] * fe_0 + g_yyyy_xxy_1[i] * pa_y[i];

        g_yyyyy_xxz_0[i] = 4.0 * g_yyy_xxz_0[i] * fbe_0 - 4.0 * g_yyy_xxz_1[i] * fz_be_0 + g_yyyy_xxz_1[i] * pa_y[i];

        g_yyyyy_xyy_0[i] = 4.0 * g_yyy_xyy_0[i] * fbe_0 - 4.0 * g_yyy_xyy_1[i] * fz_be_0 + 2.0 * g_yyyy_xy_1[i] * fe_0 + g_yyyy_xyy_1[i] * pa_y[i];

        g_yyyyy_xyz_0[i] = 4.0 * g_yyy_xyz_0[i] * fbe_0 - 4.0 * g_yyy_xyz_1[i] * fz_be_0 + g_yyyy_xz_1[i] * fe_0 + g_yyyy_xyz_1[i] * pa_y[i];

        g_yyyyy_xzz_0[i] = 4.0 * g_yyy_xzz_0[i] * fbe_0 - 4.0 * g_yyy_xzz_1[i] * fz_be_0 + g_yyyy_xzz_1[i] * pa_y[i];

        g_yyyyy_yyy_0[i] = 4.0 * g_yyy_yyy_0[i] * fbe_0 - 4.0 * g_yyy_yyy_1[i] * fz_be_0 + 3.0 * g_yyyy_yy_1[i] * fe_0 + g_yyyy_yyy_1[i] * pa_y[i];

        g_yyyyy_yyz_0[i] = 4.0 * g_yyy_yyz_0[i] * fbe_0 - 4.0 * g_yyy_yyz_1[i] * fz_be_0 + 2.0 * g_yyyy_yz_1[i] * fe_0 + g_yyyy_yyz_1[i] * pa_y[i];

        g_yyyyy_yzz_0[i] = 4.0 * g_yyy_yzz_0[i] * fbe_0 - 4.0 * g_yyy_yzz_1[i] * fz_be_0 + g_yyyy_zz_1[i] * fe_0 + g_yyyy_yzz_1[i] * pa_y[i];

        g_yyyyy_zzz_0[i] = 4.0 * g_yyy_zzz_0[i] * fbe_0 - 4.0 * g_yyy_zzz_1[i] * fz_be_0 + g_yyyy_zzz_1[i] * pa_y[i];
    }

    // Set up 160-170 components of targeted buffer : HF

    auto g_yyyyz_xxx_0 = pbuffer.data(idx_eri_0_hf + 160);

    auto g_yyyyz_xxy_0 = pbuffer.data(idx_eri_0_hf + 161);

    auto g_yyyyz_xxz_0 = pbuffer.data(idx_eri_0_hf + 162);

    auto g_yyyyz_xyy_0 = pbuffer.data(idx_eri_0_hf + 163);

    auto g_yyyyz_xyz_0 = pbuffer.data(idx_eri_0_hf + 164);

    auto g_yyyyz_xzz_0 = pbuffer.data(idx_eri_0_hf + 165);

    auto g_yyyyz_yyy_0 = pbuffer.data(idx_eri_0_hf + 166);

    auto g_yyyyz_yyz_0 = pbuffer.data(idx_eri_0_hf + 167);

    auto g_yyyyz_yzz_0 = pbuffer.data(idx_eri_0_hf + 168);

    auto g_yyyyz_zzz_0 = pbuffer.data(idx_eri_0_hf + 169);

    #pragma omp simd aligned(g_yyyy_xx_1, g_yyyy_xxx_1, g_yyyy_xxy_1, g_yyyy_xxz_1, g_yyyy_xy_1, g_yyyy_xyy_1, g_yyyy_xyz_1, g_yyyy_xz_1, g_yyyy_xzz_1, g_yyyy_yy_1, g_yyyy_yyy_1, g_yyyy_yyz_1, g_yyyy_yz_1, g_yyyy_yzz_1, g_yyyy_zz_1, g_yyyy_zzz_1, g_yyyyz_xxx_0, g_yyyyz_xxy_0, g_yyyyz_xxz_0, g_yyyyz_xyy_0, g_yyyyz_xyz_0, g_yyyyz_xzz_0, g_yyyyz_yyy_0, g_yyyyz_yyz_0, g_yyyyz_yzz_0, g_yyyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyz_xxx_0[i] = g_yyyy_xxx_1[i] * pa_z[i];

        g_yyyyz_xxy_0[i] = g_yyyy_xxy_1[i] * pa_z[i];

        g_yyyyz_xxz_0[i] = g_yyyy_xx_1[i] * fe_0 + g_yyyy_xxz_1[i] * pa_z[i];

        g_yyyyz_xyy_0[i] = g_yyyy_xyy_1[i] * pa_z[i];

        g_yyyyz_xyz_0[i] = g_yyyy_xy_1[i] * fe_0 + g_yyyy_xyz_1[i] * pa_z[i];

        g_yyyyz_xzz_0[i] = 2.0 * g_yyyy_xz_1[i] * fe_0 + g_yyyy_xzz_1[i] * pa_z[i];

        g_yyyyz_yyy_0[i] = g_yyyy_yyy_1[i] * pa_z[i];

        g_yyyyz_yyz_0[i] = g_yyyy_yy_1[i] * fe_0 + g_yyyy_yyz_1[i] * pa_z[i];

        g_yyyyz_yzz_0[i] = 2.0 * g_yyyy_yz_1[i] * fe_0 + g_yyyy_yzz_1[i] * pa_z[i];

        g_yyyyz_zzz_0[i] = 3.0 * g_yyyy_zz_1[i] * fe_0 + g_yyyy_zzz_1[i] * pa_z[i];
    }

    // Set up 170-180 components of targeted buffer : HF

    auto g_yyyzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 170);

    auto g_yyyzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 171);

    auto g_yyyzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 172);

    auto g_yyyzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 173);

    auto g_yyyzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 174);

    auto g_yyyzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 175);

    auto g_yyyzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 176);

    auto g_yyyzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 177);

    auto g_yyyzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 178);

    auto g_yyyzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 179);

    #pragma omp simd aligned(g_yyy_xxy_0, g_yyy_xxy_1, g_yyy_xyy_0, g_yyy_xyy_1, g_yyy_yyy_0, g_yyy_yyy_1, g_yyyz_xxy_1, g_yyyz_xyy_1, g_yyyz_yyy_1, g_yyyzz_xxx_0, g_yyyzz_xxy_0, g_yyyzz_xxz_0, g_yyyzz_xyy_0, g_yyyzz_xyz_0, g_yyyzz_xzz_0, g_yyyzz_yyy_0, g_yyyzz_yyz_0, g_yyyzz_yzz_0, g_yyyzz_zzz_0, g_yyzz_xxx_1, g_yyzz_xxz_1, g_yyzz_xyz_1, g_yyzz_xz_1, g_yyzz_xzz_1, g_yyzz_yyz_1, g_yyzz_yz_1, g_yyzz_yzz_1, g_yyzz_zz_1, g_yyzz_zzz_1, g_yzz_xxx_0, g_yzz_xxx_1, g_yzz_xxz_0, g_yzz_xxz_1, g_yzz_xyz_0, g_yzz_xyz_1, g_yzz_xzz_0, g_yzz_xzz_1, g_yzz_yyz_0, g_yzz_yyz_1, g_yzz_yzz_0, g_yzz_yzz_1, g_yzz_zzz_0, g_yzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzz_xxx_0[i] = 2.0 * g_yzz_xxx_0[i] * fbe_0 - 2.0 * g_yzz_xxx_1[i] * fz_be_0 + g_yyzz_xxx_1[i] * pa_y[i];

        g_yyyzz_xxy_0[i] = g_yyy_xxy_0[i] * fbe_0 - g_yyy_xxy_1[i] * fz_be_0 + g_yyyz_xxy_1[i] * pa_z[i];

        g_yyyzz_xxz_0[i] = 2.0 * g_yzz_xxz_0[i] * fbe_0 - 2.0 * g_yzz_xxz_1[i] * fz_be_0 + g_yyzz_xxz_1[i] * pa_y[i];

        g_yyyzz_xyy_0[i] = g_yyy_xyy_0[i] * fbe_0 - g_yyy_xyy_1[i] * fz_be_0 + g_yyyz_xyy_1[i] * pa_z[i];

        g_yyyzz_xyz_0[i] = 2.0 * g_yzz_xyz_0[i] * fbe_0 - 2.0 * g_yzz_xyz_1[i] * fz_be_0 + g_yyzz_xz_1[i] * fe_0 + g_yyzz_xyz_1[i] * pa_y[i];

        g_yyyzz_xzz_0[i] = 2.0 * g_yzz_xzz_0[i] * fbe_0 - 2.0 * g_yzz_xzz_1[i] * fz_be_0 + g_yyzz_xzz_1[i] * pa_y[i];

        g_yyyzz_yyy_0[i] = g_yyy_yyy_0[i] * fbe_0 - g_yyy_yyy_1[i] * fz_be_0 + g_yyyz_yyy_1[i] * pa_z[i];

        g_yyyzz_yyz_0[i] = 2.0 * g_yzz_yyz_0[i] * fbe_0 - 2.0 * g_yzz_yyz_1[i] * fz_be_0 + 2.0 * g_yyzz_yz_1[i] * fe_0 + g_yyzz_yyz_1[i] * pa_y[i];

        g_yyyzz_yzz_0[i] = 2.0 * g_yzz_yzz_0[i] * fbe_0 - 2.0 * g_yzz_yzz_1[i] * fz_be_0 + g_yyzz_zz_1[i] * fe_0 + g_yyzz_yzz_1[i] * pa_y[i];

        g_yyyzz_zzz_0[i] = 2.0 * g_yzz_zzz_0[i] * fbe_0 - 2.0 * g_yzz_zzz_1[i] * fz_be_0 + g_yyzz_zzz_1[i] * pa_y[i];
    }

    // Set up 180-190 components of targeted buffer : HF

    auto g_yyzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 180);

    auto g_yyzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 181);

    auto g_yyzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 182);

    auto g_yyzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 183);

    auto g_yyzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 184);

    auto g_yyzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 185);

    auto g_yyzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 186);

    auto g_yyzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 187);

    auto g_yyzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 188);

    auto g_yyzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 189);

    #pragma omp simd aligned(g_yyz_xxy_0, g_yyz_xxy_1, g_yyz_xyy_0, g_yyz_xyy_1, g_yyz_yyy_0, g_yyz_yyy_1, g_yyzz_xxy_1, g_yyzz_xyy_1, g_yyzz_yyy_1, g_yyzzz_xxx_0, g_yyzzz_xxy_0, g_yyzzz_xxz_0, g_yyzzz_xyy_0, g_yyzzz_xyz_0, g_yyzzz_xzz_0, g_yyzzz_yyy_0, g_yyzzz_yyz_0, g_yyzzz_yzz_0, g_yyzzz_zzz_0, g_yzzz_xxx_1, g_yzzz_xxz_1, g_yzzz_xyz_1, g_yzzz_xz_1, g_yzzz_xzz_1, g_yzzz_yyz_1, g_yzzz_yz_1, g_yzzz_yzz_1, g_yzzz_zz_1, g_yzzz_zzz_1, g_zzz_xxx_0, g_zzz_xxx_1, g_zzz_xxz_0, g_zzz_xxz_1, g_zzz_xyz_0, g_zzz_xyz_1, g_zzz_xzz_0, g_zzz_xzz_1, g_zzz_yyz_0, g_zzz_yyz_1, g_zzz_yzz_0, g_zzz_yzz_1, g_zzz_zzz_0, g_zzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzz_xxx_0[i] = g_zzz_xxx_0[i] * fbe_0 - g_zzz_xxx_1[i] * fz_be_0 + g_yzzz_xxx_1[i] * pa_y[i];

        g_yyzzz_xxy_0[i] = 2.0 * g_yyz_xxy_0[i] * fbe_0 - 2.0 * g_yyz_xxy_1[i] * fz_be_0 + g_yyzz_xxy_1[i] * pa_z[i];

        g_yyzzz_xxz_0[i] = g_zzz_xxz_0[i] * fbe_0 - g_zzz_xxz_1[i] * fz_be_0 + g_yzzz_xxz_1[i] * pa_y[i];

        g_yyzzz_xyy_0[i] = 2.0 * g_yyz_xyy_0[i] * fbe_0 - 2.0 * g_yyz_xyy_1[i] * fz_be_0 + g_yyzz_xyy_1[i] * pa_z[i];

        g_yyzzz_xyz_0[i] = g_zzz_xyz_0[i] * fbe_0 - g_zzz_xyz_1[i] * fz_be_0 + g_yzzz_xz_1[i] * fe_0 + g_yzzz_xyz_1[i] * pa_y[i];

        g_yyzzz_xzz_0[i] = g_zzz_xzz_0[i] * fbe_0 - g_zzz_xzz_1[i] * fz_be_0 + g_yzzz_xzz_1[i] * pa_y[i];

        g_yyzzz_yyy_0[i] = 2.0 * g_yyz_yyy_0[i] * fbe_0 - 2.0 * g_yyz_yyy_1[i] * fz_be_0 + g_yyzz_yyy_1[i] * pa_z[i];

        g_yyzzz_yyz_0[i] = g_zzz_yyz_0[i] * fbe_0 - g_zzz_yyz_1[i] * fz_be_0 + 2.0 * g_yzzz_yz_1[i] * fe_0 + g_yzzz_yyz_1[i] * pa_y[i];

        g_yyzzz_yzz_0[i] = g_zzz_yzz_0[i] * fbe_0 - g_zzz_yzz_1[i] * fz_be_0 + g_yzzz_zz_1[i] * fe_0 + g_yzzz_yzz_1[i] * pa_y[i];

        g_yyzzz_zzz_0[i] = g_zzz_zzz_0[i] * fbe_0 - g_zzz_zzz_1[i] * fz_be_0 + g_yzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 190-200 components of targeted buffer : HF

    auto g_yzzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 190);

    auto g_yzzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 191);

    auto g_yzzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 192);

    auto g_yzzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 193);

    auto g_yzzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 194);

    auto g_yzzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 195);

    auto g_yzzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 196);

    auto g_yzzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 197);

    auto g_yzzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 198);

    auto g_yzzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 199);

    #pragma omp simd aligned(g_yzzzz_xxx_0, g_yzzzz_xxy_0, g_yzzzz_xxz_0, g_yzzzz_xyy_0, g_yzzzz_xyz_0, g_yzzzz_xzz_0, g_yzzzz_yyy_0, g_yzzzz_yyz_0, g_yzzzz_yzz_0, g_yzzzz_zzz_0, g_zzzz_xx_1, g_zzzz_xxx_1, g_zzzz_xxy_1, g_zzzz_xxz_1, g_zzzz_xy_1, g_zzzz_xyy_1, g_zzzz_xyz_1, g_zzzz_xz_1, g_zzzz_xzz_1, g_zzzz_yy_1, g_zzzz_yyy_1, g_zzzz_yyz_1, g_zzzz_yz_1, g_zzzz_yzz_1, g_zzzz_zz_1, g_zzzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzz_xxx_0[i] = g_zzzz_xxx_1[i] * pa_y[i];

        g_yzzzz_xxy_0[i] = g_zzzz_xx_1[i] * fe_0 + g_zzzz_xxy_1[i] * pa_y[i];

        g_yzzzz_xxz_0[i] = g_zzzz_xxz_1[i] * pa_y[i];

        g_yzzzz_xyy_0[i] = 2.0 * g_zzzz_xy_1[i] * fe_0 + g_zzzz_xyy_1[i] * pa_y[i];

        g_yzzzz_xyz_0[i] = g_zzzz_xz_1[i] * fe_0 + g_zzzz_xyz_1[i] * pa_y[i];

        g_yzzzz_xzz_0[i] = g_zzzz_xzz_1[i] * pa_y[i];

        g_yzzzz_yyy_0[i] = 3.0 * g_zzzz_yy_1[i] * fe_0 + g_zzzz_yyy_1[i] * pa_y[i];

        g_yzzzz_yyz_0[i] = 2.0 * g_zzzz_yz_1[i] * fe_0 + g_zzzz_yyz_1[i] * pa_y[i];

        g_yzzzz_yzz_0[i] = g_zzzz_zz_1[i] * fe_0 + g_zzzz_yzz_1[i] * pa_y[i];

        g_yzzzz_zzz_0[i] = g_zzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 200-210 components of targeted buffer : HF

    auto g_zzzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 200);

    auto g_zzzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 201);

    auto g_zzzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 202);

    auto g_zzzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 203);

    auto g_zzzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 204);

    auto g_zzzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 205);

    auto g_zzzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 206);

    auto g_zzzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 207);

    auto g_zzzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 208);

    auto g_zzzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 209);

    #pragma omp simd aligned(g_zzz_xxx_0, g_zzz_xxx_1, g_zzz_xxy_0, g_zzz_xxy_1, g_zzz_xxz_0, g_zzz_xxz_1, g_zzz_xyy_0, g_zzz_xyy_1, g_zzz_xyz_0, g_zzz_xyz_1, g_zzz_xzz_0, g_zzz_xzz_1, g_zzz_yyy_0, g_zzz_yyy_1, g_zzz_yyz_0, g_zzz_yyz_1, g_zzz_yzz_0, g_zzz_yzz_1, g_zzz_zzz_0, g_zzz_zzz_1, g_zzzz_xx_1, g_zzzz_xxx_1, g_zzzz_xxy_1, g_zzzz_xxz_1, g_zzzz_xy_1, g_zzzz_xyy_1, g_zzzz_xyz_1, g_zzzz_xz_1, g_zzzz_xzz_1, g_zzzz_yy_1, g_zzzz_yyy_1, g_zzzz_yyz_1, g_zzzz_yz_1, g_zzzz_yzz_1, g_zzzz_zz_1, g_zzzz_zzz_1, g_zzzzz_xxx_0, g_zzzzz_xxy_0, g_zzzzz_xxz_0, g_zzzzz_xyy_0, g_zzzzz_xyz_0, g_zzzzz_xzz_0, g_zzzzz_yyy_0, g_zzzzz_yyz_0, g_zzzzz_yzz_0, g_zzzzz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzz_xxx_0[i] = 4.0 * g_zzz_xxx_0[i] * fbe_0 - 4.0 * g_zzz_xxx_1[i] * fz_be_0 + g_zzzz_xxx_1[i] * pa_z[i];

        g_zzzzz_xxy_0[i] = 4.0 * g_zzz_xxy_0[i] * fbe_0 - 4.0 * g_zzz_xxy_1[i] * fz_be_0 + g_zzzz_xxy_1[i] * pa_z[i];

        g_zzzzz_xxz_0[i] = 4.0 * g_zzz_xxz_0[i] * fbe_0 - 4.0 * g_zzz_xxz_1[i] * fz_be_0 + g_zzzz_xx_1[i] * fe_0 + g_zzzz_xxz_1[i] * pa_z[i];

        g_zzzzz_xyy_0[i] = 4.0 * g_zzz_xyy_0[i] * fbe_0 - 4.0 * g_zzz_xyy_1[i] * fz_be_0 + g_zzzz_xyy_1[i] * pa_z[i];

        g_zzzzz_xyz_0[i] = 4.0 * g_zzz_xyz_0[i] * fbe_0 - 4.0 * g_zzz_xyz_1[i] * fz_be_0 + g_zzzz_xy_1[i] * fe_0 + g_zzzz_xyz_1[i] * pa_z[i];

        g_zzzzz_xzz_0[i] = 4.0 * g_zzz_xzz_0[i] * fbe_0 - 4.0 * g_zzz_xzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xz_1[i] * fe_0 + g_zzzz_xzz_1[i] * pa_z[i];

        g_zzzzz_yyy_0[i] = 4.0 * g_zzz_yyy_0[i] * fbe_0 - 4.0 * g_zzz_yyy_1[i] * fz_be_0 + g_zzzz_yyy_1[i] * pa_z[i];

        g_zzzzz_yyz_0[i] = 4.0 * g_zzz_yyz_0[i] * fbe_0 - 4.0 * g_zzz_yyz_1[i] * fz_be_0 + g_zzzz_yy_1[i] * fe_0 + g_zzzz_yyz_1[i] * pa_z[i];

        g_zzzzz_yzz_0[i] = 4.0 * g_zzz_yzz_0[i] * fbe_0 - 4.0 * g_zzz_yzz_1[i] * fz_be_0 + 2.0 * g_zzzz_yz_1[i] * fe_0 + g_zzzz_yzz_1[i] * pa_z[i];

        g_zzzzz_zzz_0[i] = 4.0 * g_zzz_zzz_0[i] * fbe_0 - 4.0 * g_zzz_zzz_1[i] * fz_be_0 + 3.0 * g_zzzz_zz_1[i] * fe_0 + g_zzzz_zzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

