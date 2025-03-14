#include "TwoCenterElectronRepulsionPrimRecIF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_if(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_if,
                                const size_t idx_eri_0_gf,
                                const size_t idx_eri_1_gf,
                                const size_t idx_eri_1_hd,
                                const size_t idx_eri_1_hf,
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

    // Set up components of auxiliary buffer : GF

    auto g_xxxx_xxx_0 = pbuffer.data(idx_eri_0_gf);

    auto g_xxxx_xxy_0 = pbuffer.data(idx_eri_0_gf + 1);

    auto g_xxxx_xxz_0 = pbuffer.data(idx_eri_0_gf + 2);

    auto g_xxxx_xyy_0 = pbuffer.data(idx_eri_0_gf + 3);

    auto g_xxxx_xyz_0 = pbuffer.data(idx_eri_0_gf + 4);

    auto g_xxxx_xzz_0 = pbuffer.data(idx_eri_0_gf + 5);

    auto g_xxxx_yyy_0 = pbuffer.data(idx_eri_0_gf + 6);

    auto g_xxxx_yyz_0 = pbuffer.data(idx_eri_0_gf + 7);

    auto g_xxxx_yzz_0 = pbuffer.data(idx_eri_0_gf + 8);

    auto g_xxxx_zzz_0 = pbuffer.data(idx_eri_0_gf + 9);

    auto g_xxxy_xxx_0 = pbuffer.data(idx_eri_0_gf + 10);

    auto g_xxxy_xxz_0 = pbuffer.data(idx_eri_0_gf + 12);

    auto g_xxxy_xzz_0 = pbuffer.data(idx_eri_0_gf + 15);

    auto g_xxxz_xxx_0 = pbuffer.data(idx_eri_0_gf + 20);

    auto g_xxxz_xxy_0 = pbuffer.data(idx_eri_0_gf + 21);

    auto g_xxxz_xyy_0 = pbuffer.data(idx_eri_0_gf + 23);

    auto g_xxyy_xxx_0 = pbuffer.data(idx_eri_0_gf + 30);

    auto g_xxyy_xxy_0 = pbuffer.data(idx_eri_0_gf + 31);

    auto g_xxyy_xxz_0 = pbuffer.data(idx_eri_0_gf + 32);

    auto g_xxyy_xyy_0 = pbuffer.data(idx_eri_0_gf + 33);

    auto g_xxyy_xyz_0 = pbuffer.data(idx_eri_0_gf + 34);

    auto g_xxyy_xzz_0 = pbuffer.data(idx_eri_0_gf + 35);

    auto g_xxyy_yyy_0 = pbuffer.data(idx_eri_0_gf + 36);

    auto g_xxyy_yyz_0 = pbuffer.data(idx_eri_0_gf + 37);

    auto g_xxyy_yzz_0 = pbuffer.data(idx_eri_0_gf + 38);

    auto g_xxyy_zzz_0 = pbuffer.data(idx_eri_0_gf + 39);

    auto g_xxzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 50);

    auto g_xxzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 51);

    auto g_xxzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 52);

    auto g_xxzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 53);

    auto g_xxzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 54);

    auto g_xxzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 55);

    auto g_xxzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 56);

    auto g_xxzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 57);

    auto g_xxzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 58);

    auto g_xxzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 59);

    auto g_xyyy_xxy_0 = pbuffer.data(idx_eri_0_gf + 61);

    auto g_xyyy_xyy_0 = pbuffer.data(idx_eri_0_gf + 63);

    auto g_xyyy_xyz_0 = pbuffer.data(idx_eri_0_gf + 64);

    auto g_xyyy_yyy_0 = pbuffer.data(idx_eri_0_gf + 66);

    auto g_xyyy_yyz_0 = pbuffer.data(idx_eri_0_gf + 67);

    auto g_xyyy_yzz_0 = pbuffer.data(idx_eri_0_gf + 68);

    auto g_xyyy_zzz_0 = pbuffer.data(idx_eri_0_gf + 69);

    auto g_xzzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 92);

    auto g_xzzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 94);

    auto g_xzzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 95);

    auto g_xzzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 96);

    auto g_xzzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 97);

    auto g_xzzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 98);

    auto g_xzzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 99);

    auto g_yyyy_xxx_0 = pbuffer.data(idx_eri_0_gf + 100);

    auto g_yyyy_xxy_0 = pbuffer.data(idx_eri_0_gf + 101);

    auto g_yyyy_xxz_0 = pbuffer.data(idx_eri_0_gf + 102);

    auto g_yyyy_xyy_0 = pbuffer.data(idx_eri_0_gf + 103);

    auto g_yyyy_xyz_0 = pbuffer.data(idx_eri_0_gf + 104);

    auto g_yyyy_xzz_0 = pbuffer.data(idx_eri_0_gf + 105);

    auto g_yyyy_yyy_0 = pbuffer.data(idx_eri_0_gf + 106);

    auto g_yyyy_yyz_0 = pbuffer.data(idx_eri_0_gf + 107);

    auto g_yyyy_yzz_0 = pbuffer.data(idx_eri_0_gf + 108);

    auto g_yyyy_zzz_0 = pbuffer.data(idx_eri_0_gf + 109);

    auto g_yyyz_xxy_0 = pbuffer.data(idx_eri_0_gf + 111);

    auto g_yyyz_xyy_0 = pbuffer.data(idx_eri_0_gf + 113);

    auto g_yyyz_yyy_0 = pbuffer.data(idx_eri_0_gf + 116);

    auto g_yyzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 120);

    auto g_yyzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 121);

    auto g_yyzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 122);

    auto g_yyzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 123);

    auto g_yyzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 124);

    auto g_yyzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 125);

    auto g_yyzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 126);

    auto g_yyzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 127);

    auto g_yyzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 128);

    auto g_yyzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 129);

    auto g_yzzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 130);

    auto g_yzzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 132);

    auto g_yzzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 134);

    auto g_yzzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 135);

    auto g_yzzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 137);

    auto g_yzzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 138);

    auto g_yzzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 139);

    auto g_zzzz_xxx_0 = pbuffer.data(idx_eri_0_gf + 140);

    auto g_zzzz_xxy_0 = pbuffer.data(idx_eri_0_gf + 141);

    auto g_zzzz_xxz_0 = pbuffer.data(idx_eri_0_gf + 142);

    auto g_zzzz_xyy_0 = pbuffer.data(idx_eri_0_gf + 143);

    auto g_zzzz_xyz_0 = pbuffer.data(idx_eri_0_gf + 144);

    auto g_zzzz_xzz_0 = pbuffer.data(idx_eri_0_gf + 145);

    auto g_zzzz_yyy_0 = pbuffer.data(idx_eri_0_gf + 146);

    auto g_zzzz_yyz_0 = pbuffer.data(idx_eri_0_gf + 147);

    auto g_zzzz_yzz_0 = pbuffer.data(idx_eri_0_gf + 148);

    auto g_zzzz_zzz_0 = pbuffer.data(idx_eri_0_gf + 149);

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

    auto g_xxxy_xxz_1 = pbuffer.data(idx_eri_1_gf + 12);

    auto g_xxxy_xzz_1 = pbuffer.data(idx_eri_1_gf + 15);

    auto g_xxxz_xxx_1 = pbuffer.data(idx_eri_1_gf + 20);

    auto g_xxxz_xxy_1 = pbuffer.data(idx_eri_1_gf + 21);

    auto g_xxxz_xyy_1 = pbuffer.data(idx_eri_1_gf + 23);

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

    auto g_xyyy_xxy_1 = pbuffer.data(idx_eri_1_gf + 61);

    auto g_xyyy_xyy_1 = pbuffer.data(idx_eri_1_gf + 63);

    auto g_xyyy_xyz_1 = pbuffer.data(idx_eri_1_gf + 64);

    auto g_xyyy_yyy_1 = pbuffer.data(idx_eri_1_gf + 66);

    auto g_xyyy_yyz_1 = pbuffer.data(idx_eri_1_gf + 67);

    auto g_xyyy_yzz_1 = pbuffer.data(idx_eri_1_gf + 68);

    auto g_xyyy_zzz_1 = pbuffer.data(idx_eri_1_gf + 69);

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

    auto g_yyyz_xyy_1 = pbuffer.data(idx_eri_1_gf + 113);

    auto g_yyyz_yyy_1 = pbuffer.data(idx_eri_1_gf + 116);

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

    auto g_yzzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 132);

    auto g_yzzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 134);

    auto g_yzzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 135);

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

    // Set up components of auxiliary buffer : HD

    auto g_xxxxx_xx_1 = pbuffer.data(idx_eri_1_hd);

    auto g_xxxxx_xy_1 = pbuffer.data(idx_eri_1_hd + 1);

    auto g_xxxxx_xz_1 = pbuffer.data(idx_eri_1_hd + 2);

    auto g_xxxxx_yy_1 = pbuffer.data(idx_eri_1_hd + 3);

    auto g_xxxxx_yz_1 = pbuffer.data(idx_eri_1_hd + 4);

    auto g_xxxxx_zz_1 = pbuffer.data(idx_eri_1_hd + 5);

    auto g_xxxxz_xz_1 = pbuffer.data(idx_eri_1_hd + 14);

    auto g_xxxxz_yz_1 = pbuffer.data(idx_eri_1_hd + 16);

    auto g_xxxxz_zz_1 = pbuffer.data(idx_eri_1_hd + 17);

    auto g_xxxyy_xx_1 = pbuffer.data(idx_eri_1_hd + 18);

    auto g_xxxyy_xy_1 = pbuffer.data(idx_eri_1_hd + 19);

    auto g_xxxyy_xz_1 = pbuffer.data(idx_eri_1_hd + 20);

    auto g_xxxyy_yy_1 = pbuffer.data(idx_eri_1_hd + 21);

    auto g_xxxyy_yz_1 = pbuffer.data(idx_eri_1_hd + 22);

    auto g_xxxyy_zz_1 = pbuffer.data(idx_eri_1_hd + 23);

    auto g_xxxzz_xx_1 = pbuffer.data(idx_eri_1_hd + 30);

    auto g_xxxzz_xy_1 = pbuffer.data(idx_eri_1_hd + 31);

    auto g_xxxzz_xz_1 = pbuffer.data(idx_eri_1_hd + 32);

    auto g_xxxzz_yy_1 = pbuffer.data(idx_eri_1_hd + 33);

    auto g_xxxzz_yz_1 = pbuffer.data(idx_eri_1_hd + 34);

    auto g_xxxzz_zz_1 = pbuffer.data(idx_eri_1_hd + 35);

    auto g_xxyyy_xx_1 = pbuffer.data(idx_eri_1_hd + 36);

    auto g_xxyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 37);

    auto g_xxyyy_xz_1 = pbuffer.data(idx_eri_1_hd + 38);

    auto g_xxyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 39);

    auto g_xxyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 40);

    auto g_xxyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 41);

    auto g_xxzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 54);

    auto g_xxzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 55);

    auto g_xxzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 56);

    auto g_xxzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 57);

    auto g_xxzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 58);

    auto g_xxzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 59);

    auto g_xyyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 61);

    auto g_xyyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 63);

    auto g_xyyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 64);

    auto g_xyyzz_yz_1 = pbuffer.data(idx_eri_1_hd + 76);

    auto g_xzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 86);

    auto g_xzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 88);

    auto g_xzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 89);

    auto g_yyyyy_xx_1 = pbuffer.data(idx_eri_1_hd + 90);

    auto g_yyyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 91);

    auto g_yyyyy_xz_1 = pbuffer.data(idx_eri_1_hd + 92);

    auto g_yyyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 93);

    auto g_yyyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 94);

    auto g_yyyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 95);

    auto g_yyyyz_xz_1 = pbuffer.data(idx_eri_1_hd + 98);

    auto g_yyyyz_yz_1 = pbuffer.data(idx_eri_1_hd + 100);

    auto g_yyyyz_zz_1 = pbuffer.data(idx_eri_1_hd + 101);

    auto g_yyyzz_xx_1 = pbuffer.data(idx_eri_1_hd + 102);

    auto g_yyyzz_xy_1 = pbuffer.data(idx_eri_1_hd + 103);

    auto g_yyyzz_xz_1 = pbuffer.data(idx_eri_1_hd + 104);

    auto g_yyyzz_yy_1 = pbuffer.data(idx_eri_1_hd + 105);

    auto g_yyyzz_yz_1 = pbuffer.data(idx_eri_1_hd + 106);

    auto g_yyyzz_zz_1 = pbuffer.data(idx_eri_1_hd + 107);

    auto g_yyzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 108);

    auto g_yyzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 109);

    auto g_yyzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 110);

    auto g_yyzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 111);

    auto g_yyzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 112);

    auto g_yyzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 113);

    auto g_yzzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 115);

    auto g_yzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 116);

    auto g_yzzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 117);

    auto g_yzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 118);

    auto g_yzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 119);

    auto g_zzzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 120);

    auto g_zzzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 121);

    auto g_zzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 122);

    auto g_zzzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 123);

    auto g_zzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 124);

    auto g_zzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 125);

    // Set up components of auxiliary buffer : HF

    auto g_xxxxx_xxx_1 = pbuffer.data(idx_eri_1_hf);

    auto g_xxxxx_xxy_1 = pbuffer.data(idx_eri_1_hf + 1);

    auto g_xxxxx_xxz_1 = pbuffer.data(idx_eri_1_hf + 2);

    auto g_xxxxx_xyy_1 = pbuffer.data(idx_eri_1_hf + 3);

    auto g_xxxxx_xyz_1 = pbuffer.data(idx_eri_1_hf + 4);

    auto g_xxxxx_xzz_1 = pbuffer.data(idx_eri_1_hf + 5);

    auto g_xxxxx_yyy_1 = pbuffer.data(idx_eri_1_hf + 6);

    auto g_xxxxx_yyz_1 = pbuffer.data(idx_eri_1_hf + 7);

    auto g_xxxxx_yzz_1 = pbuffer.data(idx_eri_1_hf + 8);

    auto g_xxxxx_zzz_1 = pbuffer.data(idx_eri_1_hf + 9);

    auto g_xxxxy_xxx_1 = pbuffer.data(idx_eri_1_hf + 10);

    auto g_xxxxy_xxy_1 = pbuffer.data(idx_eri_1_hf + 11);

    auto g_xxxxy_xxz_1 = pbuffer.data(idx_eri_1_hf + 12);

    auto g_xxxxy_xyy_1 = pbuffer.data(idx_eri_1_hf + 13);

    auto g_xxxxy_xzz_1 = pbuffer.data(idx_eri_1_hf + 15);

    auto g_xxxxy_yyy_1 = pbuffer.data(idx_eri_1_hf + 16);

    auto g_xxxxz_xxx_1 = pbuffer.data(idx_eri_1_hf + 20);

    auto g_xxxxz_xxy_1 = pbuffer.data(idx_eri_1_hf + 21);

    auto g_xxxxz_xxz_1 = pbuffer.data(idx_eri_1_hf + 22);

    auto g_xxxxz_xyy_1 = pbuffer.data(idx_eri_1_hf + 23);

    auto g_xxxxz_xyz_1 = pbuffer.data(idx_eri_1_hf + 24);

    auto g_xxxxz_xzz_1 = pbuffer.data(idx_eri_1_hf + 25);

    auto g_xxxxz_yyz_1 = pbuffer.data(idx_eri_1_hf + 27);

    auto g_xxxxz_yzz_1 = pbuffer.data(idx_eri_1_hf + 28);

    auto g_xxxxz_zzz_1 = pbuffer.data(idx_eri_1_hf + 29);

    auto g_xxxyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 30);

    auto g_xxxyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 31);

    auto g_xxxyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 32);

    auto g_xxxyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 33);

    auto g_xxxyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 34);

    auto g_xxxyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 35);

    auto g_xxxyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 36);

    auto g_xxxyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 37);

    auto g_xxxyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 38);

    auto g_xxxyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 39);

    auto g_xxxzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 50);

    auto g_xxxzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 51);

    auto g_xxxzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 52);

    auto g_xxxzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 53);

    auto g_xxxzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 54);

    auto g_xxxzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 55);

    auto g_xxxzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 56);

    auto g_xxxzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 57);

    auto g_xxxzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 58);

    auto g_xxxzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 59);

    auto g_xxyyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 60);

    auto g_xxyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 61);

    auto g_xxyyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 62);

    auto g_xxyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 63);

    auto g_xxyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 64);

    auto g_xxyyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 65);

    auto g_xxyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 66);

    auto g_xxyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 67);

    auto g_xxyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 68);

    auto g_xxyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 69);

    auto g_xxyyz_xxy_1 = pbuffer.data(idx_eri_1_hf + 71);

    auto g_xxyyz_xyy_1 = pbuffer.data(idx_eri_1_hf + 73);

    auto g_xxyzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 80);

    auto g_xxyzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 82);

    auto g_xxyzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 85);

    auto g_xxzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 90);

    auto g_xxzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 91);

    auto g_xxzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 92);

    auto g_xxzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 93);

    auto g_xxzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 94);

    auto g_xxzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 95);

    auto g_xxzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 96);

    auto g_xxzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 97);

    auto g_xxzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 98);

    auto g_xxzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 99);

    auto g_xyyyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 100);

    auto g_xyyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 101);

    auto g_xyyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 103);

    auto g_xyyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 104);

    auto g_xyyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 106);

    auto g_xyyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 107);

    auto g_xyyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 108);

    auto g_xyyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 109);

    auto g_xyyzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 124);

    auto g_xyyzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 126);

    auto g_xyyzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 127);

    auto g_xyyzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 128);

    auto g_xyyzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 129);

    auto g_xzzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 140);

    auto g_xzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 142);

    auto g_xzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 144);

    auto g_xzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 145);

    auto g_xzzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 146);

    auto g_xzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 147);

    auto g_xzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 148);

    auto g_xzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 149);

    auto g_yyyyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 150);

    auto g_yyyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 151);

    auto g_yyyyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 152);

    auto g_yyyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 153);

    auto g_yyyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 154);

    auto g_yyyyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 155);

    auto g_yyyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 156);

    auto g_yyyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 157);

    auto g_yyyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 158);

    auto g_yyyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 159);

    auto g_yyyyz_xxy_1 = pbuffer.data(idx_eri_1_hf + 161);

    auto g_yyyyz_xxz_1 = pbuffer.data(idx_eri_1_hf + 162);

    auto g_yyyyz_xyy_1 = pbuffer.data(idx_eri_1_hf + 163);

    auto g_yyyyz_xyz_1 = pbuffer.data(idx_eri_1_hf + 164);

    auto g_yyyyz_xzz_1 = pbuffer.data(idx_eri_1_hf + 165);

    auto g_yyyyz_yyy_1 = pbuffer.data(idx_eri_1_hf + 166);

    auto g_yyyyz_yyz_1 = pbuffer.data(idx_eri_1_hf + 167);

    auto g_yyyyz_yzz_1 = pbuffer.data(idx_eri_1_hf + 168);

    auto g_yyyyz_zzz_1 = pbuffer.data(idx_eri_1_hf + 169);

    auto g_yyyzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 170);

    auto g_yyyzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 171);

    auto g_yyyzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 172);

    auto g_yyyzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 173);

    auto g_yyyzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 174);

    auto g_yyyzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 175);

    auto g_yyyzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 176);

    auto g_yyyzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 177);

    auto g_yyyzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 178);

    auto g_yyyzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 179);

    auto g_yyzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 180);

    auto g_yyzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 181);

    auto g_yyzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 182);

    auto g_yyzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 183);

    auto g_yyzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 184);

    auto g_yyzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 185);

    auto g_yyzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 186);

    auto g_yyzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 187);

    auto g_yyzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 188);

    auto g_yyzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 189);

    auto g_yzzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 190);

    auto g_yzzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 191);

    auto g_yzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 192);

    auto g_yzzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 193);

    auto g_yzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 194);

    auto g_yzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 195);

    auto g_yzzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 196);

    auto g_yzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 197);

    auto g_yzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 198);

    auto g_yzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 199);

    auto g_zzzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 200);

    auto g_zzzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 201);

    auto g_zzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 202);

    auto g_zzzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 203);

    auto g_zzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 204);

    auto g_zzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 205);

    auto g_zzzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 206);

    auto g_zzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 207);

    auto g_zzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 208);

    auto g_zzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 209);

    // Set up 0-10 components of targeted buffer : IF

    auto g_xxxxxx_xxx_0 = pbuffer.data(idx_eri_0_if);

    auto g_xxxxxx_xxy_0 = pbuffer.data(idx_eri_0_if + 1);

    auto g_xxxxxx_xxz_0 = pbuffer.data(idx_eri_0_if + 2);

    auto g_xxxxxx_xyy_0 = pbuffer.data(idx_eri_0_if + 3);

    auto g_xxxxxx_xyz_0 = pbuffer.data(idx_eri_0_if + 4);

    auto g_xxxxxx_xzz_0 = pbuffer.data(idx_eri_0_if + 5);

    auto g_xxxxxx_yyy_0 = pbuffer.data(idx_eri_0_if + 6);

    auto g_xxxxxx_yyz_0 = pbuffer.data(idx_eri_0_if + 7);

    auto g_xxxxxx_yzz_0 = pbuffer.data(idx_eri_0_if + 8);

    auto g_xxxxxx_zzz_0 = pbuffer.data(idx_eri_0_if + 9);

    #pragma omp simd aligned(g_xxxx_xxx_0, g_xxxx_xxx_1, g_xxxx_xxy_0, g_xxxx_xxy_1, g_xxxx_xxz_0, g_xxxx_xxz_1, g_xxxx_xyy_0, g_xxxx_xyy_1, g_xxxx_xyz_0, g_xxxx_xyz_1, g_xxxx_xzz_0, g_xxxx_xzz_1, g_xxxx_yyy_0, g_xxxx_yyy_1, g_xxxx_yyz_0, g_xxxx_yyz_1, g_xxxx_yzz_0, g_xxxx_yzz_1, g_xxxx_zzz_0, g_xxxx_zzz_1, g_xxxxx_xx_1, g_xxxxx_xxx_1, g_xxxxx_xxy_1, g_xxxxx_xxz_1, g_xxxxx_xy_1, g_xxxxx_xyy_1, g_xxxxx_xyz_1, g_xxxxx_xz_1, g_xxxxx_xzz_1, g_xxxxx_yy_1, g_xxxxx_yyy_1, g_xxxxx_yyz_1, g_xxxxx_yz_1, g_xxxxx_yzz_1, g_xxxxx_zz_1, g_xxxxx_zzz_1, g_xxxxxx_xxx_0, g_xxxxxx_xxy_0, g_xxxxxx_xxz_0, g_xxxxxx_xyy_0, g_xxxxxx_xyz_0, g_xxxxxx_xzz_0, g_xxxxxx_yyy_0, g_xxxxxx_yyz_0, g_xxxxxx_yzz_0, g_xxxxxx_zzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxx_xxx_0[i] = 5.0 * g_xxxx_xxx_0[i] * fbe_0 - 5.0 * g_xxxx_xxx_1[i] * fz_be_0 + 3.0 * g_xxxxx_xx_1[i] * fe_0 + g_xxxxx_xxx_1[i] * pa_x[i];

        g_xxxxxx_xxy_0[i] = 5.0 * g_xxxx_xxy_0[i] * fbe_0 - 5.0 * g_xxxx_xxy_1[i] * fz_be_0 + 2.0 * g_xxxxx_xy_1[i] * fe_0 + g_xxxxx_xxy_1[i] * pa_x[i];

        g_xxxxxx_xxz_0[i] = 5.0 * g_xxxx_xxz_0[i] * fbe_0 - 5.0 * g_xxxx_xxz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xz_1[i] * fe_0 + g_xxxxx_xxz_1[i] * pa_x[i];

        g_xxxxxx_xyy_0[i] = 5.0 * g_xxxx_xyy_0[i] * fbe_0 - 5.0 * g_xxxx_xyy_1[i] * fz_be_0 + g_xxxxx_yy_1[i] * fe_0 + g_xxxxx_xyy_1[i] * pa_x[i];

        g_xxxxxx_xyz_0[i] = 5.0 * g_xxxx_xyz_0[i] * fbe_0 - 5.0 * g_xxxx_xyz_1[i] * fz_be_0 + g_xxxxx_yz_1[i] * fe_0 + g_xxxxx_xyz_1[i] * pa_x[i];

        g_xxxxxx_xzz_0[i] = 5.0 * g_xxxx_xzz_0[i] * fbe_0 - 5.0 * g_xxxx_xzz_1[i] * fz_be_0 + g_xxxxx_zz_1[i] * fe_0 + g_xxxxx_xzz_1[i] * pa_x[i];

        g_xxxxxx_yyy_0[i] = 5.0 * g_xxxx_yyy_0[i] * fbe_0 - 5.0 * g_xxxx_yyy_1[i] * fz_be_0 + g_xxxxx_yyy_1[i] * pa_x[i];

        g_xxxxxx_yyz_0[i] = 5.0 * g_xxxx_yyz_0[i] * fbe_0 - 5.0 * g_xxxx_yyz_1[i] * fz_be_0 + g_xxxxx_yyz_1[i] * pa_x[i];

        g_xxxxxx_yzz_0[i] = 5.0 * g_xxxx_yzz_0[i] * fbe_0 - 5.0 * g_xxxx_yzz_1[i] * fz_be_0 + g_xxxxx_yzz_1[i] * pa_x[i];

        g_xxxxxx_zzz_0[i] = 5.0 * g_xxxx_zzz_0[i] * fbe_0 - 5.0 * g_xxxx_zzz_1[i] * fz_be_0 + g_xxxxx_zzz_1[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : IF

    auto g_xxxxxy_xxx_0 = pbuffer.data(idx_eri_0_if + 10);

    auto g_xxxxxy_xxy_0 = pbuffer.data(idx_eri_0_if + 11);

    auto g_xxxxxy_xxz_0 = pbuffer.data(idx_eri_0_if + 12);

    auto g_xxxxxy_xyy_0 = pbuffer.data(idx_eri_0_if + 13);

    auto g_xxxxxy_xyz_0 = pbuffer.data(idx_eri_0_if + 14);

    auto g_xxxxxy_xzz_0 = pbuffer.data(idx_eri_0_if + 15);

    auto g_xxxxxy_yyy_0 = pbuffer.data(idx_eri_0_if + 16);

    auto g_xxxxxy_yyz_0 = pbuffer.data(idx_eri_0_if + 17);

    auto g_xxxxxy_yzz_0 = pbuffer.data(idx_eri_0_if + 18);

    auto g_xxxxxy_zzz_0 = pbuffer.data(idx_eri_0_if + 19);

    #pragma omp simd aligned(g_xxxxx_xx_1, g_xxxxx_xxx_1, g_xxxxx_xxy_1, g_xxxxx_xxz_1, g_xxxxx_xy_1, g_xxxxx_xyy_1, g_xxxxx_xyz_1, g_xxxxx_xz_1, g_xxxxx_xzz_1, g_xxxxx_yy_1, g_xxxxx_yyy_1, g_xxxxx_yyz_1, g_xxxxx_yz_1, g_xxxxx_yzz_1, g_xxxxx_zz_1, g_xxxxx_zzz_1, g_xxxxxy_xxx_0, g_xxxxxy_xxy_0, g_xxxxxy_xxz_0, g_xxxxxy_xyy_0, g_xxxxxy_xyz_0, g_xxxxxy_xzz_0, g_xxxxxy_yyy_0, g_xxxxxy_yyz_0, g_xxxxxy_yzz_0, g_xxxxxy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxy_xxx_0[i] = g_xxxxx_xxx_1[i] * pa_y[i];

        g_xxxxxy_xxy_0[i] = g_xxxxx_xx_1[i] * fe_0 + g_xxxxx_xxy_1[i] * pa_y[i];

        g_xxxxxy_xxz_0[i] = g_xxxxx_xxz_1[i] * pa_y[i];

        g_xxxxxy_xyy_0[i] = 2.0 * g_xxxxx_xy_1[i] * fe_0 + g_xxxxx_xyy_1[i] * pa_y[i];

        g_xxxxxy_xyz_0[i] = g_xxxxx_xz_1[i] * fe_0 + g_xxxxx_xyz_1[i] * pa_y[i];

        g_xxxxxy_xzz_0[i] = g_xxxxx_xzz_1[i] * pa_y[i];

        g_xxxxxy_yyy_0[i] = 3.0 * g_xxxxx_yy_1[i] * fe_0 + g_xxxxx_yyy_1[i] * pa_y[i];

        g_xxxxxy_yyz_0[i] = 2.0 * g_xxxxx_yz_1[i] * fe_0 + g_xxxxx_yyz_1[i] * pa_y[i];

        g_xxxxxy_yzz_0[i] = g_xxxxx_zz_1[i] * fe_0 + g_xxxxx_yzz_1[i] * pa_y[i];

        g_xxxxxy_zzz_0[i] = g_xxxxx_zzz_1[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : IF

    auto g_xxxxxz_xxx_0 = pbuffer.data(idx_eri_0_if + 20);

    auto g_xxxxxz_xxy_0 = pbuffer.data(idx_eri_0_if + 21);

    auto g_xxxxxz_xxz_0 = pbuffer.data(idx_eri_0_if + 22);

    auto g_xxxxxz_xyy_0 = pbuffer.data(idx_eri_0_if + 23);

    auto g_xxxxxz_xyz_0 = pbuffer.data(idx_eri_0_if + 24);

    auto g_xxxxxz_xzz_0 = pbuffer.data(idx_eri_0_if + 25);

    auto g_xxxxxz_yyy_0 = pbuffer.data(idx_eri_0_if + 26);

    auto g_xxxxxz_yyz_0 = pbuffer.data(idx_eri_0_if + 27);

    auto g_xxxxxz_yzz_0 = pbuffer.data(idx_eri_0_if + 28);

    auto g_xxxxxz_zzz_0 = pbuffer.data(idx_eri_0_if + 29);

    #pragma omp simd aligned(g_xxxxx_xx_1, g_xxxxx_xxx_1, g_xxxxx_xxy_1, g_xxxxx_xxz_1, g_xxxxx_xy_1, g_xxxxx_xyy_1, g_xxxxx_xyz_1, g_xxxxx_xz_1, g_xxxxx_xzz_1, g_xxxxx_yy_1, g_xxxxx_yyy_1, g_xxxxx_yyz_1, g_xxxxx_yz_1, g_xxxxx_yzz_1, g_xxxxx_zz_1, g_xxxxx_zzz_1, g_xxxxxz_xxx_0, g_xxxxxz_xxy_0, g_xxxxxz_xxz_0, g_xxxxxz_xyy_0, g_xxxxxz_xyz_0, g_xxxxxz_xzz_0, g_xxxxxz_yyy_0, g_xxxxxz_yyz_0, g_xxxxxz_yzz_0, g_xxxxxz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxz_xxx_0[i] = g_xxxxx_xxx_1[i] * pa_z[i];

        g_xxxxxz_xxy_0[i] = g_xxxxx_xxy_1[i] * pa_z[i];

        g_xxxxxz_xxz_0[i] = g_xxxxx_xx_1[i] * fe_0 + g_xxxxx_xxz_1[i] * pa_z[i];

        g_xxxxxz_xyy_0[i] = g_xxxxx_xyy_1[i] * pa_z[i];

        g_xxxxxz_xyz_0[i] = g_xxxxx_xy_1[i] * fe_0 + g_xxxxx_xyz_1[i] * pa_z[i];

        g_xxxxxz_xzz_0[i] = 2.0 * g_xxxxx_xz_1[i] * fe_0 + g_xxxxx_xzz_1[i] * pa_z[i];

        g_xxxxxz_yyy_0[i] = g_xxxxx_yyy_1[i] * pa_z[i];

        g_xxxxxz_yyz_0[i] = g_xxxxx_yy_1[i] * fe_0 + g_xxxxx_yyz_1[i] * pa_z[i];

        g_xxxxxz_yzz_0[i] = 2.0 * g_xxxxx_yz_1[i] * fe_0 + g_xxxxx_yzz_1[i] * pa_z[i];

        g_xxxxxz_zzz_0[i] = 3.0 * g_xxxxx_zz_1[i] * fe_0 + g_xxxxx_zzz_1[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : IF

    auto g_xxxxyy_xxx_0 = pbuffer.data(idx_eri_0_if + 30);

    auto g_xxxxyy_xxy_0 = pbuffer.data(idx_eri_0_if + 31);

    auto g_xxxxyy_xxz_0 = pbuffer.data(idx_eri_0_if + 32);

    auto g_xxxxyy_xyy_0 = pbuffer.data(idx_eri_0_if + 33);

    auto g_xxxxyy_xyz_0 = pbuffer.data(idx_eri_0_if + 34);

    auto g_xxxxyy_xzz_0 = pbuffer.data(idx_eri_0_if + 35);

    auto g_xxxxyy_yyy_0 = pbuffer.data(idx_eri_0_if + 36);

    auto g_xxxxyy_yyz_0 = pbuffer.data(idx_eri_0_if + 37);

    auto g_xxxxyy_yzz_0 = pbuffer.data(idx_eri_0_if + 38);

    auto g_xxxxyy_zzz_0 = pbuffer.data(idx_eri_0_if + 39);

    #pragma omp simd aligned(g_xxxx_xxx_0, g_xxxx_xxx_1, g_xxxx_xxz_0, g_xxxx_xxz_1, g_xxxx_xzz_0, g_xxxx_xzz_1, g_xxxxy_xxx_1, g_xxxxy_xxz_1, g_xxxxy_xzz_1, g_xxxxyy_xxx_0, g_xxxxyy_xxy_0, g_xxxxyy_xxz_0, g_xxxxyy_xyy_0, g_xxxxyy_xyz_0, g_xxxxyy_xzz_0, g_xxxxyy_yyy_0, g_xxxxyy_yyz_0, g_xxxxyy_yzz_0, g_xxxxyy_zzz_0, g_xxxyy_xxy_1, g_xxxyy_xy_1, g_xxxyy_xyy_1, g_xxxyy_xyz_1, g_xxxyy_yy_1, g_xxxyy_yyy_1, g_xxxyy_yyz_1, g_xxxyy_yz_1, g_xxxyy_yzz_1, g_xxxyy_zzz_1, g_xxyy_xxy_0, g_xxyy_xxy_1, g_xxyy_xyy_0, g_xxyy_xyy_1, g_xxyy_xyz_0, g_xxyy_xyz_1, g_xxyy_yyy_0, g_xxyy_yyy_1, g_xxyy_yyz_0, g_xxyy_yyz_1, g_xxyy_yzz_0, g_xxyy_yzz_1, g_xxyy_zzz_0, g_xxyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyy_xxx_0[i] = g_xxxx_xxx_0[i] * fbe_0 - g_xxxx_xxx_1[i] * fz_be_0 + g_xxxxy_xxx_1[i] * pa_y[i];

        g_xxxxyy_xxy_0[i] = 3.0 * g_xxyy_xxy_0[i] * fbe_0 - 3.0 * g_xxyy_xxy_1[i] * fz_be_0 + 2.0 * g_xxxyy_xy_1[i] * fe_0 + g_xxxyy_xxy_1[i] * pa_x[i];

        g_xxxxyy_xxz_0[i] = g_xxxx_xxz_0[i] * fbe_0 - g_xxxx_xxz_1[i] * fz_be_0 + g_xxxxy_xxz_1[i] * pa_y[i];

        g_xxxxyy_xyy_0[i] = 3.0 * g_xxyy_xyy_0[i] * fbe_0 - 3.0 * g_xxyy_xyy_1[i] * fz_be_0 + g_xxxyy_yy_1[i] * fe_0 + g_xxxyy_xyy_1[i] * pa_x[i];

        g_xxxxyy_xyz_0[i] = 3.0 * g_xxyy_xyz_0[i] * fbe_0 - 3.0 * g_xxyy_xyz_1[i] * fz_be_0 + g_xxxyy_yz_1[i] * fe_0 + g_xxxyy_xyz_1[i] * pa_x[i];

        g_xxxxyy_xzz_0[i] = g_xxxx_xzz_0[i] * fbe_0 - g_xxxx_xzz_1[i] * fz_be_0 + g_xxxxy_xzz_1[i] * pa_y[i];

        g_xxxxyy_yyy_0[i] = 3.0 * g_xxyy_yyy_0[i] * fbe_0 - 3.0 * g_xxyy_yyy_1[i] * fz_be_0 + g_xxxyy_yyy_1[i] * pa_x[i];

        g_xxxxyy_yyz_0[i] = 3.0 * g_xxyy_yyz_0[i] * fbe_0 - 3.0 * g_xxyy_yyz_1[i] * fz_be_0 + g_xxxyy_yyz_1[i] * pa_x[i];

        g_xxxxyy_yzz_0[i] = 3.0 * g_xxyy_yzz_0[i] * fbe_0 - 3.0 * g_xxyy_yzz_1[i] * fz_be_0 + g_xxxyy_yzz_1[i] * pa_x[i];

        g_xxxxyy_zzz_0[i] = 3.0 * g_xxyy_zzz_0[i] * fbe_0 - 3.0 * g_xxyy_zzz_1[i] * fz_be_0 + g_xxxyy_zzz_1[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : IF

    auto g_xxxxyz_xxx_0 = pbuffer.data(idx_eri_0_if + 40);

    auto g_xxxxyz_xxy_0 = pbuffer.data(idx_eri_0_if + 41);

    auto g_xxxxyz_xxz_0 = pbuffer.data(idx_eri_0_if + 42);

    auto g_xxxxyz_xyy_0 = pbuffer.data(idx_eri_0_if + 43);

    auto g_xxxxyz_xyz_0 = pbuffer.data(idx_eri_0_if + 44);

    auto g_xxxxyz_xzz_0 = pbuffer.data(idx_eri_0_if + 45);

    auto g_xxxxyz_yyy_0 = pbuffer.data(idx_eri_0_if + 46);

    auto g_xxxxyz_yyz_0 = pbuffer.data(idx_eri_0_if + 47);

    auto g_xxxxyz_yzz_0 = pbuffer.data(idx_eri_0_if + 48);

    auto g_xxxxyz_zzz_0 = pbuffer.data(idx_eri_0_if + 49);

    #pragma omp simd aligned(g_xxxxy_xxy_1, g_xxxxy_xyy_1, g_xxxxy_yyy_1, g_xxxxyz_xxx_0, g_xxxxyz_xxy_0, g_xxxxyz_xxz_0, g_xxxxyz_xyy_0, g_xxxxyz_xyz_0, g_xxxxyz_xzz_0, g_xxxxyz_yyy_0, g_xxxxyz_yyz_0, g_xxxxyz_yzz_0, g_xxxxyz_zzz_0, g_xxxxz_xxx_1, g_xxxxz_xxz_1, g_xxxxz_xyz_1, g_xxxxz_xz_1, g_xxxxz_xzz_1, g_xxxxz_yyz_1, g_xxxxz_yz_1, g_xxxxz_yzz_1, g_xxxxz_zz_1, g_xxxxz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyz_xxx_0[i] = g_xxxxz_xxx_1[i] * pa_y[i];

        g_xxxxyz_xxy_0[i] = g_xxxxy_xxy_1[i] * pa_z[i];

        g_xxxxyz_xxz_0[i] = g_xxxxz_xxz_1[i] * pa_y[i];

        g_xxxxyz_xyy_0[i] = g_xxxxy_xyy_1[i] * pa_z[i];

        g_xxxxyz_xyz_0[i] = g_xxxxz_xz_1[i] * fe_0 + g_xxxxz_xyz_1[i] * pa_y[i];

        g_xxxxyz_xzz_0[i] = g_xxxxz_xzz_1[i] * pa_y[i];

        g_xxxxyz_yyy_0[i] = g_xxxxy_yyy_1[i] * pa_z[i];

        g_xxxxyz_yyz_0[i] = 2.0 * g_xxxxz_yz_1[i] * fe_0 + g_xxxxz_yyz_1[i] * pa_y[i];

        g_xxxxyz_yzz_0[i] = g_xxxxz_zz_1[i] * fe_0 + g_xxxxz_yzz_1[i] * pa_y[i];

        g_xxxxyz_zzz_0[i] = g_xxxxz_zzz_1[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : IF

    auto g_xxxxzz_xxx_0 = pbuffer.data(idx_eri_0_if + 50);

    auto g_xxxxzz_xxy_0 = pbuffer.data(idx_eri_0_if + 51);

    auto g_xxxxzz_xxz_0 = pbuffer.data(idx_eri_0_if + 52);

    auto g_xxxxzz_xyy_0 = pbuffer.data(idx_eri_0_if + 53);

    auto g_xxxxzz_xyz_0 = pbuffer.data(idx_eri_0_if + 54);

    auto g_xxxxzz_xzz_0 = pbuffer.data(idx_eri_0_if + 55);

    auto g_xxxxzz_yyy_0 = pbuffer.data(idx_eri_0_if + 56);

    auto g_xxxxzz_yyz_0 = pbuffer.data(idx_eri_0_if + 57);

    auto g_xxxxzz_yzz_0 = pbuffer.data(idx_eri_0_if + 58);

    auto g_xxxxzz_zzz_0 = pbuffer.data(idx_eri_0_if + 59);

    #pragma omp simd aligned(g_xxxx_xxx_0, g_xxxx_xxx_1, g_xxxx_xxy_0, g_xxxx_xxy_1, g_xxxx_xyy_0, g_xxxx_xyy_1, g_xxxxz_xxx_1, g_xxxxz_xxy_1, g_xxxxz_xyy_1, g_xxxxzz_xxx_0, g_xxxxzz_xxy_0, g_xxxxzz_xxz_0, g_xxxxzz_xyy_0, g_xxxxzz_xyz_0, g_xxxxzz_xzz_0, g_xxxxzz_yyy_0, g_xxxxzz_yyz_0, g_xxxxzz_yzz_0, g_xxxxzz_zzz_0, g_xxxzz_xxz_1, g_xxxzz_xyz_1, g_xxxzz_xz_1, g_xxxzz_xzz_1, g_xxxzz_yyy_1, g_xxxzz_yyz_1, g_xxxzz_yz_1, g_xxxzz_yzz_1, g_xxxzz_zz_1, g_xxxzz_zzz_1, g_xxzz_xxz_0, g_xxzz_xxz_1, g_xxzz_xyz_0, g_xxzz_xyz_1, g_xxzz_xzz_0, g_xxzz_xzz_1, g_xxzz_yyy_0, g_xxzz_yyy_1, g_xxzz_yyz_0, g_xxzz_yyz_1, g_xxzz_yzz_0, g_xxzz_yzz_1, g_xxzz_zzz_0, g_xxzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzz_xxx_0[i] = g_xxxx_xxx_0[i] * fbe_0 - g_xxxx_xxx_1[i] * fz_be_0 + g_xxxxz_xxx_1[i] * pa_z[i];

        g_xxxxzz_xxy_0[i] = g_xxxx_xxy_0[i] * fbe_0 - g_xxxx_xxy_1[i] * fz_be_0 + g_xxxxz_xxy_1[i] * pa_z[i];

        g_xxxxzz_xxz_0[i] = 3.0 * g_xxzz_xxz_0[i] * fbe_0 - 3.0 * g_xxzz_xxz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xz_1[i] * fe_0 + g_xxxzz_xxz_1[i] * pa_x[i];

        g_xxxxzz_xyy_0[i] = g_xxxx_xyy_0[i] * fbe_0 - g_xxxx_xyy_1[i] * fz_be_0 + g_xxxxz_xyy_1[i] * pa_z[i];

        g_xxxxzz_xyz_0[i] = 3.0 * g_xxzz_xyz_0[i] * fbe_0 - 3.0 * g_xxzz_xyz_1[i] * fz_be_0 + g_xxxzz_yz_1[i] * fe_0 + g_xxxzz_xyz_1[i] * pa_x[i];

        g_xxxxzz_xzz_0[i] = 3.0 * g_xxzz_xzz_0[i] * fbe_0 - 3.0 * g_xxzz_xzz_1[i] * fz_be_0 + g_xxxzz_zz_1[i] * fe_0 + g_xxxzz_xzz_1[i] * pa_x[i];

        g_xxxxzz_yyy_0[i] = 3.0 * g_xxzz_yyy_0[i] * fbe_0 - 3.0 * g_xxzz_yyy_1[i] * fz_be_0 + g_xxxzz_yyy_1[i] * pa_x[i];

        g_xxxxzz_yyz_0[i] = 3.0 * g_xxzz_yyz_0[i] * fbe_0 - 3.0 * g_xxzz_yyz_1[i] * fz_be_0 + g_xxxzz_yyz_1[i] * pa_x[i];

        g_xxxxzz_yzz_0[i] = 3.0 * g_xxzz_yzz_0[i] * fbe_0 - 3.0 * g_xxzz_yzz_1[i] * fz_be_0 + g_xxxzz_yzz_1[i] * pa_x[i];

        g_xxxxzz_zzz_0[i] = 3.0 * g_xxzz_zzz_0[i] * fbe_0 - 3.0 * g_xxzz_zzz_1[i] * fz_be_0 + g_xxxzz_zzz_1[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : IF

    auto g_xxxyyy_xxx_0 = pbuffer.data(idx_eri_0_if + 60);

    auto g_xxxyyy_xxy_0 = pbuffer.data(idx_eri_0_if + 61);

    auto g_xxxyyy_xxz_0 = pbuffer.data(idx_eri_0_if + 62);

    auto g_xxxyyy_xyy_0 = pbuffer.data(idx_eri_0_if + 63);

    auto g_xxxyyy_xyz_0 = pbuffer.data(idx_eri_0_if + 64);

    auto g_xxxyyy_xzz_0 = pbuffer.data(idx_eri_0_if + 65);

    auto g_xxxyyy_yyy_0 = pbuffer.data(idx_eri_0_if + 66);

    auto g_xxxyyy_yyz_0 = pbuffer.data(idx_eri_0_if + 67);

    auto g_xxxyyy_yzz_0 = pbuffer.data(idx_eri_0_if + 68);

    auto g_xxxyyy_zzz_0 = pbuffer.data(idx_eri_0_if + 69);

    #pragma omp simd aligned(g_xxxy_xxx_0, g_xxxy_xxx_1, g_xxxy_xxz_0, g_xxxy_xxz_1, g_xxxy_xzz_0, g_xxxy_xzz_1, g_xxxyy_xxx_1, g_xxxyy_xxz_1, g_xxxyy_xzz_1, g_xxxyyy_xxx_0, g_xxxyyy_xxy_0, g_xxxyyy_xxz_0, g_xxxyyy_xyy_0, g_xxxyyy_xyz_0, g_xxxyyy_xzz_0, g_xxxyyy_yyy_0, g_xxxyyy_yyz_0, g_xxxyyy_yzz_0, g_xxxyyy_zzz_0, g_xxyyy_xxy_1, g_xxyyy_xy_1, g_xxyyy_xyy_1, g_xxyyy_xyz_1, g_xxyyy_yy_1, g_xxyyy_yyy_1, g_xxyyy_yyz_1, g_xxyyy_yz_1, g_xxyyy_yzz_1, g_xxyyy_zzz_1, g_xyyy_xxy_0, g_xyyy_xxy_1, g_xyyy_xyy_0, g_xyyy_xyy_1, g_xyyy_xyz_0, g_xyyy_xyz_1, g_xyyy_yyy_0, g_xyyy_yyy_1, g_xyyy_yyz_0, g_xyyy_yyz_1, g_xyyy_yzz_0, g_xyyy_yzz_1, g_xyyy_zzz_0, g_xyyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyy_xxx_0[i] = 2.0 * g_xxxy_xxx_0[i] * fbe_0 - 2.0 * g_xxxy_xxx_1[i] * fz_be_0 + g_xxxyy_xxx_1[i] * pa_y[i];

        g_xxxyyy_xxy_0[i] = 2.0 * g_xyyy_xxy_0[i] * fbe_0 - 2.0 * g_xyyy_xxy_1[i] * fz_be_0 + 2.0 * g_xxyyy_xy_1[i] * fe_0 + g_xxyyy_xxy_1[i] * pa_x[i];

        g_xxxyyy_xxz_0[i] = 2.0 * g_xxxy_xxz_0[i] * fbe_0 - 2.0 * g_xxxy_xxz_1[i] * fz_be_0 + g_xxxyy_xxz_1[i] * pa_y[i];

        g_xxxyyy_xyy_0[i] = 2.0 * g_xyyy_xyy_0[i] * fbe_0 - 2.0 * g_xyyy_xyy_1[i] * fz_be_0 + g_xxyyy_yy_1[i] * fe_0 + g_xxyyy_xyy_1[i] * pa_x[i];

        g_xxxyyy_xyz_0[i] = 2.0 * g_xyyy_xyz_0[i] * fbe_0 - 2.0 * g_xyyy_xyz_1[i] * fz_be_0 + g_xxyyy_yz_1[i] * fe_0 + g_xxyyy_xyz_1[i] * pa_x[i];

        g_xxxyyy_xzz_0[i] = 2.0 * g_xxxy_xzz_0[i] * fbe_0 - 2.0 * g_xxxy_xzz_1[i] * fz_be_0 + g_xxxyy_xzz_1[i] * pa_y[i];

        g_xxxyyy_yyy_0[i] = 2.0 * g_xyyy_yyy_0[i] * fbe_0 - 2.0 * g_xyyy_yyy_1[i] * fz_be_0 + g_xxyyy_yyy_1[i] * pa_x[i];

        g_xxxyyy_yyz_0[i] = 2.0 * g_xyyy_yyz_0[i] * fbe_0 - 2.0 * g_xyyy_yyz_1[i] * fz_be_0 + g_xxyyy_yyz_1[i] * pa_x[i];

        g_xxxyyy_yzz_0[i] = 2.0 * g_xyyy_yzz_0[i] * fbe_0 - 2.0 * g_xyyy_yzz_1[i] * fz_be_0 + g_xxyyy_yzz_1[i] * pa_x[i];

        g_xxxyyy_zzz_0[i] = 2.0 * g_xyyy_zzz_0[i] * fbe_0 - 2.0 * g_xyyy_zzz_1[i] * fz_be_0 + g_xxyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : IF

    auto g_xxxyyz_xxx_0 = pbuffer.data(idx_eri_0_if + 70);

    auto g_xxxyyz_xxy_0 = pbuffer.data(idx_eri_0_if + 71);

    auto g_xxxyyz_xxz_0 = pbuffer.data(idx_eri_0_if + 72);

    auto g_xxxyyz_xyy_0 = pbuffer.data(idx_eri_0_if + 73);

    auto g_xxxyyz_xyz_0 = pbuffer.data(idx_eri_0_if + 74);

    auto g_xxxyyz_xzz_0 = pbuffer.data(idx_eri_0_if + 75);

    auto g_xxxyyz_yyy_0 = pbuffer.data(idx_eri_0_if + 76);

    auto g_xxxyyz_yyz_0 = pbuffer.data(idx_eri_0_if + 77);

    auto g_xxxyyz_yzz_0 = pbuffer.data(idx_eri_0_if + 78);

    auto g_xxxyyz_zzz_0 = pbuffer.data(idx_eri_0_if + 79);

    #pragma omp simd aligned(g_xxxyy_xx_1, g_xxxyy_xxx_1, g_xxxyy_xxy_1, g_xxxyy_xxz_1, g_xxxyy_xy_1, g_xxxyy_xyy_1, g_xxxyy_xyz_1, g_xxxyy_xz_1, g_xxxyy_xzz_1, g_xxxyy_yy_1, g_xxxyy_yyy_1, g_xxxyy_yyz_1, g_xxxyy_yz_1, g_xxxyy_yzz_1, g_xxxyy_zz_1, g_xxxyy_zzz_1, g_xxxyyz_xxx_0, g_xxxyyz_xxy_0, g_xxxyyz_xxz_0, g_xxxyyz_xyy_0, g_xxxyyz_xyz_0, g_xxxyyz_xzz_0, g_xxxyyz_yyy_0, g_xxxyyz_yyz_0, g_xxxyyz_yzz_0, g_xxxyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyz_xxx_0[i] = g_xxxyy_xxx_1[i] * pa_z[i];

        g_xxxyyz_xxy_0[i] = g_xxxyy_xxy_1[i] * pa_z[i];

        g_xxxyyz_xxz_0[i] = g_xxxyy_xx_1[i] * fe_0 + g_xxxyy_xxz_1[i] * pa_z[i];

        g_xxxyyz_xyy_0[i] = g_xxxyy_xyy_1[i] * pa_z[i];

        g_xxxyyz_xyz_0[i] = g_xxxyy_xy_1[i] * fe_0 + g_xxxyy_xyz_1[i] * pa_z[i];

        g_xxxyyz_xzz_0[i] = 2.0 * g_xxxyy_xz_1[i] * fe_0 + g_xxxyy_xzz_1[i] * pa_z[i];

        g_xxxyyz_yyy_0[i] = g_xxxyy_yyy_1[i] * pa_z[i];

        g_xxxyyz_yyz_0[i] = g_xxxyy_yy_1[i] * fe_0 + g_xxxyy_yyz_1[i] * pa_z[i];

        g_xxxyyz_yzz_0[i] = 2.0 * g_xxxyy_yz_1[i] * fe_0 + g_xxxyy_yzz_1[i] * pa_z[i];

        g_xxxyyz_zzz_0[i] = 3.0 * g_xxxyy_zz_1[i] * fe_0 + g_xxxyy_zzz_1[i] * pa_z[i];
    }

    // Set up 80-90 components of targeted buffer : IF

    auto g_xxxyzz_xxx_0 = pbuffer.data(idx_eri_0_if + 80);

    auto g_xxxyzz_xxy_0 = pbuffer.data(idx_eri_0_if + 81);

    auto g_xxxyzz_xxz_0 = pbuffer.data(idx_eri_0_if + 82);

    auto g_xxxyzz_xyy_0 = pbuffer.data(idx_eri_0_if + 83);

    auto g_xxxyzz_xyz_0 = pbuffer.data(idx_eri_0_if + 84);

    auto g_xxxyzz_xzz_0 = pbuffer.data(idx_eri_0_if + 85);

    auto g_xxxyzz_yyy_0 = pbuffer.data(idx_eri_0_if + 86);

    auto g_xxxyzz_yyz_0 = pbuffer.data(idx_eri_0_if + 87);

    auto g_xxxyzz_yzz_0 = pbuffer.data(idx_eri_0_if + 88);

    auto g_xxxyzz_zzz_0 = pbuffer.data(idx_eri_0_if + 89);

    #pragma omp simd aligned(g_xxxyzz_xxx_0, g_xxxyzz_xxy_0, g_xxxyzz_xxz_0, g_xxxyzz_xyy_0, g_xxxyzz_xyz_0, g_xxxyzz_xzz_0, g_xxxyzz_yyy_0, g_xxxyzz_yyz_0, g_xxxyzz_yzz_0, g_xxxyzz_zzz_0, g_xxxzz_xx_1, g_xxxzz_xxx_1, g_xxxzz_xxy_1, g_xxxzz_xxz_1, g_xxxzz_xy_1, g_xxxzz_xyy_1, g_xxxzz_xyz_1, g_xxxzz_xz_1, g_xxxzz_xzz_1, g_xxxzz_yy_1, g_xxxzz_yyy_1, g_xxxzz_yyz_1, g_xxxzz_yz_1, g_xxxzz_yzz_1, g_xxxzz_zz_1, g_xxxzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzz_xxx_0[i] = g_xxxzz_xxx_1[i] * pa_y[i];

        g_xxxyzz_xxy_0[i] = g_xxxzz_xx_1[i] * fe_0 + g_xxxzz_xxy_1[i] * pa_y[i];

        g_xxxyzz_xxz_0[i] = g_xxxzz_xxz_1[i] * pa_y[i];

        g_xxxyzz_xyy_0[i] = 2.0 * g_xxxzz_xy_1[i] * fe_0 + g_xxxzz_xyy_1[i] * pa_y[i];

        g_xxxyzz_xyz_0[i] = g_xxxzz_xz_1[i] * fe_0 + g_xxxzz_xyz_1[i] * pa_y[i];

        g_xxxyzz_xzz_0[i] = g_xxxzz_xzz_1[i] * pa_y[i];

        g_xxxyzz_yyy_0[i] = 3.0 * g_xxxzz_yy_1[i] * fe_0 + g_xxxzz_yyy_1[i] * pa_y[i];

        g_xxxyzz_yyz_0[i] = 2.0 * g_xxxzz_yz_1[i] * fe_0 + g_xxxzz_yyz_1[i] * pa_y[i];

        g_xxxyzz_yzz_0[i] = g_xxxzz_zz_1[i] * fe_0 + g_xxxzz_yzz_1[i] * pa_y[i];

        g_xxxyzz_zzz_0[i] = g_xxxzz_zzz_1[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : IF

    auto g_xxxzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 90);

    auto g_xxxzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 91);

    auto g_xxxzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 92);

    auto g_xxxzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 93);

    auto g_xxxzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 94);

    auto g_xxxzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 95);

    auto g_xxxzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 96);

    auto g_xxxzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 97);

    auto g_xxxzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 98);

    auto g_xxxzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 99);

    #pragma omp simd aligned(g_xxxz_xxx_0, g_xxxz_xxx_1, g_xxxz_xxy_0, g_xxxz_xxy_1, g_xxxz_xyy_0, g_xxxz_xyy_1, g_xxxzz_xxx_1, g_xxxzz_xxy_1, g_xxxzz_xyy_1, g_xxxzzz_xxx_0, g_xxxzzz_xxy_0, g_xxxzzz_xxz_0, g_xxxzzz_xyy_0, g_xxxzzz_xyz_0, g_xxxzzz_xzz_0, g_xxxzzz_yyy_0, g_xxxzzz_yyz_0, g_xxxzzz_yzz_0, g_xxxzzz_zzz_0, g_xxzzz_xxz_1, g_xxzzz_xyz_1, g_xxzzz_xz_1, g_xxzzz_xzz_1, g_xxzzz_yyy_1, g_xxzzz_yyz_1, g_xxzzz_yz_1, g_xxzzz_yzz_1, g_xxzzz_zz_1, g_xxzzz_zzz_1, g_xzzz_xxz_0, g_xzzz_xxz_1, g_xzzz_xyz_0, g_xzzz_xyz_1, g_xzzz_xzz_0, g_xzzz_xzz_1, g_xzzz_yyy_0, g_xzzz_yyy_1, g_xzzz_yyz_0, g_xzzz_yyz_1, g_xzzz_yzz_0, g_xzzz_yzz_1, g_xzzz_zzz_0, g_xzzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzz_xxx_0[i] = 2.0 * g_xxxz_xxx_0[i] * fbe_0 - 2.0 * g_xxxz_xxx_1[i] * fz_be_0 + g_xxxzz_xxx_1[i] * pa_z[i];

        g_xxxzzz_xxy_0[i] = 2.0 * g_xxxz_xxy_0[i] * fbe_0 - 2.0 * g_xxxz_xxy_1[i] * fz_be_0 + g_xxxzz_xxy_1[i] * pa_z[i];

        g_xxxzzz_xxz_0[i] = 2.0 * g_xzzz_xxz_0[i] * fbe_0 - 2.0 * g_xzzz_xxz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xz_1[i] * fe_0 + g_xxzzz_xxz_1[i] * pa_x[i];

        g_xxxzzz_xyy_0[i] = 2.0 * g_xxxz_xyy_0[i] * fbe_0 - 2.0 * g_xxxz_xyy_1[i] * fz_be_0 + g_xxxzz_xyy_1[i] * pa_z[i];

        g_xxxzzz_xyz_0[i] = 2.0 * g_xzzz_xyz_0[i] * fbe_0 - 2.0 * g_xzzz_xyz_1[i] * fz_be_0 + g_xxzzz_yz_1[i] * fe_0 + g_xxzzz_xyz_1[i] * pa_x[i];

        g_xxxzzz_xzz_0[i] = 2.0 * g_xzzz_xzz_0[i] * fbe_0 - 2.0 * g_xzzz_xzz_1[i] * fz_be_0 + g_xxzzz_zz_1[i] * fe_0 + g_xxzzz_xzz_1[i] * pa_x[i];

        g_xxxzzz_yyy_0[i] = 2.0 * g_xzzz_yyy_0[i] * fbe_0 - 2.0 * g_xzzz_yyy_1[i] * fz_be_0 + g_xxzzz_yyy_1[i] * pa_x[i];

        g_xxxzzz_yyz_0[i] = 2.0 * g_xzzz_yyz_0[i] * fbe_0 - 2.0 * g_xzzz_yyz_1[i] * fz_be_0 + g_xxzzz_yyz_1[i] * pa_x[i];

        g_xxxzzz_yzz_0[i] = 2.0 * g_xzzz_yzz_0[i] * fbe_0 - 2.0 * g_xzzz_yzz_1[i] * fz_be_0 + g_xxzzz_yzz_1[i] * pa_x[i];

        g_xxxzzz_zzz_0[i] = 2.0 * g_xzzz_zzz_0[i] * fbe_0 - 2.0 * g_xzzz_zzz_1[i] * fz_be_0 + g_xxzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : IF

    auto g_xxyyyy_xxx_0 = pbuffer.data(idx_eri_0_if + 100);

    auto g_xxyyyy_xxy_0 = pbuffer.data(idx_eri_0_if + 101);

    auto g_xxyyyy_xxz_0 = pbuffer.data(idx_eri_0_if + 102);

    auto g_xxyyyy_xyy_0 = pbuffer.data(idx_eri_0_if + 103);

    auto g_xxyyyy_xyz_0 = pbuffer.data(idx_eri_0_if + 104);

    auto g_xxyyyy_xzz_0 = pbuffer.data(idx_eri_0_if + 105);

    auto g_xxyyyy_yyy_0 = pbuffer.data(idx_eri_0_if + 106);

    auto g_xxyyyy_yyz_0 = pbuffer.data(idx_eri_0_if + 107);

    auto g_xxyyyy_yzz_0 = pbuffer.data(idx_eri_0_if + 108);

    auto g_xxyyyy_zzz_0 = pbuffer.data(idx_eri_0_if + 109);

    #pragma omp simd aligned(g_xxyy_xxx_0, g_xxyy_xxx_1, g_xxyy_xxz_0, g_xxyy_xxz_1, g_xxyy_xzz_0, g_xxyy_xzz_1, g_xxyyy_xxx_1, g_xxyyy_xxz_1, g_xxyyy_xzz_1, g_xxyyyy_xxx_0, g_xxyyyy_xxy_0, g_xxyyyy_xxz_0, g_xxyyyy_xyy_0, g_xxyyyy_xyz_0, g_xxyyyy_xzz_0, g_xxyyyy_yyy_0, g_xxyyyy_yyz_0, g_xxyyyy_yzz_0, g_xxyyyy_zzz_0, g_xyyyy_xxy_1, g_xyyyy_xy_1, g_xyyyy_xyy_1, g_xyyyy_xyz_1, g_xyyyy_yy_1, g_xyyyy_yyy_1, g_xyyyy_yyz_1, g_xyyyy_yz_1, g_xyyyy_yzz_1, g_xyyyy_zzz_1, g_yyyy_xxy_0, g_yyyy_xxy_1, g_yyyy_xyy_0, g_yyyy_xyy_1, g_yyyy_xyz_0, g_yyyy_xyz_1, g_yyyy_yyy_0, g_yyyy_yyy_1, g_yyyy_yyz_0, g_yyyy_yyz_1, g_yyyy_yzz_0, g_yyyy_yzz_1, g_yyyy_zzz_0, g_yyyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyy_xxx_0[i] = 3.0 * g_xxyy_xxx_0[i] * fbe_0 - 3.0 * g_xxyy_xxx_1[i] * fz_be_0 + g_xxyyy_xxx_1[i] * pa_y[i];

        g_xxyyyy_xxy_0[i] = g_yyyy_xxy_0[i] * fbe_0 - g_yyyy_xxy_1[i] * fz_be_0 + 2.0 * g_xyyyy_xy_1[i] * fe_0 + g_xyyyy_xxy_1[i] * pa_x[i];

        g_xxyyyy_xxz_0[i] = 3.0 * g_xxyy_xxz_0[i] * fbe_0 - 3.0 * g_xxyy_xxz_1[i] * fz_be_0 + g_xxyyy_xxz_1[i] * pa_y[i];

        g_xxyyyy_xyy_0[i] = g_yyyy_xyy_0[i] * fbe_0 - g_yyyy_xyy_1[i] * fz_be_0 + g_xyyyy_yy_1[i] * fe_0 + g_xyyyy_xyy_1[i] * pa_x[i];

        g_xxyyyy_xyz_0[i] = g_yyyy_xyz_0[i] * fbe_0 - g_yyyy_xyz_1[i] * fz_be_0 + g_xyyyy_yz_1[i] * fe_0 + g_xyyyy_xyz_1[i] * pa_x[i];

        g_xxyyyy_xzz_0[i] = 3.0 * g_xxyy_xzz_0[i] * fbe_0 - 3.0 * g_xxyy_xzz_1[i] * fz_be_0 + g_xxyyy_xzz_1[i] * pa_y[i];

        g_xxyyyy_yyy_0[i] = g_yyyy_yyy_0[i] * fbe_0 - g_yyyy_yyy_1[i] * fz_be_0 + g_xyyyy_yyy_1[i] * pa_x[i];

        g_xxyyyy_yyz_0[i] = g_yyyy_yyz_0[i] * fbe_0 - g_yyyy_yyz_1[i] * fz_be_0 + g_xyyyy_yyz_1[i] * pa_x[i];

        g_xxyyyy_yzz_0[i] = g_yyyy_yzz_0[i] * fbe_0 - g_yyyy_yzz_1[i] * fz_be_0 + g_xyyyy_yzz_1[i] * pa_x[i];

        g_xxyyyy_zzz_0[i] = g_yyyy_zzz_0[i] * fbe_0 - g_yyyy_zzz_1[i] * fz_be_0 + g_xyyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 110-120 components of targeted buffer : IF

    auto g_xxyyyz_xxx_0 = pbuffer.data(idx_eri_0_if + 110);

    auto g_xxyyyz_xxy_0 = pbuffer.data(idx_eri_0_if + 111);

    auto g_xxyyyz_xxz_0 = pbuffer.data(idx_eri_0_if + 112);

    auto g_xxyyyz_xyy_0 = pbuffer.data(idx_eri_0_if + 113);

    auto g_xxyyyz_xyz_0 = pbuffer.data(idx_eri_0_if + 114);

    auto g_xxyyyz_xzz_0 = pbuffer.data(idx_eri_0_if + 115);

    auto g_xxyyyz_yyy_0 = pbuffer.data(idx_eri_0_if + 116);

    auto g_xxyyyz_yyz_0 = pbuffer.data(idx_eri_0_if + 117);

    auto g_xxyyyz_yzz_0 = pbuffer.data(idx_eri_0_if + 118);

    auto g_xxyyyz_zzz_0 = pbuffer.data(idx_eri_0_if + 119);

    #pragma omp simd aligned(g_xxyyy_xx_1, g_xxyyy_xxx_1, g_xxyyy_xxy_1, g_xxyyy_xxz_1, g_xxyyy_xy_1, g_xxyyy_xyy_1, g_xxyyy_xyz_1, g_xxyyy_xz_1, g_xxyyy_xzz_1, g_xxyyy_yy_1, g_xxyyy_yyy_1, g_xxyyy_yyz_1, g_xxyyy_yz_1, g_xxyyy_yzz_1, g_xxyyy_zz_1, g_xxyyy_zzz_1, g_xxyyyz_xxx_0, g_xxyyyz_xxy_0, g_xxyyyz_xxz_0, g_xxyyyz_xyy_0, g_xxyyyz_xyz_0, g_xxyyyz_xzz_0, g_xxyyyz_yyy_0, g_xxyyyz_yyz_0, g_xxyyyz_yzz_0, g_xxyyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyz_xxx_0[i] = g_xxyyy_xxx_1[i] * pa_z[i];

        g_xxyyyz_xxy_0[i] = g_xxyyy_xxy_1[i] * pa_z[i];

        g_xxyyyz_xxz_0[i] = g_xxyyy_xx_1[i] * fe_0 + g_xxyyy_xxz_1[i] * pa_z[i];

        g_xxyyyz_xyy_0[i] = g_xxyyy_xyy_1[i] * pa_z[i];

        g_xxyyyz_xyz_0[i] = g_xxyyy_xy_1[i] * fe_0 + g_xxyyy_xyz_1[i] * pa_z[i];

        g_xxyyyz_xzz_0[i] = 2.0 * g_xxyyy_xz_1[i] * fe_0 + g_xxyyy_xzz_1[i] * pa_z[i];

        g_xxyyyz_yyy_0[i] = g_xxyyy_yyy_1[i] * pa_z[i];

        g_xxyyyz_yyz_0[i] = g_xxyyy_yy_1[i] * fe_0 + g_xxyyy_yyz_1[i] * pa_z[i];

        g_xxyyyz_yzz_0[i] = 2.0 * g_xxyyy_yz_1[i] * fe_0 + g_xxyyy_yzz_1[i] * pa_z[i];

        g_xxyyyz_zzz_0[i] = 3.0 * g_xxyyy_zz_1[i] * fe_0 + g_xxyyy_zzz_1[i] * pa_z[i];
    }

    // Set up 120-130 components of targeted buffer : IF

    auto g_xxyyzz_xxx_0 = pbuffer.data(idx_eri_0_if + 120);

    auto g_xxyyzz_xxy_0 = pbuffer.data(idx_eri_0_if + 121);

    auto g_xxyyzz_xxz_0 = pbuffer.data(idx_eri_0_if + 122);

    auto g_xxyyzz_xyy_0 = pbuffer.data(idx_eri_0_if + 123);

    auto g_xxyyzz_xyz_0 = pbuffer.data(idx_eri_0_if + 124);

    auto g_xxyyzz_xzz_0 = pbuffer.data(idx_eri_0_if + 125);

    auto g_xxyyzz_yyy_0 = pbuffer.data(idx_eri_0_if + 126);

    auto g_xxyyzz_yyz_0 = pbuffer.data(idx_eri_0_if + 127);

    auto g_xxyyzz_yzz_0 = pbuffer.data(idx_eri_0_if + 128);

    auto g_xxyyzz_zzz_0 = pbuffer.data(idx_eri_0_if + 129);

    #pragma omp simd aligned(g_xxyy_xxy_0, g_xxyy_xxy_1, g_xxyy_xyy_0, g_xxyy_xyy_1, g_xxyyz_xxy_1, g_xxyyz_xyy_1, g_xxyyzz_xxx_0, g_xxyyzz_xxy_0, g_xxyyzz_xxz_0, g_xxyyzz_xyy_0, g_xxyyzz_xyz_0, g_xxyyzz_xzz_0, g_xxyyzz_yyy_0, g_xxyyzz_yyz_0, g_xxyyzz_yzz_0, g_xxyyzz_zzz_0, g_xxyzz_xxx_1, g_xxyzz_xxz_1, g_xxyzz_xzz_1, g_xxzz_xxx_0, g_xxzz_xxx_1, g_xxzz_xxz_0, g_xxzz_xxz_1, g_xxzz_xzz_0, g_xxzz_xzz_1, g_xyyzz_xyz_1, g_xyyzz_yyy_1, g_xyyzz_yyz_1, g_xyyzz_yz_1, g_xyyzz_yzz_1, g_xyyzz_zzz_1, g_yyzz_xyz_0, g_yyzz_xyz_1, g_yyzz_yyy_0, g_yyzz_yyy_1, g_yyzz_yyz_0, g_yyzz_yyz_1, g_yyzz_yzz_0, g_yyzz_yzz_1, g_yyzz_zzz_0, g_yyzz_zzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzz_xxx_0[i] = g_xxzz_xxx_0[i] * fbe_0 - g_xxzz_xxx_1[i] * fz_be_0 + g_xxyzz_xxx_1[i] * pa_y[i];

        g_xxyyzz_xxy_0[i] = g_xxyy_xxy_0[i] * fbe_0 - g_xxyy_xxy_1[i] * fz_be_0 + g_xxyyz_xxy_1[i] * pa_z[i];

        g_xxyyzz_xxz_0[i] = g_xxzz_xxz_0[i] * fbe_0 - g_xxzz_xxz_1[i] * fz_be_0 + g_xxyzz_xxz_1[i] * pa_y[i];

        g_xxyyzz_xyy_0[i] = g_xxyy_xyy_0[i] * fbe_0 - g_xxyy_xyy_1[i] * fz_be_0 + g_xxyyz_xyy_1[i] * pa_z[i];

        g_xxyyzz_xyz_0[i] = g_yyzz_xyz_0[i] * fbe_0 - g_yyzz_xyz_1[i] * fz_be_0 + g_xyyzz_yz_1[i] * fe_0 + g_xyyzz_xyz_1[i] * pa_x[i];

        g_xxyyzz_xzz_0[i] = g_xxzz_xzz_0[i] * fbe_0 - g_xxzz_xzz_1[i] * fz_be_0 + g_xxyzz_xzz_1[i] * pa_y[i];

        g_xxyyzz_yyy_0[i] = g_yyzz_yyy_0[i] * fbe_0 - g_yyzz_yyy_1[i] * fz_be_0 + g_xyyzz_yyy_1[i] * pa_x[i];

        g_xxyyzz_yyz_0[i] = g_yyzz_yyz_0[i] * fbe_0 - g_yyzz_yyz_1[i] * fz_be_0 + g_xyyzz_yyz_1[i] * pa_x[i];

        g_xxyyzz_yzz_0[i] = g_yyzz_yzz_0[i] * fbe_0 - g_yyzz_yzz_1[i] * fz_be_0 + g_xyyzz_yzz_1[i] * pa_x[i];

        g_xxyyzz_zzz_0[i] = g_yyzz_zzz_0[i] * fbe_0 - g_yyzz_zzz_1[i] * fz_be_0 + g_xyyzz_zzz_1[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : IF

    auto g_xxyzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 130);

    auto g_xxyzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 131);

    auto g_xxyzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 132);

    auto g_xxyzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 133);

    auto g_xxyzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 134);

    auto g_xxyzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 135);

    auto g_xxyzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 136);

    auto g_xxyzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 137);

    auto g_xxyzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 138);

    auto g_xxyzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 139);

    #pragma omp simd aligned(g_xxyzzz_xxx_0, g_xxyzzz_xxy_0, g_xxyzzz_xxz_0, g_xxyzzz_xyy_0, g_xxyzzz_xyz_0, g_xxyzzz_xzz_0, g_xxyzzz_yyy_0, g_xxyzzz_yyz_0, g_xxyzzz_yzz_0, g_xxyzzz_zzz_0, g_xxzzz_xx_1, g_xxzzz_xxx_1, g_xxzzz_xxy_1, g_xxzzz_xxz_1, g_xxzzz_xy_1, g_xxzzz_xyy_1, g_xxzzz_xyz_1, g_xxzzz_xz_1, g_xxzzz_xzz_1, g_xxzzz_yy_1, g_xxzzz_yyy_1, g_xxzzz_yyz_1, g_xxzzz_yz_1, g_xxzzz_yzz_1, g_xxzzz_zz_1, g_xxzzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzz_xxx_0[i] = g_xxzzz_xxx_1[i] * pa_y[i];

        g_xxyzzz_xxy_0[i] = g_xxzzz_xx_1[i] * fe_0 + g_xxzzz_xxy_1[i] * pa_y[i];

        g_xxyzzz_xxz_0[i] = g_xxzzz_xxz_1[i] * pa_y[i];

        g_xxyzzz_xyy_0[i] = 2.0 * g_xxzzz_xy_1[i] * fe_0 + g_xxzzz_xyy_1[i] * pa_y[i];

        g_xxyzzz_xyz_0[i] = g_xxzzz_xz_1[i] * fe_0 + g_xxzzz_xyz_1[i] * pa_y[i];

        g_xxyzzz_xzz_0[i] = g_xxzzz_xzz_1[i] * pa_y[i];

        g_xxyzzz_yyy_0[i] = 3.0 * g_xxzzz_yy_1[i] * fe_0 + g_xxzzz_yyy_1[i] * pa_y[i];

        g_xxyzzz_yyz_0[i] = 2.0 * g_xxzzz_yz_1[i] * fe_0 + g_xxzzz_yyz_1[i] * pa_y[i];

        g_xxyzzz_yzz_0[i] = g_xxzzz_zz_1[i] * fe_0 + g_xxzzz_yzz_1[i] * pa_y[i];

        g_xxyzzz_zzz_0[i] = g_xxzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : IF

    auto g_xxzzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 140);

    auto g_xxzzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 141);

    auto g_xxzzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 142);

    auto g_xxzzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 143);

    auto g_xxzzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 144);

    auto g_xxzzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 145);

    auto g_xxzzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 146);

    auto g_xxzzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 147);

    auto g_xxzzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 148);

    auto g_xxzzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 149);

    #pragma omp simd aligned(g_xxzz_xxx_0, g_xxzz_xxx_1, g_xxzz_xxy_0, g_xxzz_xxy_1, g_xxzz_xyy_0, g_xxzz_xyy_1, g_xxzzz_xxx_1, g_xxzzz_xxy_1, g_xxzzz_xyy_1, g_xxzzzz_xxx_0, g_xxzzzz_xxy_0, g_xxzzzz_xxz_0, g_xxzzzz_xyy_0, g_xxzzzz_xyz_0, g_xxzzzz_xzz_0, g_xxzzzz_yyy_0, g_xxzzzz_yyz_0, g_xxzzzz_yzz_0, g_xxzzzz_zzz_0, g_xzzzz_xxz_1, g_xzzzz_xyz_1, g_xzzzz_xz_1, g_xzzzz_xzz_1, g_xzzzz_yyy_1, g_xzzzz_yyz_1, g_xzzzz_yz_1, g_xzzzz_yzz_1, g_xzzzz_zz_1, g_xzzzz_zzz_1, g_zzzz_xxz_0, g_zzzz_xxz_1, g_zzzz_xyz_0, g_zzzz_xyz_1, g_zzzz_xzz_0, g_zzzz_xzz_1, g_zzzz_yyy_0, g_zzzz_yyy_1, g_zzzz_yyz_0, g_zzzz_yyz_1, g_zzzz_yzz_0, g_zzzz_yzz_1, g_zzzz_zzz_0, g_zzzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzz_xxx_0[i] = 3.0 * g_xxzz_xxx_0[i] * fbe_0 - 3.0 * g_xxzz_xxx_1[i] * fz_be_0 + g_xxzzz_xxx_1[i] * pa_z[i];

        g_xxzzzz_xxy_0[i] = 3.0 * g_xxzz_xxy_0[i] * fbe_0 - 3.0 * g_xxzz_xxy_1[i] * fz_be_0 + g_xxzzz_xxy_1[i] * pa_z[i];

        g_xxzzzz_xxz_0[i] = g_zzzz_xxz_0[i] * fbe_0 - g_zzzz_xxz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xz_1[i] * fe_0 + g_xzzzz_xxz_1[i] * pa_x[i];

        g_xxzzzz_xyy_0[i] = 3.0 * g_xxzz_xyy_0[i] * fbe_0 - 3.0 * g_xxzz_xyy_1[i] * fz_be_0 + g_xxzzz_xyy_1[i] * pa_z[i];

        g_xxzzzz_xyz_0[i] = g_zzzz_xyz_0[i] * fbe_0 - g_zzzz_xyz_1[i] * fz_be_0 + g_xzzzz_yz_1[i] * fe_0 + g_xzzzz_xyz_1[i] * pa_x[i];

        g_xxzzzz_xzz_0[i] = g_zzzz_xzz_0[i] * fbe_0 - g_zzzz_xzz_1[i] * fz_be_0 + g_xzzzz_zz_1[i] * fe_0 + g_xzzzz_xzz_1[i] * pa_x[i];

        g_xxzzzz_yyy_0[i] = g_zzzz_yyy_0[i] * fbe_0 - g_zzzz_yyy_1[i] * fz_be_0 + g_xzzzz_yyy_1[i] * pa_x[i];

        g_xxzzzz_yyz_0[i] = g_zzzz_yyz_0[i] * fbe_0 - g_zzzz_yyz_1[i] * fz_be_0 + g_xzzzz_yyz_1[i] * pa_x[i];

        g_xxzzzz_yzz_0[i] = g_zzzz_yzz_0[i] * fbe_0 - g_zzzz_yzz_1[i] * fz_be_0 + g_xzzzz_yzz_1[i] * pa_x[i];

        g_xxzzzz_zzz_0[i] = g_zzzz_zzz_0[i] * fbe_0 - g_zzzz_zzz_1[i] * fz_be_0 + g_xzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : IF

    auto g_xyyyyy_xxx_0 = pbuffer.data(idx_eri_0_if + 150);

    auto g_xyyyyy_xxy_0 = pbuffer.data(idx_eri_0_if + 151);

    auto g_xyyyyy_xxz_0 = pbuffer.data(idx_eri_0_if + 152);

    auto g_xyyyyy_xyy_0 = pbuffer.data(idx_eri_0_if + 153);

    auto g_xyyyyy_xyz_0 = pbuffer.data(idx_eri_0_if + 154);

    auto g_xyyyyy_xzz_0 = pbuffer.data(idx_eri_0_if + 155);

    auto g_xyyyyy_yyy_0 = pbuffer.data(idx_eri_0_if + 156);

    auto g_xyyyyy_yyz_0 = pbuffer.data(idx_eri_0_if + 157);

    auto g_xyyyyy_yzz_0 = pbuffer.data(idx_eri_0_if + 158);

    auto g_xyyyyy_zzz_0 = pbuffer.data(idx_eri_0_if + 159);

    #pragma omp simd aligned(g_xyyyyy_xxx_0, g_xyyyyy_xxy_0, g_xyyyyy_xxz_0, g_xyyyyy_xyy_0, g_xyyyyy_xyz_0, g_xyyyyy_xzz_0, g_xyyyyy_yyy_0, g_xyyyyy_yyz_0, g_xyyyyy_yzz_0, g_xyyyyy_zzz_0, g_yyyyy_xx_1, g_yyyyy_xxx_1, g_yyyyy_xxy_1, g_yyyyy_xxz_1, g_yyyyy_xy_1, g_yyyyy_xyy_1, g_yyyyy_xyz_1, g_yyyyy_xz_1, g_yyyyy_xzz_1, g_yyyyy_yy_1, g_yyyyy_yyy_1, g_yyyyy_yyz_1, g_yyyyy_yz_1, g_yyyyy_yzz_1, g_yyyyy_zz_1, g_yyyyy_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyy_xxx_0[i] = 3.0 * g_yyyyy_xx_1[i] * fe_0 + g_yyyyy_xxx_1[i] * pa_x[i];

        g_xyyyyy_xxy_0[i] = 2.0 * g_yyyyy_xy_1[i] * fe_0 + g_yyyyy_xxy_1[i] * pa_x[i];

        g_xyyyyy_xxz_0[i] = 2.0 * g_yyyyy_xz_1[i] * fe_0 + g_yyyyy_xxz_1[i] * pa_x[i];

        g_xyyyyy_xyy_0[i] = g_yyyyy_yy_1[i] * fe_0 + g_yyyyy_xyy_1[i] * pa_x[i];

        g_xyyyyy_xyz_0[i] = g_yyyyy_yz_1[i] * fe_0 + g_yyyyy_xyz_1[i] * pa_x[i];

        g_xyyyyy_xzz_0[i] = g_yyyyy_zz_1[i] * fe_0 + g_yyyyy_xzz_1[i] * pa_x[i];

        g_xyyyyy_yyy_0[i] = g_yyyyy_yyy_1[i] * pa_x[i];

        g_xyyyyy_yyz_0[i] = g_yyyyy_yyz_1[i] * pa_x[i];

        g_xyyyyy_yzz_0[i] = g_yyyyy_yzz_1[i] * pa_x[i];

        g_xyyyyy_zzz_0[i] = g_yyyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 160-170 components of targeted buffer : IF

    auto g_xyyyyz_xxx_0 = pbuffer.data(idx_eri_0_if + 160);

    auto g_xyyyyz_xxy_0 = pbuffer.data(idx_eri_0_if + 161);

    auto g_xyyyyz_xxz_0 = pbuffer.data(idx_eri_0_if + 162);

    auto g_xyyyyz_xyy_0 = pbuffer.data(idx_eri_0_if + 163);

    auto g_xyyyyz_xyz_0 = pbuffer.data(idx_eri_0_if + 164);

    auto g_xyyyyz_xzz_0 = pbuffer.data(idx_eri_0_if + 165);

    auto g_xyyyyz_yyy_0 = pbuffer.data(idx_eri_0_if + 166);

    auto g_xyyyyz_yyz_0 = pbuffer.data(idx_eri_0_if + 167);

    auto g_xyyyyz_yzz_0 = pbuffer.data(idx_eri_0_if + 168);

    auto g_xyyyyz_zzz_0 = pbuffer.data(idx_eri_0_if + 169);

    #pragma omp simd aligned(g_xyyyy_xxx_1, g_xyyyy_xxy_1, g_xyyyy_xyy_1, g_xyyyyz_xxx_0, g_xyyyyz_xxy_0, g_xyyyyz_xxz_0, g_xyyyyz_xyy_0, g_xyyyyz_xyz_0, g_xyyyyz_xzz_0, g_xyyyyz_yyy_0, g_xyyyyz_yyz_0, g_xyyyyz_yzz_0, g_xyyyyz_zzz_0, g_yyyyz_xxz_1, g_yyyyz_xyz_1, g_yyyyz_xz_1, g_yyyyz_xzz_1, g_yyyyz_yyy_1, g_yyyyz_yyz_1, g_yyyyz_yz_1, g_yyyyz_yzz_1, g_yyyyz_zz_1, g_yyyyz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyz_xxx_0[i] = g_xyyyy_xxx_1[i] * pa_z[i];

        g_xyyyyz_xxy_0[i] = g_xyyyy_xxy_1[i] * pa_z[i];

        g_xyyyyz_xxz_0[i] = 2.0 * g_yyyyz_xz_1[i] * fe_0 + g_yyyyz_xxz_1[i] * pa_x[i];

        g_xyyyyz_xyy_0[i] = g_xyyyy_xyy_1[i] * pa_z[i];

        g_xyyyyz_xyz_0[i] = g_yyyyz_yz_1[i] * fe_0 + g_yyyyz_xyz_1[i] * pa_x[i];

        g_xyyyyz_xzz_0[i] = g_yyyyz_zz_1[i] * fe_0 + g_yyyyz_xzz_1[i] * pa_x[i];

        g_xyyyyz_yyy_0[i] = g_yyyyz_yyy_1[i] * pa_x[i];

        g_xyyyyz_yyz_0[i] = g_yyyyz_yyz_1[i] * pa_x[i];

        g_xyyyyz_yzz_0[i] = g_yyyyz_yzz_1[i] * pa_x[i];

        g_xyyyyz_zzz_0[i] = g_yyyyz_zzz_1[i] * pa_x[i];
    }

    // Set up 170-180 components of targeted buffer : IF

    auto g_xyyyzz_xxx_0 = pbuffer.data(idx_eri_0_if + 170);

    auto g_xyyyzz_xxy_0 = pbuffer.data(idx_eri_0_if + 171);

    auto g_xyyyzz_xxz_0 = pbuffer.data(idx_eri_0_if + 172);

    auto g_xyyyzz_xyy_0 = pbuffer.data(idx_eri_0_if + 173);

    auto g_xyyyzz_xyz_0 = pbuffer.data(idx_eri_0_if + 174);

    auto g_xyyyzz_xzz_0 = pbuffer.data(idx_eri_0_if + 175);

    auto g_xyyyzz_yyy_0 = pbuffer.data(idx_eri_0_if + 176);

    auto g_xyyyzz_yyz_0 = pbuffer.data(idx_eri_0_if + 177);

    auto g_xyyyzz_yzz_0 = pbuffer.data(idx_eri_0_if + 178);

    auto g_xyyyzz_zzz_0 = pbuffer.data(idx_eri_0_if + 179);

    #pragma omp simd aligned(g_xyyyzz_xxx_0, g_xyyyzz_xxy_0, g_xyyyzz_xxz_0, g_xyyyzz_xyy_0, g_xyyyzz_xyz_0, g_xyyyzz_xzz_0, g_xyyyzz_yyy_0, g_xyyyzz_yyz_0, g_xyyyzz_yzz_0, g_xyyyzz_zzz_0, g_yyyzz_xx_1, g_yyyzz_xxx_1, g_yyyzz_xxy_1, g_yyyzz_xxz_1, g_yyyzz_xy_1, g_yyyzz_xyy_1, g_yyyzz_xyz_1, g_yyyzz_xz_1, g_yyyzz_xzz_1, g_yyyzz_yy_1, g_yyyzz_yyy_1, g_yyyzz_yyz_1, g_yyyzz_yz_1, g_yyyzz_yzz_1, g_yyyzz_zz_1, g_yyyzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzz_xxx_0[i] = 3.0 * g_yyyzz_xx_1[i] * fe_0 + g_yyyzz_xxx_1[i] * pa_x[i];

        g_xyyyzz_xxy_0[i] = 2.0 * g_yyyzz_xy_1[i] * fe_0 + g_yyyzz_xxy_1[i] * pa_x[i];

        g_xyyyzz_xxz_0[i] = 2.0 * g_yyyzz_xz_1[i] * fe_0 + g_yyyzz_xxz_1[i] * pa_x[i];

        g_xyyyzz_xyy_0[i] = g_yyyzz_yy_1[i] * fe_0 + g_yyyzz_xyy_1[i] * pa_x[i];

        g_xyyyzz_xyz_0[i] = g_yyyzz_yz_1[i] * fe_0 + g_yyyzz_xyz_1[i] * pa_x[i];

        g_xyyyzz_xzz_0[i] = g_yyyzz_zz_1[i] * fe_0 + g_yyyzz_xzz_1[i] * pa_x[i];

        g_xyyyzz_yyy_0[i] = g_yyyzz_yyy_1[i] * pa_x[i];

        g_xyyyzz_yyz_0[i] = g_yyyzz_yyz_1[i] * pa_x[i];

        g_xyyyzz_yzz_0[i] = g_yyyzz_yzz_1[i] * pa_x[i];

        g_xyyyzz_zzz_0[i] = g_yyyzz_zzz_1[i] * pa_x[i];
    }

    // Set up 180-190 components of targeted buffer : IF

    auto g_xyyzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 180);

    auto g_xyyzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 181);

    auto g_xyyzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 182);

    auto g_xyyzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 183);

    auto g_xyyzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 184);

    auto g_xyyzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 185);

    auto g_xyyzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 186);

    auto g_xyyzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 187);

    auto g_xyyzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 188);

    auto g_xyyzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 189);

    #pragma omp simd aligned(g_xyyzzz_xxx_0, g_xyyzzz_xxy_0, g_xyyzzz_xxz_0, g_xyyzzz_xyy_0, g_xyyzzz_xyz_0, g_xyyzzz_xzz_0, g_xyyzzz_yyy_0, g_xyyzzz_yyz_0, g_xyyzzz_yzz_0, g_xyyzzz_zzz_0, g_yyzzz_xx_1, g_yyzzz_xxx_1, g_yyzzz_xxy_1, g_yyzzz_xxz_1, g_yyzzz_xy_1, g_yyzzz_xyy_1, g_yyzzz_xyz_1, g_yyzzz_xz_1, g_yyzzz_xzz_1, g_yyzzz_yy_1, g_yyzzz_yyy_1, g_yyzzz_yyz_1, g_yyzzz_yz_1, g_yyzzz_yzz_1, g_yyzzz_zz_1, g_yyzzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzz_xxx_0[i] = 3.0 * g_yyzzz_xx_1[i] * fe_0 + g_yyzzz_xxx_1[i] * pa_x[i];

        g_xyyzzz_xxy_0[i] = 2.0 * g_yyzzz_xy_1[i] * fe_0 + g_yyzzz_xxy_1[i] * pa_x[i];

        g_xyyzzz_xxz_0[i] = 2.0 * g_yyzzz_xz_1[i] * fe_0 + g_yyzzz_xxz_1[i] * pa_x[i];

        g_xyyzzz_xyy_0[i] = g_yyzzz_yy_1[i] * fe_0 + g_yyzzz_xyy_1[i] * pa_x[i];

        g_xyyzzz_xyz_0[i] = g_yyzzz_yz_1[i] * fe_0 + g_yyzzz_xyz_1[i] * pa_x[i];

        g_xyyzzz_xzz_0[i] = g_yyzzz_zz_1[i] * fe_0 + g_yyzzz_xzz_1[i] * pa_x[i];

        g_xyyzzz_yyy_0[i] = g_yyzzz_yyy_1[i] * pa_x[i];

        g_xyyzzz_yyz_0[i] = g_yyzzz_yyz_1[i] * pa_x[i];

        g_xyyzzz_yzz_0[i] = g_yyzzz_yzz_1[i] * pa_x[i];

        g_xyyzzz_zzz_0[i] = g_yyzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 190-200 components of targeted buffer : IF

    auto g_xyzzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 190);

    auto g_xyzzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 191);

    auto g_xyzzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 192);

    auto g_xyzzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 193);

    auto g_xyzzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 194);

    auto g_xyzzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 195);

    auto g_xyzzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 196);

    auto g_xyzzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 197);

    auto g_xyzzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 198);

    auto g_xyzzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 199);

    #pragma omp simd aligned(g_xyzzzz_xxx_0, g_xyzzzz_xxy_0, g_xyzzzz_xxz_0, g_xyzzzz_xyy_0, g_xyzzzz_xyz_0, g_xyzzzz_xzz_0, g_xyzzzz_yyy_0, g_xyzzzz_yyz_0, g_xyzzzz_yzz_0, g_xyzzzz_zzz_0, g_xzzzz_xxx_1, g_xzzzz_xxz_1, g_xzzzz_xzz_1, g_yzzzz_xxy_1, g_yzzzz_xy_1, g_yzzzz_xyy_1, g_yzzzz_xyz_1, g_yzzzz_yy_1, g_yzzzz_yyy_1, g_yzzzz_yyz_1, g_yzzzz_yz_1, g_yzzzz_yzz_1, g_yzzzz_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzz_xxx_0[i] = g_xzzzz_xxx_1[i] * pa_y[i];

        g_xyzzzz_xxy_0[i] = 2.0 * g_yzzzz_xy_1[i] * fe_0 + g_yzzzz_xxy_1[i] * pa_x[i];

        g_xyzzzz_xxz_0[i] = g_xzzzz_xxz_1[i] * pa_y[i];

        g_xyzzzz_xyy_0[i] = g_yzzzz_yy_1[i] * fe_0 + g_yzzzz_xyy_1[i] * pa_x[i];

        g_xyzzzz_xyz_0[i] = g_yzzzz_yz_1[i] * fe_0 + g_yzzzz_xyz_1[i] * pa_x[i];

        g_xyzzzz_xzz_0[i] = g_xzzzz_xzz_1[i] * pa_y[i];

        g_xyzzzz_yyy_0[i] = g_yzzzz_yyy_1[i] * pa_x[i];

        g_xyzzzz_yyz_0[i] = g_yzzzz_yyz_1[i] * pa_x[i];

        g_xyzzzz_yzz_0[i] = g_yzzzz_yzz_1[i] * pa_x[i];

        g_xyzzzz_zzz_0[i] = g_yzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 200-210 components of targeted buffer : IF

    auto g_xzzzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 200);

    auto g_xzzzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 201);

    auto g_xzzzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 202);

    auto g_xzzzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 203);

    auto g_xzzzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 204);

    auto g_xzzzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 205);

    auto g_xzzzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 206);

    auto g_xzzzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 207);

    auto g_xzzzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 208);

    auto g_xzzzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 209);

    #pragma omp simd aligned(g_xzzzzz_xxx_0, g_xzzzzz_xxy_0, g_xzzzzz_xxz_0, g_xzzzzz_xyy_0, g_xzzzzz_xyz_0, g_xzzzzz_xzz_0, g_xzzzzz_yyy_0, g_xzzzzz_yyz_0, g_xzzzzz_yzz_0, g_xzzzzz_zzz_0, g_zzzzz_xx_1, g_zzzzz_xxx_1, g_zzzzz_xxy_1, g_zzzzz_xxz_1, g_zzzzz_xy_1, g_zzzzz_xyy_1, g_zzzzz_xyz_1, g_zzzzz_xz_1, g_zzzzz_xzz_1, g_zzzzz_yy_1, g_zzzzz_yyy_1, g_zzzzz_yyz_1, g_zzzzz_yz_1, g_zzzzz_yzz_1, g_zzzzz_zz_1, g_zzzzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzz_xxx_0[i] = 3.0 * g_zzzzz_xx_1[i] * fe_0 + g_zzzzz_xxx_1[i] * pa_x[i];

        g_xzzzzz_xxy_0[i] = 2.0 * g_zzzzz_xy_1[i] * fe_0 + g_zzzzz_xxy_1[i] * pa_x[i];

        g_xzzzzz_xxz_0[i] = 2.0 * g_zzzzz_xz_1[i] * fe_0 + g_zzzzz_xxz_1[i] * pa_x[i];

        g_xzzzzz_xyy_0[i] = g_zzzzz_yy_1[i] * fe_0 + g_zzzzz_xyy_1[i] * pa_x[i];

        g_xzzzzz_xyz_0[i] = g_zzzzz_yz_1[i] * fe_0 + g_zzzzz_xyz_1[i] * pa_x[i];

        g_xzzzzz_xzz_0[i] = g_zzzzz_zz_1[i] * fe_0 + g_zzzzz_xzz_1[i] * pa_x[i];

        g_xzzzzz_yyy_0[i] = g_zzzzz_yyy_1[i] * pa_x[i];

        g_xzzzzz_yyz_0[i] = g_zzzzz_yyz_1[i] * pa_x[i];

        g_xzzzzz_yzz_0[i] = g_zzzzz_yzz_1[i] * pa_x[i];

        g_xzzzzz_zzz_0[i] = g_zzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 210-220 components of targeted buffer : IF

    auto g_yyyyyy_xxx_0 = pbuffer.data(idx_eri_0_if + 210);

    auto g_yyyyyy_xxy_0 = pbuffer.data(idx_eri_0_if + 211);

    auto g_yyyyyy_xxz_0 = pbuffer.data(idx_eri_0_if + 212);

    auto g_yyyyyy_xyy_0 = pbuffer.data(idx_eri_0_if + 213);

    auto g_yyyyyy_xyz_0 = pbuffer.data(idx_eri_0_if + 214);

    auto g_yyyyyy_xzz_0 = pbuffer.data(idx_eri_0_if + 215);

    auto g_yyyyyy_yyy_0 = pbuffer.data(idx_eri_0_if + 216);

    auto g_yyyyyy_yyz_0 = pbuffer.data(idx_eri_0_if + 217);

    auto g_yyyyyy_yzz_0 = pbuffer.data(idx_eri_0_if + 218);

    auto g_yyyyyy_zzz_0 = pbuffer.data(idx_eri_0_if + 219);

    #pragma omp simd aligned(g_yyyy_xxx_0, g_yyyy_xxx_1, g_yyyy_xxy_0, g_yyyy_xxy_1, g_yyyy_xxz_0, g_yyyy_xxz_1, g_yyyy_xyy_0, g_yyyy_xyy_1, g_yyyy_xyz_0, g_yyyy_xyz_1, g_yyyy_xzz_0, g_yyyy_xzz_1, g_yyyy_yyy_0, g_yyyy_yyy_1, g_yyyy_yyz_0, g_yyyy_yyz_1, g_yyyy_yzz_0, g_yyyy_yzz_1, g_yyyy_zzz_0, g_yyyy_zzz_1, g_yyyyy_xx_1, g_yyyyy_xxx_1, g_yyyyy_xxy_1, g_yyyyy_xxz_1, g_yyyyy_xy_1, g_yyyyy_xyy_1, g_yyyyy_xyz_1, g_yyyyy_xz_1, g_yyyyy_xzz_1, g_yyyyy_yy_1, g_yyyyy_yyy_1, g_yyyyy_yyz_1, g_yyyyy_yz_1, g_yyyyy_yzz_1, g_yyyyy_zz_1, g_yyyyy_zzz_1, g_yyyyyy_xxx_0, g_yyyyyy_xxy_0, g_yyyyyy_xxz_0, g_yyyyyy_xyy_0, g_yyyyyy_xyz_0, g_yyyyyy_xzz_0, g_yyyyyy_yyy_0, g_yyyyyy_yyz_0, g_yyyyyy_yzz_0, g_yyyyyy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyy_xxx_0[i] = 5.0 * g_yyyy_xxx_0[i] * fbe_0 - 5.0 * g_yyyy_xxx_1[i] * fz_be_0 + g_yyyyy_xxx_1[i] * pa_y[i];

        g_yyyyyy_xxy_0[i] = 5.0 * g_yyyy_xxy_0[i] * fbe_0 - 5.0 * g_yyyy_xxy_1[i] * fz_be_0 + g_yyyyy_xx_1[i] * fe_0 + g_yyyyy_xxy_1[i] * pa_y[i];

        g_yyyyyy_xxz_0[i] = 5.0 * g_yyyy_xxz_0[i] * fbe_0 - 5.0 * g_yyyy_xxz_1[i] * fz_be_0 + g_yyyyy_xxz_1[i] * pa_y[i];

        g_yyyyyy_xyy_0[i] = 5.0 * g_yyyy_xyy_0[i] * fbe_0 - 5.0 * g_yyyy_xyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_xy_1[i] * fe_0 + g_yyyyy_xyy_1[i] * pa_y[i];

        g_yyyyyy_xyz_0[i] = 5.0 * g_yyyy_xyz_0[i] * fbe_0 - 5.0 * g_yyyy_xyz_1[i] * fz_be_0 + g_yyyyy_xz_1[i] * fe_0 + g_yyyyy_xyz_1[i] * pa_y[i];

        g_yyyyyy_xzz_0[i] = 5.0 * g_yyyy_xzz_0[i] * fbe_0 - 5.0 * g_yyyy_xzz_1[i] * fz_be_0 + g_yyyyy_xzz_1[i] * pa_y[i];

        g_yyyyyy_yyy_0[i] = 5.0 * g_yyyy_yyy_0[i] * fbe_0 - 5.0 * g_yyyy_yyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_yy_1[i] * fe_0 + g_yyyyy_yyy_1[i] * pa_y[i];

        g_yyyyyy_yyz_0[i] = 5.0 * g_yyyy_yyz_0[i] * fbe_0 - 5.0 * g_yyyy_yyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_yz_1[i] * fe_0 + g_yyyyy_yyz_1[i] * pa_y[i];

        g_yyyyyy_yzz_0[i] = 5.0 * g_yyyy_yzz_0[i] * fbe_0 - 5.0 * g_yyyy_yzz_1[i] * fz_be_0 + g_yyyyy_zz_1[i] * fe_0 + g_yyyyy_yzz_1[i] * pa_y[i];

        g_yyyyyy_zzz_0[i] = 5.0 * g_yyyy_zzz_0[i] * fbe_0 - 5.0 * g_yyyy_zzz_1[i] * fz_be_0 + g_yyyyy_zzz_1[i] * pa_y[i];
    }

    // Set up 220-230 components of targeted buffer : IF

    auto g_yyyyyz_xxx_0 = pbuffer.data(idx_eri_0_if + 220);

    auto g_yyyyyz_xxy_0 = pbuffer.data(idx_eri_0_if + 221);

    auto g_yyyyyz_xxz_0 = pbuffer.data(idx_eri_0_if + 222);

    auto g_yyyyyz_xyy_0 = pbuffer.data(idx_eri_0_if + 223);

    auto g_yyyyyz_xyz_0 = pbuffer.data(idx_eri_0_if + 224);

    auto g_yyyyyz_xzz_0 = pbuffer.data(idx_eri_0_if + 225);

    auto g_yyyyyz_yyy_0 = pbuffer.data(idx_eri_0_if + 226);

    auto g_yyyyyz_yyz_0 = pbuffer.data(idx_eri_0_if + 227);

    auto g_yyyyyz_yzz_0 = pbuffer.data(idx_eri_0_if + 228);

    auto g_yyyyyz_zzz_0 = pbuffer.data(idx_eri_0_if + 229);

    #pragma omp simd aligned(g_yyyyy_xx_1, g_yyyyy_xxx_1, g_yyyyy_xxy_1, g_yyyyy_xxz_1, g_yyyyy_xy_1, g_yyyyy_xyy_1, g_yyyyy_xyz_1, g_yyyyy_xz_1, g_yyyyy_xzz_1, g_yyyyy_yy_1, g_yyyyy_yyy_1, g_yyyyy_yyz_1, g_yyyyy_yz_1, g_yyyyy_yzz_1, g_yyyyy_zz_1, g_yyyyy_zzz_1, g_yyyyyz_xxx_0, g_yyyyyz_xxy_0, g_yyyyyz_xxz_0, g_yyyyyz_xyy_0, g_yyyyyz_xyz_0, g_yyyyyz_xzz_0, g_yyyyyz_yyy_0, g_yyyyyz_yyz_0, g_yyyyyz_yzz_0, g_yyyyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyz_xxx_0[i] = g_yyyyy_xxx_1[i] * pa_z[i];

        g_yyyyyz_xxy_0[i] = g_yyyyy_xxy_1[i] * pa_z[i];

        g_yyyyyz_xxz_0[i] = g_yyyyy_xx_1[i] * fe_0 + g_yyyyy_xxz_1[i] * pa_z[i];

        g_yyyyyz_xyy_0[i] = g_yyyyy_xyy_1[i] * pa_z[i];

        g_yyyyyz_xyz_0[i] = g_yyyyy_xy_1[i] * fe_0 + g_yyyyy_xyz_1[i] * pa_z[i];

        g_yyyyyz_xzz_0[i] = 2.0 * g_yyyyy_xz_1[i] * fe_0 + g_yyyyy_xzz_1[i] * pa_z[i];

        g_yyyyyz_yyy_0[i] = g_yyyyy_yyy_1[i] * pa_z[i];

        g_yyyyyz_yyz_0[i] = g_yyyyy_yy_1[i] * fe_0 + g_yyyyy_yyz_1[i] * pa_z[i];

        g_yyyyyz_yzz_0[i] = 2.0 * g_yyyyy_yz_1[i] * fe_0 + g_yyyyy_yzz_1[i] * pa_z[i];

        g_yyyyyz_zzz_0[i] = 3.0 * g_yyyyy_zz_1[i] * fe_0 + g_yyyyy_zzz_1[i] * pa_z[i];
    }

    // Set up 230-240 components of targeted buffer : IF

    auto g_yyyyzz_xxx_0 = pbuffer.data(idx_eri_0_if + 230);

    auto g_yyyyzz_xxy_0 = pbuffer.data(idx_eri_0_if + 231);

    auto g_yyyyzz_xxz_0 = pbuffer.data(idx_eri_0_if + 232);

    auto g_yyyyzz_xyy_0 = pbuffer.data(idx_eri_0_if + 233);

    auto g_yyyyzz_xyz_0 = pbuffer.data(idx_eri_0_if + 234);

    auto g_yyyyzz_xzz_0 = pbuffer.data(idx_eri_0_if + 235);

    auto g_yyyyzz_yyy_0 = pbuffer.data(idx_eri_0_if + 236);

    auto g_yyyyzz_yyz_0 = pbuffer.data(idx_eri_0_if + 237);

    auto g_yyyyzz_yzz_0 = pbuffer.data(idx_eri_0_if + 238);

    auto g_yyyyzz_zzz_0 = pbuffer.data(idx_eri_0_if + 239);

    #pragma omp simd aligned(g_yyyy_xxy_0, g_yyyy_xxy_1, g_yyyy_xyy_0, g_yyyy_xyy_1, g_yyyy_yyy_0, g_yyyy_yyy_1, g_yyyyz_xxy_1, g_yyyyz_xyy_1, g_yyyyz_yyy_1, g_yyyyzz_xxx_0, g_yyyyzz_xxy_0, g_yyyyzz_xxz_0, g_yyyyzz_xyy_0, g_yyyyzz_xyz_0, g_yyyyzz_xzz_0, g_yyyyzz_yyy_0, g_yyyyzz_yyz_0, g_yyyyzz_yzz_0, g_yyyyzz_zzz_0, g_yyyzz_xxx_1, g_yyyzz_xxz_1, g_yyyzz_xyz_1, g_yyyzz_xz_1, g_yyyzz_xzz_1, g_yyyzz_yyz_1, g_yyyzz_yz_1, g_yyyzz_yzz_1, g_yyyzz_zz_1, g_yyyzz_zzz_1, g_yyzz_xxx_0, g_yyzz_xxx_1, g_yyzz_xxz_0, g_yyzz_xxz_1, g_yyzz_xyz_0, g_yyzz_xyz_1, g_yyzz_xzz_0, g_yyzz_xzz_1, g_yyzz_yyz_0, g_yyzz_yyz_1, g_yyzz_yzz_0, g_yyzz_yzz_1, g_yyzz_zzz_0, g_yyzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzz_xxx_0[i] = 3.0 * g_yyzz_xxx_0[i] * fbe_0 - 3.0 * g_yyzz_xxx_1[i] * fz_be_0 + g_yyyzz_xxx_1[i] * pa_y[i];

        g_yyyyzz_xxy_0[i] = g_yyyy_xxy_0[i] * fbe_0 - g_yyyy_xxy_1[i] * fz_be_0 + g_yyyyz_xxy_1[i] * pa_z[i];

        g_yyyyzz_xxz_0[i] = 3.0 * g_yyzz_xxz_0[i] * fbe_0 - 3.0 * g_yyzz_xxz_1[i] * fz_be_0 + g_yyyzz_xxz_1[i] * pa_y[i];

        g_yyyyzz_xyy_0[i] = g_yyyy_xyy_0[i] * fbe_0 - g_yyyy_xyy_1[i] * fz_be_0 + g_yyyyz_xyy_1[i] * pa_z[i];

        g_yyyyzz_xyz_0[i] = 3.0 * g_yyzz_xyz_0[i] * fbe_0 - 3.0 * g_yyzz_xyz_1[i] * fz_be_0 + g_yyyzz_xz_1[i] * fe_0 + g_yyyzz_xyz_1[i] * pa_y[i];

        g_yyyyzz_xzz_0[i] = 3.0 * g_yyzz_xzz_0[i] * fbe_0 - 3.0 * g_yyzz_xzz_1[i] * fz_be_0 + g_yyyzz_xzz_1[i] * pa_y[i];

        g_yyyyzz_yyy_0[i] = g_yyyy_yyy_0[i] * fbe_0 - g_yyyy_yyy_1[i] * fz_be_0 + g_yyyyz_yyy_1[i] * pa_z[i];

        g_yyyyzz_yyz_0[i] = 3.0 * g_yyzz_yyz_0[i] * fbe_0 - 3.0 * g_yyzz_yyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_yz_1[i] * fe_0 + g_yyyzz_yyz_1[i] * pa_y[i];

        g_yyyyzz_yzz_0[i] = 3.0 * g_yyzz_yzz_0[i] * fbe_0 - 3.0 * g_yyzz_yzz_1[i] * fz_be_0 + g_yyyzz_zz_1[i] * fe_0 + g_yyyzz_yzz_1[i] * pa_y[i];

        g_yyyyzz_zzz_0[i] = 3.0 * g_yyzz_zzz_0[i] * fbe_0 - 3.0 * g_yyzz_zzz_1[i] * fz_be_0 + g_yyyzz_zzz_1[i] * pa_y[i];
    }

    // Set up 240-250 components of targeted buffer : IF

    auto g_yyyzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 240);

    auto g_yyyzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 241);

    auto g_yyyzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 242);

    auto g_yyyzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 243);

    auto g_yyyzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 244);

    auto g_yyyzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 245);

    auto g_yyyzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 246);

    auto g_yyyzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 247);

    auto g_yyyzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 248);

    auto g_yyyzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 249);

    #pragma omp simd aligned(g_yyyz_xxy_0, g_yyyz_xxy_1, g_yyyz_xyy_0, g_yyyz_xyy_1, g_yyyz_yyy_0, g_yyyz_yyy_1, g_yyyzz_xxy_1, g_yyyzz_xyy_1, g_yyyzz_yyy_1, g_yyyzzz_xxx_0, g_yyyzzz_xxy_0, g_yyyzzz_xxz_0, g_yyyzzz_xyy_0, g_yyyzzz_xyz_0, g_yyyzzz_xzz_0, g_yyyzzz_yyy_0, g_yyyzzz_yyz_0, g_yyyzzz_yzz_0, g_yyyzzz_zzz_0, g_yyzzz_xxx_1, g_yyzzz_xxz_1, g_yyzzz_xyz_1, g_yyzzz_xz_1, g_yyzzz_xzz_1, g_yyzzz_yyz_1, g_yyzzz_yz_1, g_yyzzz_yzz_1, g_yyzzz_zz_1, g_yyzzz_zzz_1, g_yzzz_xxx_0, g_yzzz_xxx_1, g_yzzz_xxz_0, g_yzzz_xxz_1, g_yzzz_xyz_0, g_yzzz_xyz_1, g_yzzz_xzz_0, g_yzzz_xzz_1, g_yzzz_yyz_0, g_yzzz_yyz_1, g_yzzz_yzz_0, g_yzzz_yzz_1, g_yzzz_zzz_0, g_yzzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzz_xxx_0[i] = 2.0 * g_yzzz_xxx_0[i] * fbe_0 - 2.0 * g_yzzz_xxx_1[i] * fz_be_0 + g_yyzzz_xxx_1[i] * pa_y[i];

        g_yyyzzz_xxy_0[i] = 2.0 * g_yyyz_xxy_0[i] * fbe_0 - 2.0 * g_yyyz_xxy_1[i] * fz_be_0 + g_yyyzz_xxy_1[i] * pa_z[i];

        g_yyyzzz_xxz_0[i] = 2.0 * g_yzzz_xxz_0[i] * fbe_0 - 2.0 * g_yzzz_xxz_1[i] * fz_be_0 + g_yyzzz_xxz_1[i] * pa_y[i];

        g_yyyzzz_xyy_0[i] = 2.0 * g_yyyz_xyy_0[i] * fbe_0 - 2.0 * g_yyyz_xyy_1[i] * fz_be_0 + g_yyyzz_xyy_1[i] * pa_z[i];

        g_yyyzzz_xyz_0[i] = 2.0 * g_yzzz_xyz_0[i] * fbe_0 - 2.0 * g_yzzz_xyz_1[i] * fz_be_0 + g_yyzzz_xz_1[i] * fe_0 + g_yyzzz_xyz_1[i] * pa_y[i];

        g_yyyzzz_xzz_0[i] = 2.0 * g_yzzz_xzz_0[i] * fbe_0 - 2.0 * g_yzzz_xzz_1[i] * fz_be_0 + g_yyzzz_xzz_1[i] * pa_y[i];

        g_yyyzzz_yyy_0[i] = 2.0 * g_yyyz_yyy_0[i] * fbe_0 - 2.0 * g_yyyz_yyy_1[i] * fz_be_0 + g_yyyzz_yyy_1[i] * pa_z[i];

        g_yyyzzz_yyz_0[i] = 2.0 * g_yzzz_yyz_0[i] * fbe_0 - 2.0 * g_yzzz_yyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_yz_1[i] * fe_0 + g_yyzzz_yyz_1[i] * pa_y[i];

        g_yyyzzz_yzz_0[i] = 2.0 * g_yzzz_yzz_0[i] * fbe_0 - 2.0 * g_yzzz_yzz_1[i] * fz_be_0 + g_yyzzz_zz_1[i] * fe_0 + g_yyzzz_yzz_1[i] * pa_y[i];

        g_yyyzzz_zzz_0[i] = 2.0 * g_yzzz_zzz_0[i] * fbe_0 - 2.0 * g_yzzz_zzz_1[i] * fz_be_0 + g_yyzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 250-260 components of targeted buffer : IF

    auto g_yyzzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 250);

    auto g_yyzzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 251);

    auto g_yyzzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 252);

    auto g_yyzzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 253);

    auto g_yyzzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 254);

    auto g_yyzzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 255);

    auto g_yyzzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 256);

    auto g_yyzzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 257);

    auto g_yyzzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 258);

    auto g_yyzzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 259);

    #pragma omp simd aligned(g_yyzz_xxy_0, g_yyzz_xxy_1, g_yyzz_xyy_0, g_yyzz_xyy_1, g_yyzz_yyy_0, g_yyzz_yyy_1, g_yyzzz_xxy_1, g_yyzzz_xyy_1, g_yyzzz_yyy_1, g_yyzzzz_xxx_0, g_yyzzzz_xxy_0, g_yyzzzz_xxz_0, g_yyzzzz_xyy_0, g_yyzzzz_xyz_0, g_yyzzzz_xzz_0, g_yyzzzz_yyy_0, g_yyzzzz_yyz_0, g_yyzzzz_yzz_0, g_yyzzzz_zzz_0, g_yzzzz_xxx_1, g_yzzzz_xxz_1, g_yzzzz_xyz_1, g_yzzzz_xz_1, g_yzzzz_xzz_1, g_yzzzz_yyz_1, g_yzzzz_yz_1, g_yzzzz_yzz_1, g_yzzzz_zz_1, g_yzzzz_zzz_1, g_zzzz_xxx_0, g_zzzz_xxx_1, g_zzzz_xxz_0, g_zzzz_xxz_1, g_zzzz_xyz_0, g_zzzz_xyz_1, g_zzzz_xzz_0, g_zzzz_xzz_1, g_zzzz_yyz_0, g_zzzz_yyz_1, g_zzzz_yzz_0, g_zzzz_yzz_1, g_zzzz_zzz_0, g_zzzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzz_xxx_0[i] = g_zzzz_xxx_0[i] * fbe_0 - g_zzzz_xxx_1[i] * fz_be_0 + g_yzzzz_xxx_1[i] * pa_y[i];

        g_yyzzzz_xxy_0[i] = 3.0 * g_yyzz_xxy_0[i] * fbe_0 - 3.0 * g_yyzz_xxy_1[i] * fz_be_0 + g_yyzzz_xxy_1[i] * pa_z[i];

        g_yyzzzz_xxz_0[i] = g_zzzz_xxz_0[i] * fbe_0 - g_zzzz_xxz_1[i] * fz_be_0 + g_yzzzz_xxz_1[i] * pa_y[i];

        g_yyzzzz_xyy_0[i] = 3.0 * g_yyzz_xyy_0[i] * fbe_0 - 3.0 * g_yyzz_xyy_1[i] * fz_be_0 + g_yyzzz_xyy_1[i] * pa_z[i];

        g_yyzzzz_xyz_0[i] = g_zzzz_xyz_0[i] * fbe_0 - g_zzzz_xyz_1[i] * fz_be_0 + g_yzzzz_xz_1[i] * fe_0 + g_yzzzz_xyz_1[i] * pa_y[i];

        g_yyzzzz_xzz_0[i] = g_zzzz_xzz_0[i] * fbe_0 - g_zzzz_xzz_1[i] * fz_be_0 + g_yzzzz_xzz_1[i] * pa_y[i];

        g_yyzzzz_yyy_0[i] = 3.0 * g_yyzz_yyy_0[i] * fbe_0 - 3.0 * g_yyzz_yyy_1[i] * fz_be_0 + g_yyzzz_yyy_1[i] * pa_z[i];

        g_yyzzzz_yyz_0[i] = g_zzzz_yyz_0[i] * fbe_0 - g_zzzz_yyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_yz_1[i] * fe_0 + g_yzzzz_yyz_1[i] * pa_y[i];

        g_yyzzzz_yzz_0[i] = g_zzzz_yzz_0[i] * fbe_0 - g_zzzz_yzz_1[i] * fz_be_0 + g_yzzzz_zz_1[i] * fe_0 + g_yzzzz_yzz_1[i] * pa_y[i];

        g_yyzzzz_zzz_0[i] = g_zzzz_zzz_0[i] * fbe_0 - g_zzzz_zzz_1[i] * fz_be_0 + g_yzzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 260-270 components of targeted buffer : IF

    auto g_yzzzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 260);

    auto g_yzzzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 261);

    auto g_yzzzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 262);

    auto g_yzzzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 263);

    auto g_yzzzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 264);

    auto g_yzzzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 265);

    auto g_yzzzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 266);

    auto g_yzzzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 267);

    auto g_yzzzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 268);

    auto g_yzzzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 269);

    #pragma omp simd aligned(g_yzzzzz_xxx_0, g_yzzzzz_xxy_0, g_yzzzzz_xxz_0, g_yzzzzz_xyy_0, g_yzzzzz_xyz_0, g_yzzzzz_xzz_0, g_yzzzzz_yyy_0, g_yzzzzz_yyz_0, g_yzzzzz_yzz_0, g_yzzzzz_zzz_0, g_zzzzz_xx_1, g_zzzzz_xxx_1, g_zzzzz_xxy_1, g_zzzzz_xxz_1, g_zzzzz_xy_1, g_zzzzz_xyy_1, g_zzzzz_xyz_1, g_zzzzz_xz_1, g_zzzzz_xzz_1, g_zzzzz_yy_1, g_zzzzz_yyy_1, g_zzzzz_yyz_1, g_zzzzz_yz_1, g_zzzzz_yzz_1, g_zzzzz_zz_1, g_zzzzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzz_xxx_0[i] = g_zzzzz_xxx_1[i] * pa_y[i];

        g_yzzzzz_xxy_0[i] = g_zzzzz_xx_1[i] * fe_0 + g_zzzzz_xxy_1[i] * pa_y[i];

        g_yzzzzz_xxz_0[i] = g_zzzzz_xxz_1[i] * pa_y[i];

        g_yzzzzz_xyy_0[i] = 2.0 * g_zzzzz_xy_1[i] * fe_0 + g_zzzzz_xyy_1[i] * pa_y[i];

        g_yzzzzz_xyz_0[i] = g_zzzzz_xz_1[i] * fe_0 + g_zzzzz_xyz_1[i] * pa_y[i];

        g_yzzzzz_xzz_0[i] = g_zzzzz_xzz_1[i] * pa_y[i];

        g_yzzzzz_yyy_0[i] = 3.0 * g_zzzzz_yy_1[i] * fe_0 + g_zzzzz_yyy_1[i] * pa_y[i];

        g_yzzzzz_yyz_0[i] = 2.0 * g_zzzzz_yz_1[i] * fe_0 + g_zzzzz_yyz_1[i] * pa_y[i];

        g_yzzzzz_yzz_0[i] = g_zzzzz_zz_1[i] * fe_0 + g_zzzzz_yzz_1[i] * pa_y[i];

        g_yzzzzz_zzz_0[i] = g_zzzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 270-280 components of targeted buffer : IF

    auto g_zzzzzz_xxx_0 = pbuffer.data(idx_eri_0_if + 270);

    auto g_zzzzzz_xxy_0 = pbuffer.data(idx_eri_0_if + 271);

    auto g_zzzzzz_xxz_0 = pbuffer.data(idx_eri_0_if + 272);

    auto g_zzzzzz_xyy_0 = pbuffer.data(idx_eri_0_if + 273);

    auto g_zzzzzz_xyz_0 = pbuffer.data(idx_eri_0_if + 274);

    auto g_zzzzzz_xzz_0 = pbuffer.data(idx_eri_0_if + 275);

    auto g_zzzzzz_yyy_0 = pbuffer.data(idx_eri_0_if + 276);

    auto g_zzzzzz_yyz_0 = pbuffer.data(idx_eri_0_if + 277);

    auto g_zzzzzz_yzz_0 = pbuffer.data(idx_eri_0_if + 278);

    auto g_zzzzzz_zzz_0 = pbuffer.data(idx_eri_0_if + 279);

    #pragma omp simd aligned(g_zzzz_xxx_0, g_zzzz_xxx_1, g_zzzz_xxy_0, g_zzzz_xxy_1, g_zzzz_xxz_0, g_zzzz_xxz_1, g_zzzz_xyy_0, g_zzzz_xyy_1, g_zzzz_xyz_0, g_zzzz_xyz_1, g_zzzz_xzz_0, g_zzzz_xzz_1, g_zzzz_yyy_0, g_zzzz_yyy_1, g_zzzz_yyz_0, g_zzzz_yyz_1, g_zzzz_yzz_0, g_zzzz_yzz_1, g_zzzz_zzz_0, g_zzzz_zzz_1, g_zzzzz_xx_1, g_zzzzz_xxx_1, g_zzzzz_xxy_1, g_zzzzz_xxz_1, g_zzzzz_xy_1, g_zzzzz_xyy_1, g_zzzzz_xyz_1, g_zzzzz_xz_1, g_zzzzz_xzz_1, g_zzzzz_yy_1, g_zzzzz_yyy_1, g_zzzzz_yyz_1, g_zzzzz_yz_1, g_zzzzz_yzz_1, g_zzzzz_zz_1, g_zzzzz_zzz_1, g_zzzzzz_xxx_0, g_zzzzzz_xxy_0, g_zzzzzz_xxz_0, g_zzzzzz_xyy_0, g_zzzzzz_xyz_0, g_zzzzzz_xzz_0, g_zzzzzz_yyy_0, g_zzzzzz_yyz_0, g_zzzzzz_yzz_0, g_zzzzzz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzz_xxx_0[i] = 5.0 * g_zzzz_xxx_0[i] * fbe_0 - 5.0 * g_zzzz_xxx_1[i] * fz_be_0 + g_zzzzz_xxx_1[i] * pa_z[i];

        g_zzzzzz_xxy_0[i] = 5.0 * g_zzzz_xxy_0[i] * fbe_0 - 5.0 * g_zzzz_xxy_1[i] * fz_be_0 + g_zzzzz_xxy_1[i] * pa_z[i];

        g_zzzzzz_xxz_0[i] = 5.0 * g_zzzz_xxz_0[i] * fbe_0 - 5.0 * g_zzzz_xxz_1[i] * fz_be_0 + g_zzzzz_xx_1[i] * fe_0 + g_zzzzz_xxz_1[i] * pa_z[i];

        g_zzzzzz_xyy_0[i] = 5.0 * g_zzzz_xyy_0[i] * fbe_0 - 5.0 * g_zzzz_xyy_1[i] * fz_be_0 + g_zzzzz_xyy_1[i] * pa_z[i];

        g_zzzzzz_xyz_0[i] = 5.0 * g_zzzz_xyz_0[i] * fbe_0 - 5.0 * g_zzzz_xyz_1[i] * fz_be_0 + g_zzzzz_xy_1[i] * fe_0 + g_zzzzz_xyz_1[i] * pa_z[i];

        g_zzzzzz_xzz_0[i] = 5.0 * g_zzzz_xzz_0[i] * fbe_0 - 5.0 * g_zzzz_xzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xz_1[i] * fe_0 + g_zzzzz_xzz_1[i] * pa_z[i];

        g_zzzzzz_yyy_0[i] = 5.0 * g_zzzz_yyy_0[i] * fbe_0 - 5.0 * g_zzzz_yyy_1[i] * fz_be_0 + g_zzzzz_yyy_1[i] * pa_z[i];

        g_zzzzzz_yyz_0[i] = 5.0 * g_zzzz_yyz_0[i] * fbe_0 - 5.0 * g_zzzz_yyz_1[i] * fz_be_0 + g_zzzzz_yy_1[i] * fe_0 + g_zzzzz_yyz_1[i] * pa_z[i];

        g_zzzzzz_yzz_0[i] = 5.0 * g_zzzz_yzz_0[i] * fbe_0 - 5.0 * g_zzzz_yzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_yz_1[i] * fe_0 + g_zzzzz_yzz_1[i] * pa_z[i];

        g_zzzzzz_zzz_0[i] = 5.0 * g_zzzz_zzz_0[i] * fbe_0 - 5.0 * g_zzzz_zzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_zz_1[i] * fe_0 + g_zzzzz_zzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

